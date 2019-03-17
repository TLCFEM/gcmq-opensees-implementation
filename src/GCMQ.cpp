/*******************************************************************************
 * Copyright (C) 2018 Theodore Chang
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/
#include "GCMQ.h"
#include "IntegrationPlan.h"
#include <Domain.h>
#include <NDMaterial.h>
#include <Node.h>
#include <elementAPI.h>

#ifdef _USRDLL
#include <windows.h>
#define OPS_Export extern "C" __declspec(dllexport)
#elif _MACOSX
#define OPS_Export extern "C" __attribute__((visibility("default")))
#else
#define OPS_Export extern "C"
#endif

OPS_Export void* OPS_GCMQ() {
	int TAG, MAT_TAG, NL = 1;
	int NODES[4];
	double THICKNESS;

	auto NUM = 1;
	OPS_GetIntInput(&NUM, &TAG);

	NUM = 4;
	OPS_GetIntInput(&NUM, NODES);

	NUM = 1;
	OPS_GetDoubleInput(&NUM, &THICKNESS);
	OPS_GetIntInput(&NUM, &MAT_TAG);

	if(OPS_GetNumRemainingInputArgs() == 1) OPS_GetIntInput(&NUM, &NL);

	return new GCMQ(TAG, NODES, THICKNESS, MAT_TAG, NL);
}

Matrix GCMQ::mapping;

GCMQ::IntegrationPoint::IntegrationPoint(const Vector& C, const double F, NDMaterial* M)
	: coor(C)
	, factor(F)
	, m_material(M)
	, poly_stress(3, 11)
	, poly_strain(3, 11) {}

GCMQ::IntegrationPoint::~IntegrationPoint() { delete m_material; }

Matrix GCMQ::shape_disp(const double X, const double Y) {
	Matrix shape_func(2, 4);

	const auto XP = 1. + X, XM = 1. - X, YP = 1. + Y, YM = 1. - Y;

	shape_func(1, 1) = -(shape_func(1, 2) = .25 * XP);
	shape_func(1, 0) = -(shape_func(1, 3) = .25 * XM);
	shape_func(0, 3) = -(shape_func(0, 2) = .25 * YP);
	shape_func(0, 0) = -(shape_func(0, 1) = .25 * YM);

	return shape_func;
}

Matrix GCMQ::shape_stress(const Vector& C) {
	Matrix N(3, 11);

	const auto& X = C(0);
	const auto& Y = C(1);

	const auto X2 = X * X, Y2 = Y * Y, XY = X * Y;

	for(auto I = 0; I < 3; ++I) N(I, I) = 1.;

	N(2, 6) = -(N(0, 4) = N(1, 5) = Y);
	N(2, 5) = -(N(1, 3) = N(0, 6) = X);

	N(1, 7) = N(0, 8) = 2. * XY;
	N(0, 9) = N(2, 7) = -X2;
	N(1, 10) = N(2, 8) = -Y2;
	N(1, 9) = 2. * X2 - Y2;
	N(0, 10) = 2. * Y2 - X2;
	N(2, 10) = N(2, 9) = 2. * XY;

	return N;
}

Matrix GCMQ::shape_strain(const double X, const double Y, const double V) {
	Matrix N(3, 11);

	const auto X2 = X * X, Y2 = Y * Y, XY = X * Y;

	N(0, 0) = N(1, 1) = 1.;

	N(2, 2) = 2. + 2. * V;

	N(0, 1) = N(1, 0) = -V;

	N(0, 3) = -V * (N(1, 3) = X);
	N(1, 4) = -V * (N(0, 4) = Y);
	N(0, 5) = N(1, 4);
	N(0, 6) = N(1, 3);

	N(1, 5) = N(0, 4);
	N(1, 6) = N(0, 3);

	N(2, 5) = -X * N(2, 2);
	N(2, 6) = -Y * N(2, 2);

	N(1, 8) = N(0, 7) = -V * (N(1, 7) = N(0, 8) = 2. * XY);

	N(2, 7) = -X2 * N(2, 2);
	N(2, 8) = -Y2 * N(2, 2);
	N(0, 9) = V * Y2 - (2. * V + 1.) * X2;
	N(1, 9) = (2. + V) * X2 - Y2;

	N(0, 10) = (2. + V) * Y2 - X2;
	N(1, 10) = V * X2 - (2. * V + 1.) * Y2;

	N(2, 10) = N(2, 9) = 2. * XY * N(2, 2);

	return N;
}

Matrix GCMQ::form_transformation(const Matrix& C) {
	const auto jacobian = shape_disp(0., 0.) * C;

	Matrix trans(3, 3);

	trans(0, 0) = jacobian(0, 0) * jacobian(0, 0);
	trans(1, 0) = jacobian(0, 1) * jacobian(0, 1);
	trans(2, 0) = jacobian(0, 0) * jacobian(0, 1);

	trans(0, 1) = jacobian(1, 0) * jacobian(1, 0);
	trans(1, 1) = jacobian(1, 1) * jacobian(1, 1);
	trans(2, 1) = jacobian(1, 0) * jacobian(1, 1);

	trans(0, 2) = 2. * jacobian(0, 0) * jacobian(1, 0);
	trans(1, 2) = 2. * jacobian(1, 0) * jacobian(1, 1);
	trans(2, 2) = jacobian(0, 0) * jacobian(1, 1) + jacobian(0, 1) * jacobian(1, 0);

	return trans;
}

Matrix GCMQ::form_interpolation_displacement(const Matrix& pn_pxy, const Matrix& pnt_pxy) {
	Matrix poly_disp(3, 12);

	for(unsigned J = 0; J < 4; ++J) {
		poly_disp(0, 3 * J) = poly_disp(2, 3 * J + 1) = pn_pxy(0, J);
		poly_disp(2, 3 * J) = poly_disp(1, 3 * J + 1) = pn_pxy(1, J);
		poly_disp(0, 3 * J + 2) = pnt_pxy(0, J);
		poly_disp(1, 3 * J + 2) = pnt_pxy(1, J + 4);
		poly_disp(2, 3 * J + 2) = pnt_pxy(0, J + 4) + pnt_pxy(1, J);
	}

	return poly_disp;
}

Matrix GCMQ::form_interpolation_enhanced_strain(const Vector& coor) {
	Matrix poly_enhanced_strain(3, 1);

	poly_enhanced_strain(0, 0) = 3. * coor(0) * coor(0) - 1.;
	poly_enhanced_strain(1, 0) = 3. * coor(1) * coor(1) - 1.;
	poly_enhanced_strain(2, 0) = 0.;

	return poly_enhanced_strain;
}

GCMQ::GCMQ(const int tag, const int* NODES, const double TH, const int MT, const int IS)
	: Element(tag, ELE_GCMQ_TAG)
	, mat_tag(MT)
	, thickness(TH)
	, int_scheme(IS)
	, trans_mat(2, 4)
	, HT(11, 11)
	, NT(11, 12)
	, MT(11, 1)
	, N(11, 12)
	, M(11, 1)
	, initial_viwt(1, 12)
	, trial_viwt(1, 12)
	, current_viwt(1, 12)
	, trial_vif(1)
	, current_vif(1)
	, trial_alpha(11)
	, current_alpha(11)
	, trial_beta(11)
	, current_beta(11)
	, trial_zeta(1)
	, current_zeta(1)
	, pre_disp(12)
	, node_tag(4)
	, ele_coor(4, 2)
	, resistance(12)
	, initial_stiffness(12, 12)
	, stiffness(12, 12) {
	for(auto I = 0; I < 4; ++I) node_tag(I) = NODES[I];
	if(mapping.noCols() == 0) {
		Matrix t_mapping(4, 4);
		for(auto J = 0; J < 4; ++J) for(auto I = 0; I < 4; ++I) t_mapping(I, J) = .25;
		t_mapping(1, 0) = t_mapping(1, 3) = t_mapping(2, 0) = t_mapping(2, 1) = t_mapping(3, 1) = t_mapping(3, 3) = -.25;
		mapping = t_mapping;
	}
}

void GCMQ::setDomain(Domain* D) {
	if(D == nullptr) return;

	for(auto I = 0; I < 4; ++I) {
		encoding[I] = D->getNode(node_tag(I));
		const auto& t_coor = encoding[I]->getCrds();
		ele_coor(I, 0) = t_coor(0), ele_coor(I, 1) = t_coor(1);
	}
	trans_mat = trans(mapping * ele_coor);

	const auto LX1 = ele_coor(1, 1) - ele_coor(0, 1);
	const auto LX2 = ele_coor(2, 1) - ele_coor(1, 1);
	const auto LX3 = ele_coor(3, 1) - ele_coor(2, 1);
	const auto LX4 = ele_coor(0, 1) - ele_coor(3, 1);
	const auto LY1 = ele_coor(0, 0) - ele_coor(1, 0);
	const auto LY2 = ele_coor(1, 0) - ele_coor(2, 0);
	const auto LY3 = ele_coor(2, 0) - ele_coor(3, 0);
	const auto LY4 = ele_coor(3, 0) - ele_coor(0, 0);

	const auto& mat_proto = OPS_GetNDMaterial(mat_tag);

	const auto& ini_stiffness = mat_proto->getInitialTangent();

	const IntegrationPlan plan(2, int_scheme == 1 ? 2 : 3, int_scheme == 1 ? IntegrationType::IRONS : int_scheme == 2 ? IntegrationType::LOBATTO : IntegrationType::GAUSS);

	const auto jacob_trans = form_transformation(ele_coor);

	int_pt.clear(), int_pt.reserve(plan.n_rows);

	Matrix pnt(2, 8), pne(2, 2), H(11, 11), HTT(11, 11);
	Vector coor(2), disp_mode(4);
	for(unsigned I = 0; I < plan.n_rows; ++I) {
		const auto& X = plan(I, 0);
		const auto& Y = plan(I, 1);

		coor(0) = X, coor(1) = Y;

		const auto pn = shape_disp(X, Y);
		const auto jacob = pn * ele_coor;
		const auto det_jacob = jacob(0, 0) * jacob(1, 1) - jacob(0, 1) * jacob(1, 0);

		int_pt.emplace_back(coor, det_jacob * thickness * plan(I, 2), mat_proto->getCopy());

		disp_mode(3) = (disp_mode(1) = X) * (disp_mode(2) = Y);
		const auto global_coor = trans_mat * disp_mode;

		int_pt[I].poly_stress = shape_stress(global_coor);
		int_pt[I].poly_strain = solve(ini_stiffness, int_pt[I].poly_stress);

		const auto X2 = 2. * X, Y2 = 2. * Y, XP = X + 1., XM = X - 1., YP = Y + 1., YM = Y - 1.;

		pnt(0, 0) = YM * (LX4 * YP - LX1 * X2);
		pnt(0, 1) = YM * (LX2 * YP + LX1 * X2);
		pnt(0, 2) = YP * (LX3 * X2 - LX2 * YM);
		pnt(0, 3) = -YP * (LX3 * X2 + LX4 * YM);
		pnt(0, 4) = YM * (LY4 * YP - LY1 * X2);
		pnt(0, 5) = YM * (LY2 * YP + LY1 * X2);
		pnt(0, 6) = YP * (LY3 * X2 - LY2 * YM);
		pnt(0, 7) = -YP * (LY3 * X2 + LY4 * YM);
		pnt(1, 0) = XM * (LX4 * Y2 - LX1 * XP);
		pnt(1, 1) = XP * (LX1 * XM + LX2 * Y2);
		pnt(1, 2) = XP * (LX3 * XM - LX2 * Y2);
		pnt(1, 3) = -XM * (LX3 * XP + LX4 * Y2);
		pnt(1, 4) = XM * (LY4 * Y2 - LY1 * XP);
		pnt(1, 5) = XP * (LY1 * XM + LY2 * Y2);
		pnt(1, 6) = XP * (LY3 * XM - LY2 * Y2);
		pnt(1, 7) = -XM * (LY3 * XP + LY4 * Y2);

		pne(0, 0) = X * X - Y * Y + X;
		pne(1, 0) = Y * Y - X * X + Y;
		pne(0, 1) = 6. * X * Y + 3. * Y * Y - 1. + X;
		pne(1, 1) = 6. * X * Y + 3. * X * X - 1. + Y;

		M.addMatrixTransposeProduct(1., int_pt[I].poly_stress, jacob_trans * form_interpolation_enhanced_strain(coor), int_pt[I].factor);
		N.addMatrixTransposeProduct(1., int_pt[I].poly_stress, form_interpolation_displacement(solve(jacob, pn), solve(jacob, pnt / 16.)), int_pt[I].factor);
		H.addMatrixTransposeProduct(1., int_pt[I].poly_stress, int_pt[I].poly_strain, int_pt[I].factor);
		HTT.addMatrixTripleProduct(1., int_pt[I].poly_strain, ini_stiffness, int_pt[I].factor);
	}

	HT = trans(H);

	NT = solve(H, N), MT = solve(H, M);

	const auto T = HTT * MT;

	const auto W = trans(NT) * T;

	initial_stiffness.addMatrixTripleProduct(0., NT, HTT, 1.);
	stiffness = initial_stiffness -= W * (trial_viwt = current_viwt = initial_viwt = solve(trans(MT) * T, trans(W)));

	DomainComponent::setDomain(D);
}

int GCMQ::update() {
	auto idx = 0;
	Vector incre_disp(12);
	for(auto& I : encoding) {
		auto& t_disp = I->getIncrDeltaDisp();
		for(unsigned pos = 0; pos < 3; ++pos) incre_disp(idx++) = t_disp(pos);
	}

	const auto incre_zeta = -1. * trial_viwt * incre_disp - trial_vif;

	trial_zeta += incre_zeta;
	trial_beta += NT * incre_disp + MT * incre_zeta;

	Vector local_stress(11);
	Matrix local_stiffness(11, 11);
	for(const auto& t_pt : int_pt) {
		if(t_pt.m_material->setTrialStrain(t_pt.poly_strain * trial_beta) != 0) return -1;
		const auto tmp_a = trans(t_pt.poly_strain) * t_pt.factor;
		local_stress += tmp_a * t_pt.m_material->getStress();
		local_stiffness += tmp_a * t_pt.m_material->getTangent() * t_pt.poly_strain;
		// local_stiffness.addMatrixTransposeProduct(1., t_pt.poly_strain, t_pt.m_material->getTangent(), t_pt.factor);
	}

	const auto T = trans(NT) * local_stiffness;

	const auto V = trans(MT) * local_stiffness * MT;

	const auto W = T * MT;

	trial_alpha = solve(HT, local_stress);
	trial_viwt = solve(V, trans(W));
	trial_vif = solve(V, trans(M) * trial_alpha);

	resistance = trans(N) * trial_alpha - W * trial_vif;

	stiffness = T * (NT - MT * trial_viwt);

	return 0;
}

int GCMQ::commitState() {
	current_zeta = trial_zeta;
	current_beta = trial_beta;
	current_alpha = trial_alpha;
	current_vif = trial_vif;
	current_viwt = trial_viwt;

	auto code = 0;
	for(const auto& I : int_pt) code += I.m_material->commitState();
	return code;
}

int GCMQ::revertToLastCommit() {
	trial_zeta = current_zeta;
	trial_beta = current_beta;
	trial_alpha = current_alpha;
	trial_vif = current_vif;
	trial_viwt = current_viwt;

	auto code = 0;
	for(const auto& I : int_pt) code += I.m_material->revertToLastCommit();
	return code;
}

int GCMQ::revertToStart() {
	current_zeta.Zero();
	trial_zeta.Zero();
	current_beta.Zero();
	trial_beta.Zero();
	current_alpha.Zero();
	trial_alpha.Zero();
	current_vif.Zero();
	trial_vif.Zero();

	current_viwt = trial_viwt = initial_viwt;

	auto code = 0;
	for(const auto& I : int_pt) code += I.m_material->revertToStart();
	return code;
}

const Matrix& GCMQ::getTangentStiff() { return stiffness; }

const Matrix& GCMQ::getInitialStiff() { return initial_stiffness; }

const Vector& GCMQ::getResistingForce() { return resistance; }

int GCMQ::sendSelf(int, Channel&) { return 0; }

int GCMQ::recvSelf(int, Channel&, FEM_ObjectBroker&) { return 0; }

void GCMQ::Print(OPS_Stream&, int) {}

int GCMQ::getNumDOF() { return 12; }

int GCMQ::getNumExternalNodes() const { return 4; }

Node** GCMQ::getNodePtrs() { return encoding; }

const ID& GCMQ::getExternalNodes() { return node_tag; }

Matrix trans(const Matrix& A) {
	Matrix B(A.noCols(), A.noRows());
	for(auto I = 0; I < A.noRows(); ++I) for(auto J = 0; J < A.noCols(); ++J) B(J, I) = A(I, J);

	return B;
}

Matrix solve(const Matrix& A, const Matrix& B) {
	auto ACOPY = A;
	auto X = B;

	auto N = ACOPY.noRows();
	auto NRHS = X.noCols();
	auto IPIV = new int[N];
	auto INFO = 0;

	dgesv_(&N, &NRHS, &ACOPY(0, 0), &N, IPIV, &X(0, 0), &N, &INFO);

	delete[] IPIV;

	return X;
}

Vector solve(const Matrix& A, const Vector& B) {
	auto ACOPY = A;
	auto X = B;

	auto N = ACOPY.noRows();
	auto NRHS = 1;
	auto IPIV = new int[N];
	auto INFO = 0;

	dgesv_(&N, &NRHS, &ACOPY(0, 0), &N, IPIV, &X(0), &N, &INFO);

	delete[] IPIV;

	return X;
}
