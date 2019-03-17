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
/**
 * @class GCMQ
 * @brief A generalized conforming mixed quadrilateral element.
 *
 * @author tlc
 * @date 01/03/2019
 * @version 0.2.0
 * @file GCMQ.h
 * @{
 */

#ifndef GCMQ_H
#define GCMQ_H

#ifndef ELE_GCMQ_TAG
#define ELE_GCMQ_TAG 20018
#endif

#include <Element.h>
#include <Matrix.h>
#include <NDMaterial.h>
#include <vector>

class GCMQ final : public Element {
	static Matrix mapping;

	struct IntegrationPoint final {
		Vector coor;
		double factor; // integration weight * thickness * determinant of Jacobian
		NDMaterial* m_material;
		Matrix poly_stress, poly_strain;
		IntegrationPoint(const Vector&, double, NDMaterial*);
		~IntegrationPoint();
	};

	const int mat_tag;
	const double thickness;
	const int int_scheme;

	std::vector<IntegrationPoint> int_pt;

	Matrix trans_mat;

	Matrix HT, NT, MT, N, M;

	Matrix initial_viwt, trial_viwt, current_viwt;
	Vector trial_vif, current_vif;
	Vector trial_alpha, current_alpha; // stress
	Vector trial_beta, current_beta;   // strain
	Vector trial_zeta, current_zeta;   // enhanced strain

	Vector pre_disp;

	ID node_tag;
	Node* encoding[4]{};
	Matrix ele_coor;
	Vector resistance;
	Matrix initial_stiffness, stiffness;

	static Matrix shape_disp(double, double);
	static Matrix shape_stress(const Vector&);
	static Matrix shape_strain(double, double, double);
	static Matrix form_transformation(const Matrix&);
	static Matrix form_interpolation_displacement(const Matrix&, const Matrix&);
	static Matrix form_interpolation_enhanced_strain(const Vector&);

public:
	GCMQ(int, const int*, double, int, int);

	void setDomain(Domain*) override;

	int update() override;
	int commitState() override;
	int revertToLastCommit() override;
	int revertToStart() override;

	const Matrix& getTangentStiff() override;
	const Matrix& getInitialStiff() override;
	const Vector& getResistingForce() override;

	int getNumDOF() override;
	int getNumExternalNodes() const override;
	Node** getNodePtrs() override;
	const ID& getExternalNodes() override;

	int sendSelf(int, Channel&) override;
	int recvSelf(int, Channel&, FEM_ObjectBroker&) override;

	void Print(OPS_Stream&, int) override;
};

extern "C" int dgesv_(int* N, int* NRHS, double* A, int* LDA, int* iPiv, double* B, int* LDB, int* INFO);
Matrix trans(const Matrix&);
Matrix solve(const Matrix&, const Matrix&);
Vector solve(const Matrix&, const Vector&);

#endif

//! @}
