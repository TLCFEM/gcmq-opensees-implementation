////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2017-2019 Theodore Chang
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
////////////////////////////////////////////////////////////////////////////////

#include "IntegrationPlan.h"
#include <cmath>
#include <stdexcept>

IntegrationPlan::IntegrationPlan(const unsigned intDimension, const unsigned intOrder, const IntegrationType intType)
	: n_rows(unsigned(round(pow(intOrder, intDimension))))
	, n_cols(intDimension + 1)
	, n_elem(n_rows * n_cols) {
	if(intType == IntegrationType::IRONS) {
		const_cast<unsigned&>(n_cols) = intDimension + 1;

		if(intDimension == 3) {
			const_cast<unsigned&>(n_rows) = intOrder == 2 ? 6 : 14;
			const_cast<unsigned&>(n_elem) = n_rows * n_cols;

			int_pts = new double*[n_rows];
			for(unsigned I = 0; I < n_rows; ++I) int_pts[I] = new double[n_cols];

			if(intOrder == 2) {
				const auto WB = 4. / 3.;
				for(auto I = 0; I < 6; ++I) int_pts[I][3] = WB;
				int_pts[0][0] = int_pts[2][1] = int_pts[4][2] = -(int_pts[1][0] = int_pts[3][1] = int_pts[5][2] = 1.);
				int_pts[0][1] = int_pts[0][2] = int_pts[1][1] = int_pts[1][2] = int_pts[2][0] = int_pts[2][2] = int_pts[3][0] = int_pts[3][2] = int_pts[4][0] = int_pts[4][1] = int_pts[5][0] = int_pts[5][1] = 0.;
			} else if(intOrder == 3) {
				const auto WB = .886426593;
				const auto WC = .335180055;
				const auto LB = .795822426;
				const auto LC = .758786911;
				for(auto I = 0; I < 6; ++I) int_pts[I][3] = WB;
				for(auto I = 6; I < 14; ++I) int_pts[I][3] = WC;
				int_pts[0][0] = int_pts[2][1] = int_pts[4][2] = -(int_pts[1][0] = int_pts[3][1] = int_pts[5][2] = LB);
				int_pts[0][1] = int_pts[0][2] = int_pts[1][1] = int_pts[1][2] = int_pts[2][0] = int_pts[2][2] = int_pts[3][0] = int_pts[3][2] = int_pts[4][0] = int_pts[4][1] = int_pts[5][0] = int_pts[5][1] = 0.;
				int_pts[6][0] = int_pts[6][1] = int_pts[6][2] = int_pts[7][1] = int_pts[7][2] = int_pts[8][2] = int_pts[9][0] = int_pts[9][2] = int_pts[10][0] = int_pts[10][1] = int_pts[11][1] = int_pts[13][0] = -LC;
				int_pts[7][0] = int_pts[8][0] = int_pts[8][1] = int_pts[9][1] = int_pts[10][2] = int_pts[11][0] = int_pts[11][2] = int_pts[12][0] = int_pts[12][1] = int_pts[12][2] = int_pts[13][1] = int_pts[13][2] = LC;
			}

			return;
		}

		if(intDimension == 2) {
			const_cast<unsigned&>(n_rows) = 5;
			const_cast<unsigned&>(n_elem) = n_rows * n_cols;

			int_pts = new double*[n_rows];
			for(unsigned I = 0; I < n_rows; ++I) int_pts[I] = new double[n_cols];

			if(intOrder == 2) {
				const auto WB = 2. / 3.;
				for(auto I = 0; I < 5; ++I) {
					for(auto J = 0; J < 2; ++J) int_pts[I][J] = 0.;
					int_pts[I][2] = WB;
				}
				int_pts[0][2] *= 2.;
				int_pts[1][0] = int_pts[2][1] = -1.;
				int_pts[3][0] = int_pts[4][1] = 1.;
			}

			return;
		}
	}

	const auto PTL = new double[intOrder];
	const auto PTW = new double[intOrder];

	// GAUSS INTEGRATION
	if(intType == IntegrationType::GAUSS) {
		switch(intOrder) {
		case 1: {
			PTL[0] = 0.; // POINTS LOCATION.
			PTW[0] = 2.; // POINTS WEIGHT.
			break;
		}
		case 2: {
			const auto TMPA = 1. / sqrt(3.);
			PTL[0] = -TMPA;
			PTL[1] = TMPA;

			PTW[0] = 1.;
			PTW[1] = 1.;
			break;
		}
		case 3: {
			const auto TMPA = sqrt(.6);
			PTL[0] = -TMPA;
			PTL[1] = 0.;
			PTL[2] = TMPA;

			PTW[0] = 5. / 9.;
			PTW[1] = 8. / 9.;
			PTW[2] = 5. / 9.;
			break;
		}
		case 4: {
			const auto TMPA = 3. / 7.;
			const auto TMPB = 2. / 7. * sqrt(1.2);
			const auto TMPC = sqrt(TMPA - TMPB);
			const auto TMPD = sqrt(TMPA + TMPB);
			const auto TMPE = sqrt(30.) / 36.;
			PTL[0] = -TMPD;
			PTL[1] = -TMPC;
			PTL[2] = TMPC;
			PTL[3] = TMPD;

			PTW[0] = -TMPE + .5;
			PTW[1] = TMPE + .5;
			PTW[2] = TMPE + .5;
			PTW[3] = -TMPE + .5;
			break;
		}
		case 5: {
			const auto TMPA = 2. * sqrt(10. / 7.);
			const auto TMPB = sqrt(5. - TMPA);
			const auto TMPC = sqrt(5. + TMPA);
			const auto TMPD = 13. * sqrt(70.);
			PTL[0] = -TMPC / 3.;
			PTL[1] = -TMPB / 3.;
			PTL[2] = 0.;
			PTL[3] = TMPB / 3.;
			PTL[4] = TMPC / 3.;

			PTW[0] = (322. - TMPD) / 900.;
			PTW[1] = (322. + TMPD) / 900.;
			PTW[2] = 512. / 900.;
			PTW[3] = (322. + TMPD) / 900.;
			PTW[4] = (322. - TMPD) / 900.;
			break;
		}
		case 6: {
			PTL[0] = -.932469514203152028;
			PTL[1] = -.661209386466264514;
			PTL[2] = -.238619186083196909;
			PTL[3] = .238619186083196909;
			PTL[4] = .66120938646626451;
			PTL[5] = .932469514203152028;

			PTW[0] = .171324492379170345;
			PTW[1] = .360761573048138608;
			PTW[2] = .46791393457269105;
			PTW[3] = .46791393457269105;
			PTW[4] = .360761573048138608;
			PTW[5] = .171324492379170345;
			break;
		}
		case 7: {
			PTL[0] = -.949107912342758525;
			PTL[1] = -.74153118559939444;
			PTL[2] = -.405845151377397167;
			PTL[3] = .0;
			PTL[4] = .405845151377397167;
			PTL[5] = .74153118559939444;
			PTL[6] = .949107912342758525;

			PTW[0] = .129484966168869693;
			PTW[1] = .279705391489276668;
			PTW[2] = .381830050505118945;
			PTW[3] = .417959183673469388;
			PTW[4] = .38183005050511895;
			PTW[5] = .279705391489276668;
			PTW[6] = .129484966168869693;
			break;
		}
		case 8: {
			PTL[0] = -.960289856497536232;
			PTL[1] = -.79666647741362674;
			PTL[2] = -.525532409916328986;
			PTL[3] = -.183434642495649805;
			PTL[4] = .183434642495649805;
			PTL[5] = .525532409916328986;
			PTL[6] = .79666647741362674;
			PTL[7] = .960289856497536232;

			PTW[0] = .101228536290376259;
			PTW[1] = .22238103445337447;
			PTW[2] = .313706645877887287;
			PTW[3] = .362683783378361983;
			PTW[4] = .36268378337836198;
			PTW[5] = .313706645877887287;
			PTW[6] = .222381034453374471;
			PTW[7] = .10122853629037626;
			break;
		}
		case 9: {
			PTL[0] = -.96816023950762609;
			PTL[1] = -.836031107326635794;
			PTL[2] = -.613371432700590397;
			PTL[3] = -.32425342340380893;
			PTL[4] = .0;
			PTL[5] = .32425342340380893;
			PTL[6] = .613371432700590397;
			PTL[7] = .836031107326635794;
			PTL[8] = .96816023950762609;

			PTW[0] = .081274388361574412;
			PTW[1] = .1806481606948574;
			PTW[2] = .260610696402935462;
			PTW[3] = .31234707704000284;
			PTW[4] = .33023935500125976;
			PTW[5] = .31234707704000284;
			PTW[6] = .26061069640293546;
			PTW[7] = .180648160694857404;
			PTW[8] = .081274388361574412;
			break;
		}
		case 10: {
			PTL[0] = -.97390652851717172;
			PTL[1] = -.86506336668898451;
			PTL[2] = -.679409568299024406;
			PTL[3] = -.433395394129247191;
			PTL[4] = -.14887433898163121;
			PTL[5] = .148874338981631211;
			PTL[6] = .433395394129247191;
			PTL[7] = .679409568299024406;
			PTL[8] = .865063366688984511;
			PTL[9] = .97390652851717172;

			PTW[0] = .066671344308688138;
			PTW[1] = .14945134915058059;
			PTW[2] = .219086362515982044;
			PTW[3] = .26926671930999636;
			PTW[4] = .29552422471475287;
			PTW[5] = .29552422471475287;
			PTW[6] = .269266719309996355;
			PTW[7] = .21908636251598204;
			PTW[8] = .14945134915058059;
			PTW[9] = .066671344308688138;
			break;
		}
		case 11: {
			PTL[0] = -0.978228658146056993;
			PTL[1] = -0.887062599768095299;
			PTL[2] = -0.730152005574049324;
			PTL[3] = -0.519096129206811816;
			PTL[4] = -0.26954315595234497;
			PTL[5] = 0.;
			PTL[6] = 0.269543155952344972;
			PTL[7] = 0.519096129206811816;
			PTL[8] = 0.730152005574049324;
			PTL[9] = 0.887062599768095299;
			PTL[10] = 0.978228658146056993;

			PTW[0] = 0.055668567116173666;
			PTW[1] = 0.125580369464904625;
			PTW[2] = 0.186290210927734251;
			PTW[3] = 0.23319376459199048;
			PTW[4] = 0.26280454451024666;
			PTW[5] = 0.27292508677790063;
			PTW[6] = 0.262804544510246662;
			PTW[7] = 0.23319376459199048;
			PTW[8] = 0.18629021092773425;
			PTW[9] = 0.12558036946490462;
			PTW[10] = 0.055668567116173666;
			break;
		}
		case 12: {
			PTL[0] = -0.981560634246719251;
			PTL[1] = -0.904117256370474857;
			PTL[2] = -0.769902674194304687;
			PTL[3] = -0.587317954286617447;
			PTL[4] = -0.367831498998180194;
			PTL[5] = -0.125233408511468916;
			PTL[6] = 0.125233408511468916;
			PTL[7] = 0.367831498998180194;
			PTL[8] = 0.587317954286617447;
			PTL[9] = 0.769902674194304687;
			PTL[10] = 0.904117256370474857;
			PTL[11] = 0.981560634246719251;

			PTW[0] = 0.047175336386511827;
			PTW[1] = 0.10693932599531843;
			PTW[2] = 0.160078328543346226;
			PTW[3] = 0.203167426723065922;
			PTW[4] = 0.23349253653835481;
			PTW[5] = 0.249147045813402785;
			PTW[6] = 0.24914704581340279;
			PTW[7] = 0.233492536538354809;
			PTW[8] = 0.203167426723065922;
			PTW[9] = 0.16007832854334623;
			PTW[10] = 0.106939325995318431;
			PTW[11] = 0.047175336386511827;
			break;
		}
		case 13: {
			PTL[0] = -0.98418305471858815;
			PTL[1] = -0.917598399222977965;
			PTL[2] = -0.801578090733309913;
			PTL[3] = -0.642349339440340221;
			PTL[4] = -0.448492751036446853;
			PTL[5] = -0.230458315955134794;
			PTL[6] = 0;
			PTL[7] = 0.23045831595513479;
			PTL[8] = 0.44849275103644685;
			PTL[9] = 0.642349339440340221;
			PTL[10] = 0.801578090733309913;
			PTL[11] = 0.917598399222977965;
			PTL[12] = 0.98418305471858815;

			PTW[0] = 0.04048400476531588;
			PTW[1] = 0.092121499837728448;
			PTW[2] = 0.138873510219787239;
			PTW[3] = 0.178145980761945738;
			PTW[4] = 0.207816047536888502;
			PTW[5] = 0.22628318026289724;
			PTW[6] = 0.23255155323087391;
			PTW[7] = 0.22628318026289724;
			PTW[8] = 0.207816047536888502;
			PTW[9] = 0.17814598076194574;
			PTW[10] = 0.138873510219787239;
			PTW[11] = 0.09212149983772845;
			PTW[12] = 0.04048400476531588;
			break;
		}
		case 14: {
			PTL[0] = -0.986283808696812339;
			PTL[1] = -0.928434883663573517;
			PTL[2] = -0.827201315069764993;
			PTL[3] = -0.68729290481168547;
			PTL[4] = -0.515248636358154092;
			PTL[5] = -0.31911236892788976;
			PTL[6] = -0.108054948707343662;
			PTL[7] = 0.108054948707343662;
			PTL[8] = 0.31911236892788976;
			PTL[9] = 0.515248636358154092;
			PTL[10] = 0.68729290481168547;
			PTL[11] = 0.827201315069764993;
			PTL[12] = 0.928434883663573517;
			PTL[13] = 0.986283808696812339;

			PTW[0] = 0.035119460331751863;
			PTW[1] = 0.08015808715976021;
			PTW[2] = 0.121518570687903185;
			PTW[3] = 0.157203167158193535;
			PTW[4] = 0.185538397477937814;
			PTW[5] = 0.205198463721295604;
			PTW[6] = 0.21526385346315779;
			PTW[7] = 0.21526385346315779;
			PTW[8] = 0.205198463721295604;
			PTW[9] = 0.18553839747793781;
			PTW[10] = 0.157203167158193535;
			PTW[11] = 0.121518570687903185;
			PTW[12] = 0.08015808715976021;
			PTW[13] = 0.035119460331751863;
			break;
		}
		case 15: {
			PTL[0] = -0.987992518020485429;
			PTL[1] = -0.9372733924007059;
			PTL[2] = -0.848206583410427216;
			PTL[3] = -0.724417731360170047;
			PTL[4] = -0.570972172608538848;
			PTL[5] = -0.39415134707756337;
			PTL[6] = -0.201194093997434522;
			PTL[7] = 0.;
			PTL[8] = 0.201194093997434522;
			PTL[9] = 0.39415134707756337;
			PTL[10] = 0.570972172608538848;
			PTL[11] = 0.724417731360170047;
			PTL[12] = 0.848206583410427216;
			PTL[13] = 0.937273392400705904;
			PTL[14] = 0.987992518020485429;

			PTW[0] = 0.0307532419961172684;
			PTW[1] = 0.0703660474881081247;
			PTW[2] = 0.107159220467171935;
			PTW[3] = 0.139570677926154314;
			PTW[4] = 0.166269205816993934;
			PTW[5] = 0.18616100001556221;
			PTW[6] = 0.198431485327111577;
			PTW[7] = 0.202578241925561273;
			PTW[8] = 0.19843148532711158;
			PTW[9] = 0.18616100001556221;
			PTW[10] = 0.166269205816993934;
			PTW[11] = 0.13957067792615431;
			PTW[12] = 0.10715922046717194;
			PTW[13] = 0.07036604748810812;
			PTW[14] = 0.030753241996117268;
			break;
		}
		case 16: {
			PTL[0] = -0.989400934991649933;
			PTL[1] = -0.944575023073232576;
			PTL[2] = -0.865631202387831744;
			PTL[3] = -0.75540440835500303;
			PTL[4] = -0.61787624440264375;
			PTL[5] = -0.458016777657227386;
			PTL[6] = -0.281603550779258913;
			PTL[7] = -0.09501250983763744;
			PTL[8] = 0.09501250983763744;
			PTL[9] = 0.281603550779258913;
			PTL[10] = 0.458016777657227386;
			PTL[11] = 0.617876244402643748;
			PTL[12] = 0.755404408355003034;
			PTL[13] = 0.86563120238783174;
			PTL[14] = 0.944575023073232576;
			PTL[15] = 0.989400934991649933;

			PTW[0] = 0.027152459411754095;
			PTW[1] = 0.062253523938647893;
			PTW[2] = 0.095158511682492785;
			PTW[3] = 0.12462897125553387;
			PTW[4] = 0.149595988816576732;
			PTW[5] = 0.169156519395002538;
			PTW[6] = 0.18260341504492359;
			PTW[7] = 0.1894506104550685;
			PTW[8] = 0.189450610455068496;
			PTW[9] = 0.182603415044923589;
			PTW[10] = 0.169156519395002538;
			PTW[11] = 0.149595988816576732;
			PTW[12] = 0.12462897125553387;
			PTW[13] = 0.09515851168249278;
			PTW[14] = 0.062253523938647893;
			PTW[15] = 0.02715245941175409;
			break;
		}
		case 17: {
			PTL[0] = -0.990575475314417336;
			PTL[1] = -0.950675521768767761;
			PTL[2] = -0.880239153726985902;
			PTL[3] = -0.78151400389680141;
			PTL[4] = -0.657671159216690766;
			PTL[5] = -0.512690537086476968;
			PTL[6] = -0.351231763453876315;
			PTL[7] = -0.178484181495847856;
			PTL[8] = 0.;
			PTL[9] = 0.178484181495847856;
			PTL[10] = 0.351231763453876315;
			PTL[11] = 0.512690537086476968;
			PTL[12] = 0.657671159216690766;
			PTL[13] = 0.781514003896801407;
			PTL[14] = 0.880239153726985902;
			PTL[15] = 0.950675521768767761;
			PTL[16] = 0.99057547531441734;

			PTW[0] = 0.024148302868547932;
			PTW[1] = 0.055459529373987201;
			PTW[2] = 0.085036148317179181;
			PTW[3] = 0.11188384719340397;
			PTW[4] = 0.135136368468525473;
			PTW[5] = 0.15404576107681029;
			PTW[6] = 0.168004102156450045;
			PTW[7] = 0.176562705366992646;
			PTW[8] = 0.179446470356206526;
			PTW[9] = 0.17656270536699265;
			PTW[10] = 0.168004102156450045;
			PTW[11] = 0.154045761076810288;
			PTW[12] = 0.13513636846852547;
			PTW[13] = 0.111883847193403971;
			PTW[14] = 0.085036148317179181;
			PTW[15] = 0.055459529373987201;
			PTW[16] = 0.024148302868547932;
			break;
		}
		case 18: {
			PTL[0] = -0.991565168420930947;
			PTL[1] = -0.955823949571397755;
			PTL[2] = -0.892602466497555739;
			PTL[3] = -0.803704958972523116;
			PTL[4] = -0.69168704306035321;
			PTL[5] = -0.559770831073947535;
			PTL[6] = -0.411751161462842646;
			PTL[7] = -0.25188622569150551;
			PTL[8] = -0.084775013041735301;
			PTL[9] = 0.084775013041735301;
			PTL[10] = 0.25188622569150551;
			PTL[11] = 0.411751161462842646;
			PTL[12] = 0.559770831073947535;
			PTL[13] = 0.691687043060353208;
			PTL[14] = 0.803704958972523116;
			PTL[15] = 0.892602466497555739;
			PTL[16] = 0.955823949571397755;
			PTL[17] = 0.991565168420930947;

			PTW[0] = 0.02161601352648331;
			PTW[1] = 0.0497145488949698;
			PTW[2] = 0.076425730254889057;
			PTW[3] = 0.100942044106287166;
			PTW[4] = 0.12255520671147846;
			PTW[5] = 0.140642914670650651;
			PTW[6] = 0.154684675126265245;
			PTW[7] = 0.164276483745832723;
			PTW[8] = 0.169142382963143592;
			PTW[9] = 0.169142382963143592;
			PTW[10] = 0.16427648374583272;
			PTW[11] = 0.15468467512626524;
			PTW[12] = 0.140642914670650651;
			PTW[13] = 0.12255520671147846;
			PTW[14] = 0.100942044106287166;
			PTW[15] = 0.07642573025488906;
			PTW[16] = 0.049714548894969796;
			PTW[17] = 0.02161601352648331;
			break;
		}
		case 19: {
			PTL[0] = -0.992406843843584403;
			PTL[1] = -0.960208152134830031;
			PTL[2] = -0.903155903614817902;
			PTL[3] = -0.82271465653714283;
			PTL[4] = -0.720966177335229379;
			PTL[5] = -0.600545304661681024;
			PTL[6] = -0.464570741375960946;
			PTL[7] = -0.316564099963629832;
			PTL[8] = -0.160358645640225376;
			PTL[9] = 0;
			PTL[10] = 0.160358645640225376;
			PTL[11] = 0.316564099963629832;
			PTL[12] = 0.464570741375960946;
			PTL[13] = 0.600545304661681024;
			PTL[14] = 0.720966177335229379;
			PTL[15] = 0.822714656537142825;
			PTL[16] = 0.903155903614817902;
			PTL[17] = 0.960208152134830031;
			PTL[18] = 0.992406843843584403;

			PTW[0] = 0.019461788229726477;
			PTW[1] = 0.0448142267656996;
			PTW[2] = 0.069044542737641227;
			PTW[3] = 0.091490021622449999;
			PTW[4] = 0.111566645547333995;
			PTW[5] = 0.128753962539336228;
			PTW[6] = 0.142606702173606612;
			PTW[7] = 0.152766042065859667;
			PTW[8] = 0.15896884339395435;
			PTW[9] = 0.161054449848783696;
			PTW[10] = 0.158968843393954348;
			PTW[11] = 0.152766042065859667;
			PTW[12] = 0.14260670217360661;
			PTW[13] = 0.128753962539336228;
			PTW[14] = 0.111566645547333995;
			PTW[15] = 0.091490021622449999;
			PTW[16] = 0.06904454273764123;
			PTW[17] = 0.0448142267656996;
			PTW[18] = 0.019461788229726477;
			break;
		}
		case 20:
		default: {
			PTL[0] = -0.993128599185094925;
			PTL[1] = -0.963971927277913791;
			PTL[2] = -0.912234428251325906;
			PTL[3] = -0.839116971822218823;
			PTL[4] = -0.746331906460150793;
			PTL[5] = -0.636053680726515026;
			PTL[6] = -0.510867001950827098;
			PTL[7] = -0.373706088715419561;
			PTL[8] = -0.227785851141645078;
			PTL[9] = -0.0765265211334973338;
			PTL[10] = 0.076526521133497334;
			PTL[11] = 0.227785851141645078;
			PTL[12] = 0.373706088715419561;
			PTL[13] = 0.510867001950827098;
			PTL[14] = 0.636053680726515026;
			PTL[15] = 0.746331906460150793;
			PTL[16] = 0.839116971822218823;
			PTL[17] = 0.912234428251325906;
			PTL[18] = 0.963971927277913791;
			PTL[19] = 0.993128599185094925;

			PTW[0] = 0.017614007139152118;
			PTW[1] = 0.040601429800386941;
			PTW[2] = 0.062672048334109064;
			PTW[3] = 0.083276741576704749;
			PTW[4] = 0.10193011981724044;
			PTW[5] = 0.118194531961518417;
			PTW[6] = 0.131688638449176627;
			PTW[7] = 0.142096109318382051;
			PTW[8] = 0.149172986472603747;
			PTW[9] = 0.152753387130725851;
			PTW[10] = 0.15275338713072585;
			PTW[11] = 0.149172986472603747;
			PTW[12] = 0.142096109318382051;
			PTW[13] = 0.13168863844917663;
			PTW[14] = 0.118194531961518417;
			PTW[15] = 0.10193011981724044;
			PTW[16] = 0.083276741576704749;
			PTW[17] = 0.062672048334109064;
			PTW[18] = 0.040601429800386941;
			PTW[19] = 0.0176140071391521183;
			break;
		}
		}
	}
		// LOBATTO INTEGRATION
	else if(intType == IntegrationType::LOBATTO) {
		switch(intOrder) {
		case 3: {
			PTL[0] = -1.;
			PTL[1] = 0.;
			PTL[2] = 1.;
			PTW[0] = 1. / 3.;
			PTW[1] = 4. / 3.;
			PTW[2] = 1. / 3.;
			break;
		}
		case 4: {
			const auto TMPA = sqrt(.2);
			PTL[0] = -1.;
			PTL[1] = -TMPA;
			PTL[2] = TMPA;
			PTL[3] = 1.;
			PTW[0] = 1. / 6.;
			PTW[1] = 5. / 6.;
			PTW[2] = 5. / 6.;
			PTW[3] = 1. / 6.;
			break;
		}
		case 5: {
			PTL[0] = -1.;
			PTL[1] = -sqrt(3. / 7.);
			PTL[2] = 0.;
			PTL[3] = sqrt(3. / 7.);
			PTL[4] = 1.;
			PTW[0] = .1;
			PTW[1] = 49. / 90.;
			PTW[2] = 64. / 90.;
			PTW[3] = 49. / 90.;
			PTW[4] = .1;
			break;
		}
		case 6: {
			const auto TMPA = 2. * sqrt(7.);
			const auto TMPB = sqrt((7. - TMPA) / 21.);
			const auto TMPC = sqrt((7. + TMPA) / 21.);
			const auto TMPD = 14. / 30.;
			const auto TMPE = sqrt(7.) / 30.;
			PTL[0] = -1.;
			PTL[1] = -TMPC;
			PTL[2] = -TMPB;
			PTL[3] = TMPB;
			PTL[4] = TMPC;
			PTL[5] = 1.;
			PTW[0] = 1. / 15.;
			PTW[1] = TMPD - TMPE;
			PTW[2] = TMPD + TMPE;
			PTW[3] = TMPD + TMPE;
			PTW[4] = TMPD - TMPE;
			PTW[5] = 1. / 15.;
			break;
		}
		case 7: {
			PTL[0] = -1.;
			PTL[1] = -0.83022389627856693;
			PTL[2] = -0.468848793470714214;
			PTL[3] = 0.;
			PTL[4] = 0.468848793470714214;
			PTL[5] = 0.83022389627856693;
			PTL[6] = 1.;

			PTW[0] = 0.0476190476190476191;
			PTW[1] = 0.276826047361565948;
			PTW[2] = 0.431745381209862623;
			PTW[3] = 0.48761904761904762;
			PTW[4] = 0.431745381209862623;
			PTW[5] = 0.276826047361565948;
			PTW[6] = 0.0476190476190476191;

			break;
		}
		case 8: {
			PTL[0] = -1.;
			PTL[1] = -0.871740148509606615;
			PTL[2] = -0.591700181433142302;
			PTL[3] = -0.209299217902478869;
			PTL[4] = 0.209299217902478869;
			PTL[5] = 0.591700181433142302;
			PTL[6] = 0.871740148509606615;
			PTL[7] = 1.;

			PTW[0] = 0.0357142857142857143;
			PTW[1] = 0.210704227143506039;
			PTW[2] = 0.34112269248350436;
			PTW[3] = 0.412458794658703882;
			PTW[4] = 0.412458794658703882;
			PTW[5] = 0.341122692483504365;
			PTW[6] = 0.21070422714350604;
			PTW[7] = 0.0357142857142857143;

			break;
		}
		case 9: {
			PTL[0] = -1.;
			PTL[1] = -0.899757995411460157;
			PTL[2] = -0.677186279510737753;
			PTL[3] = -0.363117463826178159;
			PTL[4] = 0.;
			PTL[5] = 0.363117463826178159;
			PTL[6] = 0.677186279510737753;
			PTL[7] = 0.899757995411460157;
			PTL[8] = 1.;

			PTW[0] = 0.0277777777777777778;
			PTW[1] = 0.165495361560805525;
			PTW[2] = 0.274538712500161735;
			PTW[3] = 0.346428510973046345;
			PTW[4] = 0.371519274376417234;
			PTW[5] = 0.34642851097304635;
			PTW[6] = 0.274538712500161735;
			PTW[7] = 0.16549536156080553;
			PTW[8] = 0.0277777777777777778;

			break;
		}
		case 10: {
			PTL[0] = -1.;
			PTL[1] = -0.919533908166458814;
			PTL[2] = -0.738773865105505075;
			PTL[3] = -0.477924949810444496;
			PTL[4] = -0.165278957666387025;
			PTL[5] = 0.165278957666387025;
			PTL[6] = 0.477924949810444496;
			PTL[7] = 0.738773865105505075;
			PTL[8] = 0.919533908166458814;
			PTL[9] = 1.;

			PTW[0] = 0.0222222222222222222;
			PTW[1] = 0.13330599085107011;
			PTW[2] = 0.22488934206312645;
			PTW[3] = 0.292042683679683758;
			PTW[4] = 0.327539761183897457;
			PTW[5] = 0.327539761183897457;
			PTW[6] = 0.292042683679683758;
			PTW[7] = 0.22488934206312645;
			PTW[8] = 0.133305990851070111;
			PTW[9] = 0.0222222222222222222;

			break;
		}
		case 11: {
			PTL[0] = -1.;
			PTL[1] = -0.934001430408059134;
			PTL[2] = -0.78448347366314442;
			PTL[3] = -0.565235326996205007;
			PTL[4] = -0.295758135586939391;
			PTL[5] = 0.;
			PTL[6] = 0.295758135586939391;
			PTL[7] = 0.565235326996205007;
			PTL[8] = 0.784483473663144419;
			PTL[9] = 0.934001430408059134;
			PTL[10] = 1.;

			PTW[0] = 0.0181818181818181818;
			PTW[1] = 0.109612273266994865;
			PTW[2] = 0.187169881780305204;
			PTW[3] = 0.248048104264028314;
			PTW[4] = 0.286879124779008089;
			PTW[5] = 0.30021759545569069;
			PTW[6] = 0.286879124779008089;
			PTW[7] = 0.248048104264028314;
			PTW[8] = 0.187169881780305204;
			PTW[9] = 0.109612273266994865;
			PTW[10] = 0.0181818181818181818;

			break;
		}
		case 12: {
			PTL[0] = -1.;
			PTL[1] = -0.944899272222882223;
			PTL[2] = -0.819279321644006678;
			PTL[3] = -0.632876153031860678;
			PTL[4] = -0.399530940965348932;
			PTL[5] = -0.136552932854927555;
			PTL[6] = 0.136552932854927555;
			PTL[7] = 0.399530940965348932;
			PTL[8] = 0.632876153031860678;
			PTL[9] = 0.819279321644006678;
			PTL[10] = 0.944899272222882223;
			PTL[11] = 1.;

			PTW[0] = 0.0151515151515151515;
			PTW[1] = 0.09168451741319613;
			PTW[2] = 0.15797470556437012;
			PTW[3] = 0.212508417761021145;
			PTW[4] = 0.25127560319920128;
			PTW[5] = 0.271405240910696177;
			PTW[6] = 0.271405240910696177;
			PTW[7] = 0.25127560319920128;
			PTW[8] = 0.21250841776102115;
			PTW[9] = 0.157974705564370115;
			PTW[10] = 0.091684517413196131;
			PTW[11] = 0.0151515151515151515;

			break;
		}
		case 13: {
			PTL[0] = -1.;
			PTL[1] = -0.953309846642163912;
			PTL[2] = -0.84634756465187232;
			PTL[3] = -0.686188469081757426;
			PTL[4] = -0.482909821091336202;
			PTL[5] = -0.249286930106239993;
			PTL[6] = 0.;
			PTL[7] = 0.249286930106239993;
			PTL[8] = 0.482909821091336202;
			PTL[9] = 0.686188469081757426;
			PTL[10] = 0.84634756465187232;
			PTL[11] = 0.953309846642163912;
			PTL[12] = 1.;

			PTW[0] = 0.0128205128205128205;
			PTW[1] = 0.07780168674681893;
			PTW[2] = 0.134981926689608349;
			PTW[3] = 0.18364686520355009;
			PTW[4] = 0.220767793566110086;
			PTW[5] = 0.24401579030667636;
			PTW[6] = 0.251930849333446736;
			PTW[7] = 0.244015790306676357;
			PTW[8] = 0.220767793566110086;
			PTW[9] = 0.183646865203550092;
			PTW[10] = 0.13498192668960835;
			PTW[11] = 0.07780168674681893;
			PTW[12] = 0.0128205128205128205;

			break;
		}
		case 14: {
			PTL[0] = -1.;
			PTL[1] = -0.959935045267260901;
			PTL[2] = -0.867801053830347251;
			PTL[3] = -0.728868599091326141;
			PTL[4] = -0.550639402928647055;
			PTL[5] = -0.342724013342712845;
			PTL[6] = -0.116331868883703868;
			PTL[7] = 0.116331868883703868;
			PTL[8] = 0.342724013342712845;
			PTL[9] = 0.550639402928647055;
			PTL[10] = 0.728868599091326141;
			PTL[11] = 0.867801053830347251;
			PTL[12] = 0.959935045267260901;
			PTL[13] = 1.;

			PTW[0] = 0.010989010989010989;
			PTW[1] = 0.06683728449768128;
			PTW[2] = 0.116586655898711652;
			PTW[3] = 0.160021851762952142;
			PTW[4] = 0.19482614937341612;
			PTW[5] = 0.219126253009770755;
			PTW[6] = 0.231612794468457059;
			PTW[7] = 0.23161279446845706;
			PTW[8] = 0.219126253009770755;
			PTW[9] = 0.19482614937341612;
			PTW[10] = 0.16002185176295214;
			PTW[11] = 0.116586655898711652;
			PTW[12] = 0.066837284497681285;
			PTW[13] = 0.010989010989010989;

			break;
		}
		case 15: {
			PTL[0] = -1.;
			PTL[1] = -0.965245926503838573;
			PTL[2] = -0.885082044222976299;
			PTL[3] = -0.763519689951815201;
			PTL[4] = -0.606253205469845711;
			PTL[5] = -0.420638054713672481;
			PTL[6] = -0.215353955363794238;
			PTL[7] = 0.;
			PTL[8] = 0.215353955363794238;
			PTL[9] = 0.420638054713672481;
			PTL[10] = 0.606253205469845711;
			PTL[11] = 0.763519689951815201;
			PTL[12] = 0.885082044222976299;
			PTL[13] = 0.965245926503838573;
			PTL[14] = 1.;

			PTW[0] = 0.00952380952380952381;
			PTW[1] = 0.05802989302860125;
			PTW[2] = 0.101660070325718068;
			PTW[3] = 0.14051169980242811;
			PTW[4] = 0.172789647253600949;
			PTW[5] = 0.19698723596461336;
			PTW[6] = 0.21197358592682092;
			PTW[7] = 0.21704811634881565;
			PTW[8] = 0.21197358592682092;
			PTW[9] = 0.19698723596461336;
			PTW[10] = 0.17278964725360095;
			PTW[11] = 0.14051169980242811;
			PTW[12] = 0.101660070325718068;
			PTW[13] = 0.058029893028601249;
			PTW[14] = 0.00952380952380952381;

			break;
		}
		case 16: {
			PTL[0] = -1.;
			PTL[1] = -0.969568046270217933;
			PTL[2] = -0.899200533093472093;
			PTL[3] = -0.792008291861815064;
			PTL[4] = -0.65238870288249309;
			PTL[5] = -0.486059421887137612;
			PTL[6] = -0.299830468900763208;
			PTL[7] = -0.101326273521949448;
			PTL[8] = 0.101326273521949448;
			PTL[9] = 0.299830468900763208;
			PTL[10] = 0.486059421887137612;
			PTL[11] = 0.65238870288249309;
			PTL[12] = 0.79200829186181506;
			PTL[13] = 0.899200533093472093;
			PTL[14] = 0.969568046270217933;
			PTL[15] = 1.;

			PTW[0] = 0.00833333333333333333;
			PTW[1] = 0.05085036100591991;
			PTW[2] = 0.089393697325930801;
			PTW[3] = 0.124255382132514098;
			PTW[4] = 0.154026980807164281;
			PTW[5] = 0.177491913391704125;
			PTW[6] = 0.19369002382520358;
			PTW[7] = 0.201958308178229872;
			PTW[8] = 0.201958308178229872;
			PTW[9] = 0.19369002382520358;
			PTW[10] = 0.177491913391704125;
			PTW[11] = 0.154026980807164281;
			PTW[12] = 0.124255382132514098;
			PTW[13] = 0.089393697325930801;
			PTW[14] = 0.050850361005919905;
			PTW[15] = 0.00833333333333333333;

			break;
		}
		case 17: {
			PTL[0] = -1.;
			PTL[1] = -0.973132176631418314;
			PTL[2] = -0.910879995915573596;
			PTL[3] = -0.815696251221770307;
			PTL[4] = -0.691028980627684705;
			PTL[5] = -0.541385399330101539;
			PTL[6] = -0.372174433565477042;
			PTL[7] = -0.189511973518317388;
			PTL[8] = 0.;
			PTL[9] = 0.189511973518317388;
			PTL[10] = 0.372174433565477042;
			PTL[11] = 0.541385399330101539;
			PTL[12] = 0.691028980627684705;
			PTL[13] = 0.815696251221770307;
			PTL[14] = 0.910879995915573596;
			PTL[15] = 0.973132176631418314;
			PTL[16] = 1.;

			PTW[0] = 0.00735294117647058824;
			PTW[1] = 0.04492194054325421;
			PTW[2] = 0.079198270503687119;
			PTW[3] = 0.110592909007028161;
			PTW[4] = 0.137987746201926559;
			PTW[5] = 0.16039466199762154;
			PTW[6] = 0.17700425351565787;
			PTW[7] = 0.18721633967761924;
			PTW[8] = 0.190661874753469433;
			PTW[9] = 0.187216339677619236;
			PTW[10] = 0.17700425351565787;
			PTW[11] = 0.16039466199762154;
			PTW[12] = 0.137987746201926559;
			PTW[13] = 0.110592909007028161;
			PTW[14] = 0.079198270503687119;
			PTW[15] = 0.04492194054325421;
			PTW[16] = 0.00735294117647058824;

			break;
		}
		case 18: {
			PTL[0] = -1.;
			PTL[1] = -0.976105557412198543;
			PTL[2] = -0.920649185347533874;
			PTL[3] = -0.835593535218090214;
			PTL[4] = -0.723679329283242681;
			PTL[5] = -0.588504834318661761;
			PTL[6] = -0.434415036912123975;
			PTL[7] = -0.26636265287828098;
			PTL[8] = -0.089749093484652111;
			PTL[9] = 0.089749093484652111;
			PTL[10] = 0.266362652878280984;
			PTL[11] = 0.434415036912123975;
			PTL[12] = 0.588504834318661761;
			PTL[13] = 0.723679329283242681;
			PTL[14] = 0.835593535218090214;
			PTL[15] = 0.920649185347533874;
			PTL[16] = 0.976105557412198543;
			PTL[17] = 1.;

			PTW[0] = 0.00653594771241830065;
			PTW[1] = 0.039970628810914066;
			PTW[2] = 0.070637166885633665;
			PTW[3] = 0.0990162717175028;
			PTW[4] = 0.1242105331329671;
			PTW[5] = 0.145411961573802268;
			PTW[6] = 0.161939517237602489;
			PTW[7] = 0.173262109489456226;
			PTW[8] = 0.17901586343970308;
			PTW[9] = 0.17901586343970308;
			PTW[10] = 0.173262109489456226;
			PTW[11] = 0.16193951723760249;
			PTW[12] = 0.145411961573802268;
			PTW[13] = 0.1242105331329671;
			PTW[14] = 0.099016271717502802;
			PTW[15] = 0.070637166885633665;
			PTW[16] = 0.03997062881091407;
			PTW[17] = 0.00653594771241830065;

			break;
		}
		case 19: {
			PTL[0] = -1.;
			PTL[1] = -0.978611766222080095;
			PTL[2] = -0.928901528152586244;
			PTL[3] = -0.852460577796646093;
			PTL[4] = -0.751494202552613014;
			PTL[5] = -0.628908137265220498;
			PTL[6] = -0.488229285680713503;
			PTL[7] = -0.33350484782449861;
			PTL[8] = -0.169186023409281571;
			PTL[9] = 0.;
			PTL[10] = 0.169186023409281571;
			PTL[11] = 0.33350484782449861;
			PTL[12] = 0.488229285680713503;
			PTL[13] = 0.628908137265220498;
			PTL[14] = 0.751494202552613014;
			PTL[15] = 0.852460577796646093;
			PTL[16] = 0.928901528152586244;
			PTL[17] = 0.978611766222080095;
			PTL[18] = 1.;

			PTW[0] = 0.00584795321637426901;
			PTW[1] = 0.035793365186176477;
			PTW[2] = 0.063381891762629737;
			PTW[3] = 0.0891317570992070845;
			PTW[4] = 0.11231534147730504;
			PTW[5] = 0.132267280448750777;
			PTW[6] = 0.14841394259593889;
			PTW[7] = 0.160290924044061242;
			PTW[8] = 0.167556584527142867;
			PTW[9] = 0.170001919284827235;
			PTW[10] = 0.16755658452714287;
			PTW[11] = 0.160290924044061242;
			PTW[12] = 0.14841394259593889;
			PTW[13] = 0.132267280448750777;
			PTW[14] = 0.112315341477305044;
			PTW[15] = 0.089131757099207084;
			PTW[16] = 0.063381891762629737;
			PTW[17] = 0.035793365186176477;
			PTW[18] = 0.00584795321637426901;

			break;
		}
		case 20:
		default: {
			PTL[0] = -1.;
			PTL[1] = -0.980743704893914172;
			PTL[2] = -0.935934498812665436;
			PTL[3] = -0.866877978089950141;
			PTL[4] = -0.77536826095205587;
			PTL[5] = -0.66377640229031129;
			PTL[6] = -0.53499286403188626;
			PTL[7] = -0.392353183713909299;
			PTL[8] = -0.239551705922986495;
			PTL[9] = -0.080545937238821838;
			PTL[10] = 0.080545937238821838;
			PTL[11] = 0.239551705922986495;
			PTL[12] = 0.392353183713909299;
			PTL[13] = 0.534992864031886262;
			PTL[14] = 0.66377640229031129;
			PTL[15] = 0.77536826095205587;
			PTL[16] = 0.866877978089950141;
			PTL[17] = 0.935934498812665436;
			PTL[18] = 0.980743704893914172;
			PTL[19] = 1.;

			PTW[0] = 0.00526315789473684211;
			PTW[1] = 0.032237123188488941;
			PTW[2] = 0.057181802127566826;
			PTW[3] = 0.080631763996119603;
			PTW[4] = 0.10199149969945082;
			PTW[5] = 0.12070922762867473;
			PTW[6] = 0.136300482358724185;
			PTW[7] = 0.148361554070916826;
			PTW[8] = 0.156580102647475487;
			PTW[9] = 0.16074328638784575;
			PTW[10] = 0.160743286387845749;
			PTW[11] = 0.156580102647475487;
			PTW[12] = 0.14836155407091683;
			PTW[13] = 0.13630048235872418;
			PTW[14] = 0.12070922762867473;
			PTW[15] = 0.101991499699450816;
			PTW[16] = 0.0806317639961196;
			PTW[17] = 0.057181802127566826;
			PTW[18] = 0.032237123188488941;
			PTW[19] = 0.00526315789473684211;

			break;
		}
		}
	} else if(intType == IntegrationType::RADAU) {
		switch(intOrder) {
		case 2: {
			PTL[0] = -1.;
			PTL[1] = 1. / 3.;

			PTW[0] = .5;
			PTW[1] = 1.5;
			break;
		}
		case 3: {
			PTL[0] = -1.;
			PTL[1] = -.2898979485;
			PTL[2] = .6898979485;

			PTW[0] = 2. / 9.;
			PTW[1] = 1.0249716523;
			PTW[2] = .7528061254;
			break;
		}
		case 4: {
			PTL[0] = -1.;
			PTL[1] = -.5753189235;
			PTL[2] = .1810662711;
			PTL[3] = .8228240809;

			PTW[0] = .125;
			PTW[1] = .6576886399;
			PTW[2] = .7763869376;
			PTW[3] = .4409244223;
			break;
		}
		case 5: {
			PTL[0] = -1.;
			PTL[1] = -.7204802713;
			PTL[2] = -.1671808647;
			PTL[3] = -.44631398727;
			PTL[4] = .8857916077;

			PTW[0] = .08;
			PTW[1] = .4462078021;
			PTW[2] = .6236530459;
			PTW[3] = .5627120302;
			PTW[4] = .2874271215;
			break;
		}
		case 6: {
			PTL[0] = -1.;
			PTL[1] = -.8029298284;
			PTL[2] = .3909285467;
			PTL[3] = .1240503795;
			PTL[4] = .6039731642;
			PTL[5] = .9203802858;

			PTW[0] = 1. / 18.;
			PTW[1] = .3196407532;
			PTW[2] = .4853871884;
			PTW[3] = .5209267831;
			PTW[4] = .4169013343;
			PTW[5] = .2015883852;
			break;
		}
		case 7: {
			PTL[0] = -1.;
			PTL[1] = -.8538913426;
			PTL[2] = -.5384677240;
			PTL[3] = -.1173430375;
			PTL[4] = .3260306194;
			PTL[5] = .7038428006;
			PTL[6] = .9413671456;

			PTW[0] = 2. / 49.;
			PTW[1] = .2392274892;
			PTW[2] = .3809498736;
			PTW[3] = .4471098290;
			PTW[4] = .4247037790;
			PTW[5] = .3182042314;
			PTW[6] = .1489884711;
			break;
		}
		case 8: {
			PTL[0] = -1.;
			PTL[1] = -.8874748789;
			PTL[2] = -.6395186165;
			PTL[3] = -.2947505657;
			PTL[4] = .0943072526;
			PTL[5] = .4684203544;
			PTL[6] = .7706418936;
			PTL[7] = .9550412271;

			PTW[0] = 1. / 32.;
			PTW[1] = .1853581548;
			PTW[2] = .3041306206;
			PTW[3] = .3765175453;
			PTW[4] = .3915721674;
			PTW[5] = .3470147956;
			PTW[6] = .2496479013;
			PTW[7] = .1145088147;
			break;
		}
		case 9: {
			PTL[0] = -1.;
			PTL[1] = -.9107320894;
			PTL[2] = -.7112674859;
			PTL[3] = -.4263504857;
			PTL[4] = -.0903733696;
			PTL[5] = .2561356708;
			PTL[6] = .5713830412;
			PTL[7] = .8173527842;
			PTL[8] = .9644401697;

			PTW[0] = 2. / 81.;
			PTW[1] = .1476540190;
			PTW[2] = .2471893782;
			PTW[3] = .3168437756;
			PTW[4] = .3482730027;
			PTW[5] = .3376939669;
			PTW[6] = .2863866963;
			PTW[7] = .2005532980;
			PTW[8] = .0907145049;
			break;
		}
		case 10:
		default: {
			PTL[0] = -1.;
			PTL[1] = -.9274843742;
			PTL[2] = -.7638420424;
			PTL[3] = -.5256460303;
			PTL[4] = -.2362344693;
			PTL[5] = .0760591978;
			PTL[6] = .3806648401;
			PTL[7] = .6477666876;
			PTL[8] = .8512252205;
			PTL[9] = .9711751807;

			PTW[0] = .02;
			PTW[1] = .1202966705;
			PTW[2] = .2042701318;
			PTW[3] = .2681948378;
			PTW[4] = .3058592877;
			PTW[5] = .3135824572;
			PTW[6] = .2906101648;
			PTW[7] = .2391934317;
			PTW[8] = .1643760127;
			PTW[9] = .0736170054;
			break;
		}
		}
	}
	auto IDX = 0;
	int_pts = new double*[n_rows];
	if(intDimension == 1)
		for(unsigned i = 0; i < intOrder; ++i) {
			int_pts[IDX] = new double[n_cols];
			int_pts[IDX][0] = PTL[i];
			int_pts[IDX++][1] = PTW[i];
		}
	else if(intDimension == 2)
		for(unsigned i = 0; i < intOrder; ++i)
			for(unsigned j = 0; j < intOrder; ++j) {
				int_pts[IDX] = new double[n_cols];
				int_pts[IDX][0] = PTL[i];
				int_pts[IDX][1] = PTL[j];
				int_pts[IDX++][2] = PTW[i] * PTW[j];
			}
	else if(intDimension == 3)
		for(unsigned i = 0; i < intOrder; ++i)
			for(unsigned j = 0; j < intOrder; ++j)
				for(unsigned k = 0; k < intOrder; ++k) {
					int_pts[IDX] = new double[n_cols];
					int_pts[IDX][0] = PTL[i];
					int_pts[IDX][1] = PTL[j];
					int_pts[IDX][2] = PTL[k];
					int_pts[IDX++][3] = PTW[i] * PTW[j] * PTW[k];
				}
	else throw std::invalid_argument("not supported");

	delete[] PTL;
	delete[] PTW;
}

IntegrationPlan::~IntegrationPlan() {
	for(unsigned i = 0; i < n_rows; ++i) if(int_pts[i] != nullptr) delete[] int_pts[i];
	delete[] int_pts;
}

double** IntegrationPlan::get_integration_scheme() const { return int_pts; }

double IntegrationPlan::operator()(const unsigned i, const unsigned j) const {
	if(i < n_rows && j < n_cols) return int_pts[i][j];
	throw std::invalid_argument("not supported");
}

/*
void IntegrationPlan::print() const {
	for(unsigned i = 0; i < n_rows; ++i) {
		printf("Node %u\t", i + 1);
		for(unsigned j = 0; j < n_cols - 1; ++j) printf("%+.6E\t", int_pts[i][j]);
		printf("Weight\t%+.6E\n", int_pts[i][n_cols - 1]);
	}
	printf("\n");
}
*/
