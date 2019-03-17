/*******************************************************************************
 * Copyright (C) 2017-2019 Theodore Chang
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
 * @class IntegrationPlan
 * @brief An IntegrationPlan class.
 *
 * This is a simple class provides users the integration points coordinates and
 * corresponding weights.
 *
 * Invoke the object with following definition:
 *
 *     IntegrationPlan TEST(intDimension, intOrder, intType);
 *
 * `intDimension` is the dimension of the problem.
 * `intOrder` is the integration order.
 * `intType` is the integration type. 1 is Gauss and 2 is Lobatto.
 *
 * The data is stored in `double** int_pts`.
 * The data can be acquired by using method `get_integration_scheme()`.
 * Single element could be accessed by `(int, int)` operator.
 * The data is arranged as follows.
 *
 *     |-------------------------------------|
 *     | POINT |   COORDINATES   |   WEIGHT  |
 *     |-------|-----|-----|-----|-----|-----|  <---
 *     |   1   | X_1 | Y_1 | ... |  WEIGHT_1 |    |
 *     |   2   | X_2 | Y_2 | ... |  WEIGHT_2 |    |
 *     |   3   | X_3 | Y_3 | ... |  WEIGHT_3 |    |
 *     |   .   |  .  |  .  |  .  |     .     |    |---ROWS=intOrder^intDimension
 *     |   .   |  .  |  .  |  .  |     .     |    |
 *     |   .   |  .  |  .  |  .  |     .     |    |
 *     |   N   | X_N | Y_N | ... |  WEIGHT_N |    |
 *     |-------------------------------------|  <---
 *             ^                             ^
 *             |-----------------------------|
 *                            |
 *                   COLS=intDimension+1
 *
 * @author tlc
 * @date 05/09/2017
 * @version 0.1.0
 * @file IntegrationPlan.h
 * @addtogroup Utility
 * @{
 */

#ifndef INTEGRATIONPLAN_H
#define INTEGRATIONPLAN_H

enum class IntegrationType { GAUSS, LOBATTO, RADAU, IRONS };

class IntegrationPlan final {
	double** int_pts;
public:
	const unsigned n_rows;
	const unsigned n_cols;
	const unsigned n_elem;

	explicit IntegrationPlan(unsigned = 2, unsigned = 2, IntegrationType = IntegrationType::GAUSS);
	IntegrationPlan(const IntegrationPlan&) = delete;
	IntegrationPlan(IntegrationPlan&&) = delete;
	IntegrationPlan& operator=(const IntegrationPlan&) = delete;
	IntegrationPlan& operator=(IntegrationPlan&&) = delete;
	~IntegrationPlan();

	double** get_integration_scheme() const;

	double operator()(unsigned, unsigned) const;

	// void print() const;
};

#endif

//! @}
