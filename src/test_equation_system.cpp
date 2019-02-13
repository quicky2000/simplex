/*    This file is part of simplex
      Copyright (C) 2019  Julien Thevenon ( julien_thevenon at yahoo.fr )

      This program is free software: you can redistribute it and/or modify
      it under the terms of the GNU General Public License as published by
      the Free Software Foundation, either version 3 of the License, or
      (at your option) any later version.

      This program is distributed in the hope that it will be useful,
      but WITHOUT ANY WARRANTY; without even the implied warranty of
      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
      GNU General Public License for more details.

      You should have received a copy of the GNU General Public License
      along with this program.  If not, see <http://www.gnu.org/licenses/>
*/

#include "equation_system.h"
#include "quicky_test.h"

#ifdef SIMPLEX_SELF_TEST
//-----------------------------------------------------------------------------
bool test_equation_system()
{
    bool l_ok =true;
    my_square_matrix<double> l_matrix(3);
    l_matrix.set_data(0,0,1);
    l_matrix.set_data(0,1,1);
    l_matrix.set_data(0,2,1);
    l_matrix.set_data(1,0,-1);
    l_matrix.set_data(1,1,1);
    l_matrix.set_data(1,2,1);
    l_matrix.set_data(2,0,-1);
    l_matrix.set_data(2,1,-1);
    l_matrix.set_data(2,2,1);

    my_matrix<double> l_coef(3,1);
    l_coef.set_data(0,0,6);
    l_coef.set_data(1,0,4);
    l_coef.set_data(2,0,0);

    SystemEquation<double> l_system(l_matrix,l_coef);
    my_matrix<double> l_result = l_system.solve();

    assert(l_result.get_height());

    l_ok &= quicky_utils::quicky_test::check_expected(l_result.get_data(0, 0), 1.0, "Variable[0]");
    l_ok &= quicky_utils::quicky_test::check_expected(l_result.get_data(1, 0), 2.0, "Variable[1]");
    l_ok &= quicky_utils::quicky_test::check_expected(l_result.get_data(2, 0), 3.0, "Variable[2]");

    std::cout << "Equation Matrix :" << std::endl;
    std::cout << l_matrix.to_string();

    std::cout << "Coef matrix:" << std::endl;
    std::cout << l_coef.to_string();

    std::cout << "Equation system:" << std::endl;
    std::cout << l_system.to_string();


    for(unsigned int i = 0; i < l_result.get_height();++i)
    {
        std::cout << "Result[" << std::to_string(i) << "] : " << l_result.get_data(i, 0) << std::endl;
    }

    l_matrix.set_data(2,0,0);
    l_matrix.set_data(2,1,0);
    l_matrix.set_data(2,2,0);
    SystemEquation<double> l_system2 = SystemEquation<double>(l_matrix,l_coef);
    l_result = l_system2.solve();
    l_ok &= quicky_utils::quicky_test::check_expected(l_result.get_height(), 0u, "equation_system::solve() No result");
    return l_ok;
}
#endif // SIMPLEX_SELF_TEST

//EOF
