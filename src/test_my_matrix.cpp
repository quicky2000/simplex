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

#include "my_matrix.h"
#include "quicky_test.h"
#include <iostream>

#ifdef SIMPLEX_SELF_TEST
namespace simplex
{
    //-----------------------------------------------------------------------------
    template <typename T>
    bool
    check_max(const std::tuple<T, unsigned int, unsigned int> & p_max
             ,const std::tuple<T, unsigned int, unsigned int> & p_ref
             ,const std::string & p_method_name
             )
    {
        bool l_ok = true;
        l_ok &= quicky_utils::quicky_test::check_expected(std::get<0>(p_max), std::get<0>(p_ref), p_method_name + " value");
        l_ok &= quicky_utils::quicky_test::check_expected(std::get<1>(p_max), std::get<1>(p_ref), p_method_name + " row index");
        l_ok &= quicky_utils::quicky_test::check_expected(std::get<2>(p_max), std::get<2>(p_ref), p_method_name + " column index");
        return l_ok;
    }

    //-----------------------------------------------------------------------------
    template <typename T>
    bool
    check_max(const std::tuple<T, unsigned int> & p_max
             ,const std::tuple<T, unsigned int> & p_ref
             ,const std::string & p_method_name
             )
    {
        bool l_ok = true;
        l_ok &= quicky_utils::quicky_test::check_expected(std::get<0>(p_max), std::get<0>(p_ref), p_method_name + " value");
        l_ok &= quicky_utils::quicky_test::check_expected(std::get<1>(p_max), std::get<1>(p_ref), p_method_name + " row/column index");
        return l_ok;
    }

    //-------------------------------------------------------------------------
    bool
    test_my_matrix()
    {
        bool l_ok = true;
        my_matrix<double> l_matrix(3, 4);
        my_matrix<double> l_matrix2(4, 3);

        l_matrix.set_data(0, 0, -17);
        l_matrix.set_data(0, 1, 15.2);
        l_matrix.set_data(0, 2, 10);
        l_matrix.set_data(0, 3, 3);
        l_matrix.set_data(1, 0, 5);
        l_matrix.set_data(1, 1, 5);
        l_matrix.set_data(1, 2, 6);
        l_matrix.set_data(1, 3, 7);
        l_matrix.set_data(2, 0, 5);
        l_matrix.set_data(2, 1, 9);
        l_matrix.set_data(2, 2, 2);
        l_matrix.set_data(2, 3, -11);

        l_matrix2.set_data(0, 0, 1);
        l_matrix2.set_data(0, 1, 2);
        l_matrix2.set_data(0, 2, 3);
        l_matrix2.set_data(1, 0, 5);
        l_matrix2.set_data(1, 1, 7);
        l_matrix2.set_data(1, 2, 11);
        l_matrix2.set_data(2, 0, 13);
        l_matrix2.set_data(2, 1, 17);
        l_matrix2.set_data(2, 2, 19);
        l_matrix2.set_data(3, 0, 23);
        l_matrix2.set_data(3, 1, 29);
        l_matrix2.set_data(3, 2, 31);

        std::cout << "matrix  :" << std::endl << l_matrix.to_string() << std::endl;
        std::cout << "matrix2 :" << std::endl << l_matrix2.to_string() << std::endl;

        l_ok &= quicky_utils::quicky_test::check_expected(l_matrix == l_matrix2, false, "my_matrix::operator==()");
        l_ok &= quicky_utils::quicky_test::check_expected(l_matrix == l_matrix, true, "my_matrix::operator==()");

        l_ok &= quicky_utils::quicky_test::check_expected(l_matrix.get_data(0, 0), -17.0, "my_matrix::get_data()");
        l_ok &= quicky_utils::quicky_test::check_expected(l_matrix.get_data(2, 0), 5.0, "my_matrix::get_data()");
        l_ok &= quicky_utils::quicky_test::check_expected(l_matrix.get_data(2, 3), -11.0, "my_matrix::get_data()");
        l_ok &= quicky_utils::quicky_test::check_expected(l_matrix.get_data(0, 3), 3.0, "my_matrix::get_data()");
        l_ok &= quicky_utils::quicky_test::check_expected(l_matrix.get_width(), 4u, "my_matrix::get_width()");
        l_ok &= quicky_utils::quicky_test::check_expected(l_matrix.get_height(), 3u, "my_matrix::get_height()");
        l_ok &= check_max(l_matrix.max(), std::make_tuple(15.2, 0u, 1u), "my_matrix::max");
        l_ok &= check_max(l_matrix.max_abs(), std::make_tuple(17.0, 0u, 0u), "my_matrix::max_abs");
        l_ok &= check_max(l_matrix.max_sub_matrix(0, 1), std::make_tuple(15.2, 0u, 1u), "my_matrix::max_sub_matrix");
        l_ok &= check_max(l_matrix.max_sub_matrix(1, 0), std::make_tuple(9.0, 2u, 1u), "my_matrix::max_sub_matrix");
        l_ok &= check_max(l_matrix.max_abs_sub_matrix(0, 1), std::make_tuple(15.2, 0u, 1u), "my_matrix::max_abs_sub_matrix");
        l_ok &= check_max(l_matrix.max_abs_sub_matrix(1, 0), std::make_tuple(11.0, 2u, 3u), "my_matrix::max_abs_sub_matrix");
        l_ok &= check_max(l_matrix.max_column(0), std::make_tuple(5.0, 1u), "my_matrix::max_column");
        l_ok &= check_max(l_matrix.max_column(3), std::make_tuple(7.0, 1u), "my_matrix::max_column");
        l_ok &= check_max(l_matrix.max_abs_column(0), std::make_tuple(17.0, 0u), "my_matrix::max_abs_column");
        l_ok &= check_max(l_matrix.max_abs_column(3), std::make_tuple(11.0, 2u), "my_matrix::max_abs_column");
        l_ok &= check_max(l_matrix.max_abs_sub_column(1, 0), std::make_tuple(5.0, 1u), "my_matrix::max_abs_sub_column");
        l_ok &= check_max(l_matrix.max_sub_column(1, 2), std::make_tuple(6.0, 1u), "my_matrix::max_sub_column");

        {
            my_matrix<double> l_matrix3 = l_matrix.extract_matrix(1, 2);
            my_matrix<double> l_matrix_ref(2, 3);
            l_matrix_ref.set_data(0, 0, -17.0);
            l_matrix_ref.set_data(0, 1, 15.2);
            l_matrix_ref.set_data(0, 2, 3.0);
            l_matrix_ref.set_data(1, 0, 5.0);
            l_matrix_ref.set_data(1, 1, 9.0);
            l_matrix_ref.set_data(1, 2, -11.0);
            l_ok &= quicky_utils::quicky_test::check_expected(l_matrix_ref == l_matrix3, true, "my_matrix::extract_matrix()");
        }
        {
            my_matrix<double> l_matrix_ref(3, 4);
            l_matrix_ref.set_data(0, 0, -17);
            l_matrix_ref.set_data(0, 1, 15.2);
            l_matrix_ref.set_data(0, 2, 10);
            l_matrix_ref.set_data(0, 3, 3);
            l_matrix_ref.set_data(1, 0, 5);
            l_matrix_ref.set_data(1, 1, 9);
            l_matrix_ref.set_data(1, 2, 2);
            l_matrix_ref.set_data(1, 3, -11);
            l_matrix_ref.set_data(2, 0, 5);
            l_matrix_ref.set_data(2, 1, 5);
            l_matrix_ref.set_data(2, 2, 6);
            l_matrix_ref.set_data(2, 3, 7);
            my_matrix<double> l_matrix_copy(l_matrix);
            l_matrix_copy.swap_line(1, 2);
            l_ok &= quicky_utils::quicky_test::check_expected(l_matrix_ref == l_matrix_copy, true, "my_matrix::swap_line()");
        }
        {
            my_matrix<double> l_matrix_ref(3, 4);
            l_matrix_ref.set_data(0, 0, -17);
            l_matrix_ref.set_data(0, 1, 10);
            l_matrix_ref.set_data(0, 2, 15.2);
            l_matrix_ref.set_data(0, 3, 3);
            l_matrix_ref.set_data(1, 0, 5);
            l_matrix_ref.set_data(1, 1, 6);
            l_matrix_ref.set_data(1, 2, 5);
            l_matrix_ref.set_data(1, 3, 7);
            l_matrix_ref.set_data(2, 0, 5);
            l_matrix_ref.set_data(2, 1, 2);
            l_matrix_ref.set_data(2, 2, 9);
            l_matrix_ref.set_data(2, 3, -11);
            my_matrix<double> l_matrix_copy(l_matrix);
            l_matrix_copy.swap_column(1, 2);
            l_ok &= quicky_utils::quicky_test::check_expected(l_matrix_ref == l_matrix_copy, true, "my_matrix::swap_column()");
        }
        {
            my_matrix<double> l_op1(2, 3);
            l_op1.set_data(0, 0, 1.0);
            l_op1.set_data(0, 1, 2.0);
            l_op1.set_data(0, 2, 3.0);
            l_op1.set_data(1, 0, 4.0);
            l_op1.set_data(1, 1, 5.0);
            l_op1.set_data(1, 2, 6.0);

            my_matrix<double> l_op2(3, 1);
            l_op2.set_data(0, 0, 0.5);
            l_op2.set_data(1, 0, 1.5);
            l_op2.set_data(2, 0, 2.5);

            my_matrix<double> l_ref_result(2, 1);
            l_ref_result.set_data(0, 0, 11.0);
            l_ref_result.set_data(1, 0, 24.5);
            my_matrix<double> l_mult = l_op1.mult(l_op2);
            l_ok &= quicky_utils::quicky_test::check_expected(l_ref_result == l_mult, true, "my_matrix::mult()");
        }
        return l_ok;
    }
}
#endif // SIMPLEX_SELF_TEST
//EOF
