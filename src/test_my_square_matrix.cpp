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

#include "my_square_matrix.h"
#include "quicky_test.h"

#ifdef SIMPLEX_SELF_TEST
namespace simplex
{
    //-----------------------------------------------------------------------------
    bool
    test_square_matrix()
    {
        bool l_ok = true;
        my_square_matrix<double> l_matrix(3);
        l_matrix.set_data(0, 0, -1.0);
        l_matrix.set_data(0, 1, 2.0);
        l_matrix.set_data(0, 2, 5.0);
        l_matrix.set_data(1, 0, 1.0);
        l_matrix.set_data(1, 1, 2.0);
        l_matrix.set_data(1, 2, 3.0);
        l_matrix.set_data(2, 0, -2.0);
        l_matrix.set_data(2, 1, 8.0);
        l_matrix.set_data(2, 2, 10.0);
        l_ok &= quicky_utils::quicky_test::check_expected(l_matrix.get_determ(), 32.0, "my_square_matrix::get_determ()");

        l_matrix.set_data(1, 0, 0.0);
        l_matrix.set_data(1, 1, 4.0);
        l_matrix.set_data(1, 2, 8.0);
        l_matrix.set_data(2, 0, 0.0);
        l_matrix.set_data(2, 1, 4.0);
        l_matrix.set_data(2, 2, 0.0);
        l_ok &= quicky_utils::quicky_test::check_expected(l_matrix.get_determ(), 32.0, "my_square_matrix::get_determ()");

        my_square_matrix<double> l_extracted_ref(2);
        l_extracted_ref.set_data(0, 0, 0.0);
        l_extracted_ref.set_data(0, 1, 8.0);
        l_extracted_ref.set_data(1, 0, 0.0);
        l_extracted_ref.set_data(1, 1, 0.0);
        my_square_matrix<double> l_extracted = l_matrix.extract_square_matrix(0, 1);
        l_ok &= quicky_utils::quicky_test::check_expected(l_extracted == l_extracted_ref, true, "my_square_matrix::extract_square_matrix()");

        l_matrix.set_data(0, 0, 100);
        l_matrix.set_data(0, 1, 0);
        l_matrix.set_data(0, 2, 0);
        l_matrix.set_data(1, 0, 0);
        l_matrix.set_data(1, 1, 100);
        l_matrix.set_data(1, 2, 0);
        l_matrix.set_data(2, 0, 0);
        l_matrix.set_data(2, 1, 0);
        l_matrix.set_data(2, 2, 100);
        l_ok &= quicky_utils::quicky_test::check_expected(l_matrix.get_determ(), 1000000.0, "my_square_matrix::get_determ()");

        {
            my_square_matrix<double> l_small_matrix(2);
            l_small_matrix.set_data(0, 0, 0);
            l_small_matrix.set_data(0, 1, 1);
            l_small_matrix.set_data(1, 0, 2);
            l_small_matrix.set_data(1, 1, 3);
            my_square_matrix<double> l_transposed = l_small_matrix.get_transposed();

            my_square_matrix<double> l_reference(2);
            l_reference.set_data(0, 0, 0);
            l_reference.set_data(1, 0, 1);
            l_reference.set_data(0, 1, 2);
            l_reference.set_data(1, 1, 3);

            l_ok &= quicky_utils::quicky_test::check_expected(l_transposed == l_reference, true, "my_square_matrix::get_transposed()");
        }
        return l_ok;
    }
}
#endif // SIMPLEX_SELF_TEST

//EOF
