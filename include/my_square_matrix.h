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

#ifndef _MY_SQUARE_MATRIX_H_
#define _MY_SQUARE_MATRIX_H_

#include "my_matrix.h"

class my_square_matrix: public my_matrix
{
  public:
    my_square_matrix(unsigned int dimension);

    my_square_matrix * get_transposed() const;

    double get_determ() const;

    my_square_matrix *extract_square_matrix(unsigned int p_excluded_row_index
                                           ,unsigned int p_excluded_column_index
                                           ) const;


};

//-------------------------------------------------------------------------
my_square_matrix::my_square_matrix(unsigned int dimension)
:my_matrix(dimension,dimension)
{
}

//-------------------------------------------------------------------------
my_square_matrix * my_square_matrix::get_transposed() const
{
    unsigned int l_height = this->get_height();
    my_square_matrix * l_transposed = new my_square_matrix(l_height);

    for(unsigned int l_row_index = 0; l_row_index < l_height; l_row_index++)
    {
        for(unsigned int l_column_index = 0; l_column_index < l_height; l_column_index++)
        {
            l_transposed->set_data(l_row_index, l_column_index, get_data(l_column_index,l_row_index));
        }
    }
    return(l_transposed);

}

//-------------------------------------------------------------------------
double my_square_matrix::get_determ() const
{
    unsigned int l_height = this->get_height();
    if(2 == l_height)
    {
        return(get_data(0,0) * get_data(1,1) - get_data(0,1) * get_data(1,0));
    }

    double l_determ = 0;
    bool l_column = false;
    unsigned int l_max = 0;
    // Search line and column having maximum number of zero to speed up computation
    int l_line_max = 0;
    int l_column_max = 0;
    unsigned int * l_number_line = new unsigned int[l_height];
    unsigned int * l_number_column = new unsigned int[l_height];

    // Array's initialization
    for(unsigned int l_index = 0; l_index < l_height; ++l_index)
    {
        l_number_line[l_index] = 0;
        l_number_column[l_index] = 0;
    }

    // Number of zero computation
    for(unsigned int l_row_index = 0; l_row_index < l_height ; l_row_index++)
    {
        for (unsigned int l_column_index = 0; l_column_index < l_height ; ++l_column_index)
        {
            if(0 == get_data(l_row_index,l_column_index))
            {
                l_number_line[l_row_index]++;
                l_number_column[l_column_index]++;
            }
        }
    }

    // Maximum's search
    for(unsigned int l_row_index = 0; l_row_index < l_height; l_row_index++)
    {
        if(l_number_line[l_row_index] > l_max)
        {
            l_max = l_number_line[l_row_index];
            l_line_max = l_row_index;
        }
    }
    delete[] l_number_line;

    for(unsigned int l_column_index = 0; l_column_index < l_height; ++l_column_index)
    {
        if(l_number_column[l_column_index] > l_max)
        {
            l_max=l_number_column[l_column_index];
            l_column_max=l_column_index;
            l_column=true;
        }
    }
    delete[] l_number_column;

    if(l_max == l_height)
    {
        return 0;
    }

    // Recursive computation
    if(l_column)
    {
        for(unsigned int l_row_index = 0; l_row_index < l_height; ++l_row_index)
        {
            my_square_matrix * l_extracted_matrix = extract_square_matrix(l_row_index, l_column_max);
            l_determ += get_data(l_row_index, l_column_max) * pow(-1, l_row_index) * pow(-1, l_column_max) * l_extracted_matrix->get_determ();
            delete l_extracted_matrix;
        }
    }
    else
    {
        for(unsigned int l_column_index = 0; l_column_index < l_height; ++l_column_index)
        {
            my_square_matrix * l_extracted_matrix = extract_square_matrix(l_line_max, l_column_index);
            l_determ += get_data(l_line_max, l_column_index) * pow(-1, l_column_index) * pow(-1, l_line_max) * l_extracted_matrix->get_determ();
            delete l_extracted_matrix;
        }
    }
    return(l_determ);

}

//-----------------------------------------------------------------------------
my_square_matrix * my_square_matrix::extract_square_matrix(unsigned int p_excluded_row_index
                                                          ,unsigned int p_excluded_column_index
                                                          ) const
{
    unsigned int l_height = this->get_height();
    my_square_matrix * extractedMatrix = new my_square_matrix(l_height - 1);
    unsigned int l_extracted_row_index = 0;
    unsigned int l_extracted_column_index = 0;

    for(unsigned int l_row_index = 0; l_row_index < l_height; ++l_row_index)
    {
        for(unsigned int l_column_index = 0; l_column_index < l_height; ++l_column_index)
        {
            if(l_row_index != p_excluded_row_index && l_column_index != p_excluded_column_index)
            {
                extractedMatrix->set_data(l_extracted_row_index, l_extracted_column_index, get_data(l_row_index, l_column_index));
            }
            if(l_column_index != p_excluded_column_index)
            {
                ++l_extracted_column_index;
            }
        }
        l_extracted_column_index = 0;
        if(l_row_index != p_excluded_row_index)
        {
            ++l_extracted_row_index;
        }
    }
    return extractedMatrix;
}

#ifdef SIMPLEX_SELF_TEST
bool test_square_matrix()
{
    bool l_ok = true;
    my_square_matrix l_matrix(3);
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

    my_square_matrix l_extracted_ref(2);
    l_extracted_ref.set_data(0, 0, 0.0);
    l_extracted_ref.set_data(0, 1, 8.0);
    l_extracted_ref.set_data(1, 0, 0.0);
    l_extracted_ref.set_data(1, 1, 0.0);
    my_square_matrix * l_extracted = l_matrix.extract_square_matrix(0, 1);
    l_ok &= quicky_utils::quicky_test::check_expected(*l_extracted == l_extracted_ref, true, "my_square_matrix::extract_square_matrix()");
    delete l_extracted;

    l_matrix.set_data(0,0,100);
    l_matrix.set_data(0,1,0);
    l_matrix.set_data(0,2,0);
    l_matrix.set_data(1,0,0);
    l_matrix.set_data(1,1,100);
    l_matrix.set_data(1,2,0);
    l_matrix.set_data(2,0,0);
    l_matrix.set_data(2,1,0);
    l_matrix.set_data(2,2,100);
    l_ok &= quicky_utils::quicky_test::check_expected(l_matrix.get_determ(), 1000000.0, "my_square_matrix::get_determ()");

    std::cout << "Matrix :" << std::endl;
    std::cout << l_matrix.to_string();
    std::cout << std::endl;

    std::cout << "Determ: " << l_matrix.get_determ() << std::endl;
    return l_ok;
}
#endif // SIMPLEX_SELF_TEST


#endif // _MY_SQUARE_MATRIX_H_
// EOF
