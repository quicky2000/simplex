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
    my_square_matrix(int dimension);

    my_square_matrix * get_transposed();

    double getDeterm();

    my_square_matrix *extract_square_matrix(unsigned int p_excluded_row_index,
                                        unsigned int p_excluded_column_index
                                       );


};

//-------------------------------------------------------------------------
my_square_matrix::my_square_matrix(int dimension)
:my_matrix(dimension,dimension)
{
}

//-------------------------------------------------------------------------
my_square_matrix * my_square_matrix::get_transposed()
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
double my_square_matrix::getDeterm()
{
    unsigned int l_height = this->get_height();
    if(2 == l_height)
    {
        return(get_data(0,0) * get_data(1,1) - get_data(0,1) * get_data(1,0));
    }

    double determ = 0;
    bool column = false;
    unsigned int max = 0;
    // Search line and column having maximum number of zero to speed up computation
    int lineMax = 0;
    int columnMax = 0;
    unsigned int * numberLine = new unsigned int[l_height];
    unsigned int * numberColumn = new unsigned int[l_height];

    // Array's initialization
    for(unsigned int l_index = 0; l_index < l_height; ++l_index)
    {
        numberLine[l_index] = 0;
        numberColumn[l_index] = 0;
    }

    // Number of zero computation
    for(unsigned int l_row_index = 0; l_row_index < l_height ; l_row_index++)
    {
        for (unsigned int l_column_index = 0; l_column_index < l_height ; ++l_column_index)
        {
            if(0 == get_data(l_row_index,l_column_index))
            {
                numberLine[l_row_index]++;
                numberColumn[l_column_index]++;
            }
        }
    }

    // Maximum's search
    for(unsigned int l_row_index = 0; l_row_index < l_height; l_row_index++)
    {
        if(numberLine[l_row_index] > max)
        {
            max = numberLine[l_row_index];
            lineMax = l_row_index;
        }
    }

    for(unsigned int l_column_index = 0; l_column_index < l_height; ++l_column_index)
    {
        if(numberColumn[l_column_index] > max)
        {
            max=numberColumn[l_column_index];
            columnMax=l_column_index;
            column=true;
        }
    }

    if(max == l_height)
    {
        return 0;
    }

    // Recursive computation
    if(column)
    {
        for(unsigned int l_row_index = 0; l_row_index < l_height; ++l_row_index)
        {
            determ += get_data(l_row_index, columnMax) * pow(-1, l_row_index) * pow(-1, columnMax) * (extract_square_matrix(l_row_index,
                                                                                                                            columnMax
                                                                                                                           ))->getDeterm();
        }
    }
    else
    {
        for(unsigned int l_column_index = 0; l_column_index < l_height; ++l_column_index)
        {
            determ += get_data(lineMax, l_column_index) * pow(-1, l_column_index) * pow(-1, lineMax) * (extract_square_matrix(lineMax,
                                                                                                                              l_column_index
                                                                                                                             ))->getDeterm();
        }
    }
    return(determ);

}

//-----------------------------------------------------------------------------
my_square_matrix * my_square_matrix::extract_square_matrix(unsigned int p_excluded_row_index,
                                                   unsigned int p_excluded_column_index
                                                  )
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
void test_square_matrix()
{
    my_square_matrix matrice(3);

    matrice.set_data(0,0,100);
    matrice.set_data(0,1,0);
    matrice.set_data(0,2,0);
    matrice.set_data(1,0,0);
    matrice.set_data(1,1,100);
    matrice.set_data(1,2,0);
    matrice.set_data(2,0,0);
    matrice.set_data(2,1,0);
    matrice.set_data(2,2,100);

    std::cout << "Matrix :" << std::endl;
    std::cout << matrice.to_string();
    std::cout << std::endl;

    std::cout << "Determ: " << matrice.getDeterm() << std::endl;

}
#endif // SIMPLEX_SELF_TEST


#endif // _MY_SQUARE_MATRIX_H_
// EOF
