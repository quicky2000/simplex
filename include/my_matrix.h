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

#ifndef _MY_MATRIX_H_
#define _MY_MATRIX_H_

#include "quicky_exception.h"
#include <string>
#include <cassert>
#include <cmath>
#include <iostream>

class my_matrix
{
  public:
	my_matrix(unsigned int p_height
	         ,unsigned int p_width
             );
	~my_matrix();

    void set_data(unsigned int p_row_index
                 ,unsigned int p_column_index
                 ,double p_value
                 );
	double get_data(unsigned int p_row_index
	               ,unsigned int p_column_index
                   ) const;
	unsigned int get_width() const;
	unsigned int get_height() const;
    void initialize(double a);

    /**
     * Create extracted matrix by excluding 1 line and 1 column
     * @param p_line_index index of line to exclude
     * @param p_column_index index of column to exclude
     * @return extracted matrix
     */
    my_matrix *
    extract_matrix(unsigned int p_line_index
                  ,unsigned int p_column_index
                  ) const;

    double* max() const;
    double* maxAbs() const;
    double* max_sub_matrix(unsigned int p_min_height,
                           unsigned int p_min_width
                          ) const;
    double* max_abs_sub_matrix(unsigned int p_min_height,
                               unsigned int p_min_width
                              )const ;
    double* max_column(unsigned int p_column) const;
    double* max_sub_column(unsigned int p_row_index
                          ,unsigned int p_column_index
                          ) const;
    double* max_abs_column(unsigned int p_column_index) const;
    double* max_abs_sub_column(unsigned int p_row_index
                              ,unsigned int p_column_index
                              ) const;
    void swap_line(unsigned int p_row_index_1,
                   unsigned int p_row_index_2
                  );
    void swap_column(unsigned int p_column_index_1,
                     unsigned int p_column_index_2
                    );
    my_matrix *
    mult(my_matrix & p_matrix);
    std::string to_string() const;

  private:
  protected:
    unsigned int m_width;
    unsigned int m_height;
    double * m_data;
};

//-------------------------------------------------------------------------
my_matrix::my_matrix(unsigned int p_height
                    ,unsigned int p_width
                    )
        :m_width(p_width)
        ,m_height(p_height)
        ,m_data(new double[p_height * p_width])
{
}

//-------------------------------------------------------------------------
void my_matrix::set_data(unsigned int p_row_index
                        ,unsigned int p_column_index
                        ,double p_value
                        )
{
    assert(p_row_index < m_height);
    assert(p_column_index < m_width);
    m_data[p_row_index * m_width + p_column_index] = p_value;
}

//-------------------------------------------------------------------------
double my_matrix::get_data(unsigned int p_row_index
                          ,unsigned int p_column_index
                          ) const
{
    assert(p_row_index < m_height);
    assert(p_column_index < m_width);
    return m_data[p_row_index * m_width + p_column_index];
}

//-------------------------------------------------------------------------
unsigned int my_matrix::get_width() const
{
    return m_width;
}

//-------------------------------------------------------------------------
unsigned int my_matrix::get_height() const
{
    return m_height;
}

//-------------------------------------------------------------------------
void my_matrix::initialize(double a)
{
    for(unsigned int i = 0; i < m_height; i++)
    {
        for(unsigned int j = 0; j < m_width; j++)
        {
            m_data[i * m_width + j] = a;
        }
    }
}

//-------------------------------------------------------------------------
my_matrix *
my_matrix::extract_matrix(unsigned int p_line_index
                         ,unsigned int p_column_index
                         ) const
{
    my_matrix * l_extracted_matrix = new my_matrix(m_height - 1, m_width - 1);
    unsigned int p_extracted_row_index = 0;
    unsigned int p_extracted_column_index = 0;

    for(unsigned int p_current_row_index = 0; p_current_row_index < m_height; p_current_row_index++)
    {
        for(unsigned int p_current_column_index = 0; p_current_column_index < m_width; p_current_column_index++)
        {
            if(p_current_row_index != p_line_index && p_current_column_index != p_column_index)
            {
                l_extracted_matrix->set_data(p_extracted_row_index, p_extracted_column_index, get_data(p_current_row_index,p_current_column_index));
            }
            if(p_current_column_index !=  p_column_index)
            {
                p_extracted_column_index++;
            }
        }
        p_extracted_column_index = 0;
        if(p_current_row_index != p_line_index)
        {
            p_extracted_row_index++;
        }
    }
    return(l_extracted_matrix);
}

//-------------------------------------------------------------------------
double* my_matrix::max() const
{
    return(max_sub_matrix(0, 0));
}

//-------------------------------------------------------------------------
double* my_matrix::maxAbs() const
{
    return(max_abs_sub_matrix(0, 0));
}

//-------------------------------------------------------------------------
double* my_matrix::max_sub_matrix(unsigned int p_min_height
                                 ,unsigned int p_min_width
                                 ) const
{
    assert(p_min_height < m_height);
    assert(p_min_width < m_width);

    double* l_result = new double[3];

    l_result[0] = get_data(p_min_height, p_min_width);
    l_result[1] = p_min_height;
    l_result[2] = p_min_width;

    for(unsigned int l_row_index = p_min_height; l_row_index < m_height; ++l_row_index)
    {
        for(unsigned l_column_index = p_min_width; l_column_index < m_width; ++l_column_index)
        {
            if(l_result[0] < get_data(l_row_index , l_column_index))
            {
                l_result[0] = get_data(l_row_index,l_column_index);
                l_result[1] = l_row_index;
                l_result[2] = l_column_index;
            }
        }
    }
    return(l_result);
}

//-------------------------------------------------------------------------
double* my_matrix::max_abs_sub_matrix(unsigned int p_min_height
                                     ,unsigned int p_min_width
                                     ) const
{
    assert(p_min_height < m_height);
    assert(p_min_width < m_width);

    double* l_result = new double[3];

    l_result[0] = std::abs(get_data(p_min_height, p_min_width));
    l_result[1] = p_min_height;
    l_result[2] = p_min_width;

    for(unsigned int l_row_index = p_min_height; l_row_index < m_height; l_row_index++)
    {
        for(unsigned int l_column_index = p_min_width; l_column_index < m_width; l_column_index++)
        {
            if(l_result[0] < std::abs(get_data(l_row_index, l_column_index)))
            {
                l_result[0] = std::abs(get_data(l_row_index, l_column_index));
                l_result[1] = l_row_index;
                l_result[2] = l_column_index;
            }
        }
    }
    return(l_result);
}

//-------------------------------------------------------------------------
double* my_matrix::max_column(unsigned int p_column) const
{
    return(max_sub_column(0, p_column));
}

//-------------------------------------------------------------------------
double* my_matrix::max_sub_column(unsigned int p_row_index
                                 ,unsigned int p_column_index
                                 ) const
{
    assert(p_row_index < m_height);
    assert(p_column_index < m_width);

    double * l_result = new double[2];
    l_result[0] = get_data(p_row_index, p_column_index);
    l_result[1] = p_row_index;

    for(unsigned int l_row_index = p_row_index; l_row_index < m_height; l_row_index++)
    {
        if(l_result[0] < get_data(l_row_index, p_column_index))
        {
            l_result[0] = get_data(l_row_index, p_column_index);
            l_result[1] = l_row_index;
        }
    }
    return(l_result);
}

//-------------------------------------------------------------------------
double* my_matrix::max_abs_column(unsigned int p_column_index) const
{
    return(max_abs_sub_column(0, p_column_index));
}

//-------------------------------------------------------------------------
double* my_matrix::max_abs_sub_column(unsigned int p_row_index
                                     ,unsigned int p_column_index
                                     ) const
{
    assert(p_row_index < m_height);
    assert(p_column_index < m_width);

    double* l_result = new double[2];
    l_result [0] = get_data(p_row_index, p_column_index);
    l_result [1] = p_row_index;

    for(unsigned int l_row_index = p_row_index ; l_row_index < m_height; l_row_index++)
    {
        if(l_result [0] < std::abs(get_data(l_row_index, p_column_index)))
        {
            l_result [0] = std::abs(get_data(l_row_index, p_column_index));
            l_result [1] = l_row_index;
        }
    }
    return(l_result );
}

//-------------------------------------------------------------------------
void my_matrix::swap_line(unsigned int p_row_index_1
                         ,unsigned int p_row_index_2
                         )
{
    assert(p_row_index_1 < m_height);
    assert(p_row_index_2 < m_height);

    double l_tmp_data;

    for(unsigned int l_column_index = 0 ; l_column_index < m_width; l_column_index++)
    {
        l_tmp_data= get_data(p_row_index_1, l_column_index);
        set_data(p_row_index_1, l_column_index, get_data(p_row_index_2, l_column_index));
        set_data(p_row_index_2, l_column_index, l_tmp_data);
    }
}

//-------------------------------------------------------------------------
void my_matrix::swap_column(unsigned int p_column_index_1
                           ,unsigned int p_column_index_2
                           )
{
    assert(p_column_index_1 < m_width);
    assert(p_column_index_2 < m_width);

    double l_tmp_data;

    for(unsigned int l_row_index = 0; l_row_index < m_height; l_row_index++)
    {
        l_tmp_data= get_data(l_row_index, p_column_index_1);
        set_data(l_row_index, p_column_index_1, get_data(l_row_index, p_column_index_2));
        set_data(l_row_index, p_column_index_2, l_tmp_data);
    }
}

//-------------------------------------------------------------------------
my_matrix *
my_matrix::mult(my_matrix & p_matrix)
{
    if(m_width != p_matrix.get_height())
    {
        throw quicky_exception::quicky_logic_exception("my_matrix.class: fatal error ! you try to multiplicate two matrix with incompatible sizes", __LINE__, __FILE__);
    }
    my_matrix* l_result = new my_matrix(m_height, p_matrix.get_width());
    double l_total;

    for(unsigned int l_row_index = 0; l_row_index < m_height; l_row_index++)
    {
        for(unsigned int l_column_index = 0; l_column_index < p_matrix.get_width(); l_column_index++)
        {
            l_total = 0;
            for(unsigned int l_mixed_index = 0; l_mixed_index < m_width; l_mixed_index++)
            {
                l_total+= get_data(l_row_index, l_mixed_index)* p_matrix.get_data(l_mixed_index, l_column_index);
            }
            l_result->set_data(l_row_index, l_column_index, l_total);
        }
    }
    return l_result;
}

//-------------------------------------------------------------------------
std::string
my_matrix::to_string() const
{
    std::string l_string("Width=");
    l_string += std::to_string(m_width) + "\n";
    l_string += "Height=" + std::to_string(m_height) + "\n";
    for(unsigned int l_row_index = 0; l_row_index < m_height; l_row_index++)
    {
        for(unsigned int l_column_index = 0;l_column_index < m_width; l_column_index++)
        {
            l_string += std::to_string(get_data(l_row_index, l_column_index)) + "\t";
        }
        l_string += "\n";
    }

    return l_string;
}

//-------------------------------------------------------------------------
my_matrix::~my_matrix()
{
    delete[] m_data;
}

#ifdef SIMPLEX_SELF_TEST
//-------------------------------------------------------------------------
void test_my_matrix()
{
    my_matrix l_matrix(3,4);
    my_matrix l_matrix2(4,3);

    l_matrix.set_data(0, 0, 1);
    l_matrix.set_data(0, 1, 2);
    l_matrix.set_data(0, 2, 2);
    l_matrix.set_data(0, 3, 3);
    l_matrix.set_data(1, 0, 5);
    l_matrix.set_data(1, 1, 5);
    l_matrix.set_data(1, 2, 6);
    l_matrix.set_data(1, 3, 7);
    l_matrix.set_data(2, 0, 5);
    l_matrix.set_data(2, 1, 9);
    l_matrix.set_data(2, 2, 10);
    l_matrix.set_data(2, 3, 11);

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

    std::cout << "matrix * matrix2 :" << std::endl << l_matrix.mult(l_matrix2)->to_string() << std::endl;
    std::cout << "matrix2 * matrix :" << std::endl << l_matrix2.mult(l_matrix)->to_string() << std::endl;
}
#endif // SIMPLEX_SELF_TEST

#endif // _MY_MATRIX_H_
// EOF
