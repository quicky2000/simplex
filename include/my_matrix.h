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
#include <cstring>
#include <tuple>

class my_matrix
{
  public:
    my_matrix();
	my_matrix(unsigned int p_height
	         ,unsigned int p_width
             );
	~my_matrix();

	my_matrix(const my_matrix & p_matrix);
	my_matrix(my_matrix && p_matrix);

    my_matrix & operator=(const my_matrix & p_matrix);
    my_matrix & operator=(my_matrix && p_matrix);

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
    my_matrix
    extract_matrix(unsigned int p_line_index
                  ,unsigned int p_column_index
                  ) const;

    /**
     * Return max of matrix values and its coordinates
     * First tuple member is max value
     * Seconde tuple member is row index
     * Third tuple member is column index
     * @return max of matrix values
     */
    std::tuple<double, unsigned int, unsigned int> max() const;

    /**
     * Return max of matrix absolute values and its coordinates
     * First tuple member is max value
     * Seconde tuple member is row index
     * Third tuple member is column index
     * @return max of matrix values
     */
    std::tuple<double, unsigned int, unsigned int>
    max_abs() const;

    /**
     * Return max of sub matrix values and its coordinates
     * First tuple member is max value
     * Seconde tuple member is row index
     * Third tuple member is column index
     * @param p_min_height start row index
     * @param p_min_width start column index
     * @return max of matrix values
     */
    std::tuple<double, unsigned int, unsigned int>
    max_sub_matrix(unsigned int p_min_height
                  ,unsigned int p_min_width
                  ) const;

    /**
     * Return max of sub matrix absolute values and its coordinates
     * First tuple member is max value
     * Seconde tuple member is row index
     * Third tuple member is column index
     * @param p_min_height start row index
     * @param p_min_width start column index
     * @return max of matrix values
     */
    std::tuple<double, unsigned int, unsigned int>
    max_abs_sub_matrix(unsigned int p_min_height
                      ,unsigned int p_min_width
                      ) const ;

    /**
     * Return max of matrix column whose index is passed as parameter and it's
     * row index
     * First tuple member is max value
     * Seconde tuple member is row index
     * @param p_column_index column index
     * @return max value of column
     */
    std::tuple<double, unsigned int>
    max_column(unsigned int p_column_index) const;

    /**
     * Return max of matrix subcolumn whose index is passed as parameter and it's
     * row index
     * First tuple member is max value
     * Seconde tuple member is row index
     * @param p_row_index starting row index
     * @param p_column_index column index
     * @return max value of column
     */
    std::tuple<double, unsigned int>
    max_sub_column(unsigned int p_row_index
                  ,unsigned int p_column_index
                  ) const;

    /**
     * Return max of matrix column absolute values whose index is passed as
     * parameter and it's row index
     * First tuple member is max value
     * Seconde tuple member is row index
     * @param p_column_index column index
     * @return max value of column
     */
    std::tuple<double, unsigned int>
    max_abs_column(unsigned int p_column_index) const;

    /**
     * Return max of matrix subcolumn absolute values whose index is passed as
     * parameter and it's row index
     * First tuple member is max value
     * Seconde tuple member is row index
     * @param p_row_index starting row index
     * @param p_column_index column index
     * @return max value of column
     */
    std::tuple<double, unsigned int>
    max_abs_sub_column(unsigned int p_row_index
                      ,unsigned int p_column_index
                      ) const;

    void swap_line(unsigned int p_row_index_1
                  ,unsigned int p_row_index_2
                  );
    void swap_column(unsigned int p_column_index_1
                    ,unsigned int p_column_index_2
                    );
    my_matrix
    mult(const my_matrix & p_matrix);

    std::string to_string() const;

    bool operator==(const my_matrix & p_matrix) const;

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
my_matrix
my_matrix::extract_matrix(unsigned int p_line_index,
                          unsigned int p_column_index
                         ) const
{
    my_matrix l_extracted_matrix(m_height - 1, m_width - 1);
    unsigned int p_extracted_row_index = 0;
    unsigned int p_extracted_column_index = 0;

    for(unsigned int p_current_row_index = 0; p_current_row_index < m_height; p_current_row_index++)
    {
        for(unsigned int p_current_column_index = 0; p_current_column_index < m_width; p_current_column_index++)
        {
            if(p_current_row_index != p_line_index && p_current_column_index != p_column_index)
            {
                l_extracted_matrix.set_data(p_extracted_row_index, p_extracted_column_index, get_data(p_current_row_index,p_current_column_index));
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
std::tuple<double, unsigned int, unsigned int> my_matrix::max() const
{
    return(max_sub_matrix(0, 0));
}

//-------------------------------------------------------------------------
std::tuple<double, unsigned int, unsigned int>
my_matrix::max_abs() const
{
    return(max_abs_sub_matrix(0, 0));
}

//-------------------------------------------------------------------------
std::tuple<double, unsigned int, unsigned int>
my_matrix::max_sub_matrix(unsigned int p_min_height
                         ,unsigned int p_min_width
                         ) const
{
    assert(p_min_height < m_height);
    assert(p_min_width < m_width);

    std::tuple<double, unsigned int, unsigned int> l_result(get_data(p_min_height, p_min_width)
                                                           ,p_min_height
                                                           ,p_min_width
                                                           );

    for(unsigned int l_row_index = p_min_height; l_row_index < m_height; ++l_row_index)
    {
        for(unsigned l_column_index = p_min_width; l_column_index < m_width; ++l_column_index)
        {
            if(std::get<0>(l_result) < get_data(l_row_index , l_column_index))
            {
                l_result = std::make_tuple(get_data(l_row_index,l_column_index)
                                          ,l_row_index
                                          ,l_column_index
                                          );
            }
        }
    }
    return(l_result);
}

//-------------------------------------------------------------------------
std::tuple<double, unsigned int, unsigned int>
my_matrix::max_abs_sub_matrix(unsigned int p_min_height
                             ,unsigned int p_min_width
                             ) const
{
    assert(p_min_height < m_height);
    assert(p_min_width < m_width);

    std::tuple<double, unsigned int, unsigned int> l_result(std::abs(get_data(p_min_height, p_min_width))
                                                           ,p_min_height
                                                           ,p_min_width
                                                           );

    for(unsigned int l_row_index = p_min_height; l_row_index < m_height; l_row_index++)
    {
        for(unsigned int l_column_index = p_min_width; l_column_index < m_width; l_column_index++)
        {
            if(std::get<0>(l_result) < std::abs(get_data(l_row_index, l_column_index)))
            {
                l_result = std::make_tuple(std::abs(get_data(l_row_index, l_column_index))
                                          ,l_row_index
                                          ,l_column_index
                                          );
            }
        }
    }
    return(l_result);
}

//-------------------------------------------------------------------------
std::tuple<double, unsigned int>
my_matrix::max_column(unsigned int p_column_index) const
{
    return(max_sub_column(0, p_column_index));
}

//-------------------------------------------------------------------------
std::tuple<double, unsigned int>
my_matrix::max_sub_column(unsigned int p_row_index
                         ,unsigned int p_column_index
                         ) const
{
    assert(p_row_index < m_height);
    assert(p_column_index < m_width);

    std::tuple<double, unsigned int> l_result(get_data(p_row_index, p_column_index), p_row_index);

    for(unsigned int l_row_index = p_row_index; l_row_index < m_height; l_row_index++)
    {
        if(std::get<0>(l_result) < get_data(l_row_index, p_column_index))
        {
            l_result = std::make_tuple(get_data(l_row_index, p_column_index)
                                      ,l_row_index
                                      );
        }
    }
    return(l_result);
}

//-------------------------------------------------------------------------
std::tuple<double, unsigned int>
my_matrix::max_abs_column(unsigned int p_column_index) const
{
    return(max_abs_sub_column(0, p_column_index));
}

//-------------------------------------------------------------------------
std::tuple<double, unsigned int>
my_matrix::max_abs_sub_column(unsigned int p_row_index
                             ,unsigned int p_column_index
                             ) const
{
    assert(p_row_index < m_height);
    assert(p_column_index < m_width);

    std::tuple<double, unsigned int> l_result(std::abs(get_data(p_row_index, p_column_index))
                                             ,p_row_index
                                             );

    for(unsigned int l_row_index = p_row_index ; l_row_index < m_height; l_row_index++)
    {
        if(std::get<0>(l_result) < std::abs(get_data(l_row_index, p_column_index)))
        {
            l_result = std::make_tuple(std::abs(get_data(l_row_index, p_column_index))
                                      ,l_row_index
                                      );
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
my_matrix
my_matrix::mult(const my_matrix & p_matrix)
{
    if(m_width != p_matrix.get_height())
    {
        throw quicky_exception::quicky_logic_exception("my_matrix.class: fatal error ! you try to multiplicate two matrix with incompatible sizes", __LINE__, __FILE__);
    }
    my_matrix l_result(m_height, p_matrix.get_width());
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
            l_result.set_data(l_row_index, l_column_index, l_total);
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

my_matrix::my_matrix(const my_matrix & p_matrix)
:m_width(p_matrix.m_width)
,m_height(p_matrix.m_height)
,m_data(new double[m_width * m_height])
{
    memcpy(m_data, p_matrix.m_data, sizeof(double) * m_width * m_height);
}

//-------------------------------------------------------------------------
bool
my_matrix::operator==(const my_matrix & p_matrix) const
{
    if(m_height != p_matrix.m_height || m_width != p_matrix.m_width)
    {
        return false;
    }
    for(unsigned int l_row_index = 0; l_row_index < m_height; ++l_row_index)
    {
        for(unsigned int l_column_index = 0; l_column_index < m_width; ++l_column_index)
        {
            if(get_data(l_row_index, l_column_index) != p_matrix.get_data(l_row_index, l_column_index))
            {
                return false;
            }
        }
    }
    return true;
}

//-----------------------------------------------------------------------------
my_matrix::my_matrix(my_matrix && p_matrix)
:m_width(p_matrix.m_width)
,m_height(p_matrix.m_height)
,m_data(p_matrix.m_data)
{
    p_matrix.m_width = 0;
    p_matrix.m_height = 0;
    p_matrix.m_data = NULL;
}

//-----------------------------------------------------------------------------
my_matrix &
my_matrix::operator=(const my_matrix & p_matrix)
{
    m_width = p_matrix.m_width;
    m_height = p_matrix.m_height;
    delete[] m_data;
    m_data = new double[m_width * m_height];
    memcpy(m_data, p_matrix.m_data, m_width * m_height * sizeof(double));
    return *this;
}

//-----------------------------------------------------------------------------
my_matrix &
my_matrix::operator=(my_matrix && p_matrix)
{
    m_width = p_matrix.m_width;
    m_height = p_matrix.m_height;
    delete[] m_data;
    m_data = p_matrix.m_data;
    p_matrix.m_width = 0;
    p_matrix.m_height = 0;
    p_matrix.m_data = NULL;
    return *this;
}

//-----------------------------------------------------------------------------
my_matrix::my_matrix()
:m_width(0)
,m_height(0)
,m_data(NULL)
{

}

#ifdef SIMPLEX_SELF_TEST
//-----------------------------------------------------------------------------
bool check_max(const std::tuple<double, unsigned int, unsigned int> & p_max
              ,const std::tuple<double, unsigned int, unsigned int> & p_ref
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
bool check_max(const std::tuple<double, unsigned int> & p_max
              ,const std::tuple<double, unsigned int> & p_ref
              ,const std::string & p_method_name
              )
{
    bool l_ok = true;
    l_ok &= quicky_utils::quicky_test::check_expected(std::get<0>(p_max), std::get<0>(p_ref), p_method_name + " value");
    l_ok &= quicky_utils::quicky_test::check_expected(std::get<1>(p_max), std::get<1>(p_ref), p_method_name + " row/column index");
    return l_ok;
}

//-------------------------------------------------------------------------
bool test_my_matrix()
{
    bool l_ok = true;
    my_matrix l_matrix(3,4);
    my_matrix l_matrix2(4,3);

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

    l_ok &= quicky_utils::quicky_test::check_expected(l_matrix.get_data(0, 0), -17.0 ,"my_matrix::get_data()");
    l_ok &= quicky_utils::quicky_test::check_expected(l_matrix.get_data(2, 0), 5.0 ,"my_matrix::get_data()");
    l_ok &= quicky_utils::quicky_test::check_expected(l_matrix.get_data(2, 3), -11.0 ,"my_matrix::get_data()");
    l_ok &= quicky_utils::quicky_test::check_expected(l_matrix.get_data(0, 3), 3.0 ,"my_matrix::get_data()");
    l_ok &= quicky_utils::quicky_test::check_expected(l_matrix.get_width(), 4u ,"my_matrix::get_width()");
    l_ok &= quicky_utils::quicky_test::check_expected(l_matrix.get_height(), 3u, "my_matrix::get_height()");
    l_ok &= check_max(l_matrix.max(), std::make_tuple(15.2, 0u, 1u), "my_matrix::max");
    l_ok &= check_max(l_matrix.max_abs(), std::make_tuple(17.0, 0u, 0u), "my_matrix::max_abs");
    l_ok &= check_max(l_matrix.max_sub_matrix(0,1), std::make_tuple(15.2, 0u, 1u), "my_matrix::max_sub_matrix");
    l_ok &= check_max(l_matrix.max_sub_matrix(1,0), std::make_tuple(9.0, 2u, 1u), "my_matrix::max_sub_matrix");
    l_ok &= check_max(l_matrix.max_abs_sub_matrix(0,1), std::make_tuple(15.2, 0u, 1u), "my_matrix::max_abs_sub_matrix");
    l_ok &= check_max(l_matrix.max_abs_sub_matrix(1,0), std::make_tuple(11.0, 2u, 3u), "my_matrix::max_abs_sub_matrix");
    l_ok &= check_max(l_matrix.max_column(0), std::make_tuple(5.0, 1u), "my_matrix::max_column");
    l_ok &= check_max(l_matrix.max_column(3), std::make_tuple(7.0, 1u), "my_matrix::max_column");
    l_ok &= check_max(l_matrix.max_abs_column(0), std::make_tuple(17.0, 0u), "my_matrix::max_abs_column");
    l_ok &= check_max(l_matrix.max_abs_column(3), std::make_tuple(11.0, 2u), "my_matrix::max_abs_column");
    l_ok &= check_max(l_matrix.max_abs_sub_column(1, 0), std::make_tuple(5.0, 1u), "my_matrix::max_abs_sub_column");
    l_ok &= check_max(l_matrix.max_sub_column(1, 2), std::make_tuple(6.0, 1u), "my_matrix::max_sub_column");

    {
        my_matrix l_matrix3 = l_matrix.extract_matrix(1,2);
        my_matrix l_matrix_ref(2,3);
        l_matrix_ref.set_data(0, 0, -17.0);
        l_matrix_ref.set_data(0, 1, 15.2);
        l_matrix_ref.set_data(0, 2, 3.0);
        l_matrix_ref.set_data(1, 0, 5.0);
        l_matrix_ref.set_data(1, 1, 9.0);
        l_matrix_ref.set_data(1, 2, -11.0);
        l_ok &= quicky_utils::quicky_test::check_expected(l_matrix_ref == l_matrix3, true, "my_matrix::extract_matrix()");
    }
    {
        my_matrix l_matrix_ref(3,4);
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
        my_matrix l_matrix_copy(l_matrix);
        l_matrix_copy.swap_line(1, 2);
        l_ok &= quicky_utils::quicky_test::check_expected(l_matrix_ref == l_matrix_copy, true, "my_matrix::swap_line()");
    }
    {
        my_matrix l_matrix_ref(3,4);
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
        my_matrix l_matrix_copy(l_matrix);
        l_matrix_copy.swap_column(1, 2);
        l_ok &= quicky_utils::quicky_test::check_expected(l_matrix_ref == l_matrix_copy, true, "my_matrix::swap_column()");
    }
    {
        my_matrix l_op1(2,3);
        l_op1.set_data(0, 0, 1.0);
        l_op1.set_data(0, 1, 2.0);
        l_op1.set_data(0, 2, 3.0);
        l_op1.set_data(1, 0, 4.0);
        l_op1.set_data(1, 1, 5.0);
        l_op1.set_data(1, 2, 6.0);

        my_matrix l_op2(3,1);
        l_op2.set_data(0, 0, 0.5);
        l_op2.set_data(1, 0, 1.5);
        l_op2.set_data(2, 0, 2.5);

        my_matrix l_ref_result(2,1);
        l_ref_result.set_data(0, 0, 11.0);
        l_ref_result.set_data(1, 0, 24.5);
        my_matrix l_mult = l_op1.mult(l_op2);
        l_ok &= quicky_utils::quicky_test::check_expected(l_ref_result == l_mult, true, "my_matrix::mult()");
    }
    return l_ok;
}
#endif // SIMPLEX_SELF_TEST

#endif // _MY_MATRIX_H_
// EOF
