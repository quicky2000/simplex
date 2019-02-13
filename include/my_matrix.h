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
#include <cstring>
#include <tuple>

template <typename T>
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
                 ,const T & p_value
                 );
	const T &
    get_data(unsigned int p_row_index,
             unsigned int p_column_index
            ) const;
	unsigned int get_width() const;
	unsigned int get_height() const;
    void initialize(const T & a);

    /**
     * Create extracted matrix by excluding 1 line and 1 column
     * @param p_line_index index of line to exclude
     * @param p_column_index index of column to exclude
     * @return extracted matrix
     */
    my_matrix<T>
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
    std::tuple<T, unsigned int, unsigned int> max() const;

    /**
     * Return max of matrix absolute values and its coordinates
     * First tuple member is max value
     * Seconde tuple member is row index
     * Third tuple member is column index
     * @return max of matrix values
     */
    std::tuple<T, unsigned int, unsigned int>
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
    std::tuple<T, unsigned int, unsigned int>
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
    std::tuple<T, unsigned int, unsigned int>
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
    std::tuple<T, unsigned int>
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
    std::tuple<T, unsigned int>
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
    std::tuple<T, unsigned int>
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
    std::tuple<T, unsigned int>
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

    typedef T coef_type_t;

  private:
    unsigned int m_width;
    unsigned int m_height;
    T * m_data;
};

//-------------------------------------------------------------------------
template <typename T>
my_matrix<T>::my_matrix(unsigned int p_height
                       ,unsigned int p_width
                       )
        :m_width(p_width)
        ,m_height(p_height)
        ,m_data(new T[p_height * p_width])
{
}

//-------------------------------------------------------------------------
template <typename T>
void my_matrix<T>::set_data(unsigned int p_row_index
                           ,unsigned int p_column_index
                           ,const T & p_value
                           )
{
    assert(p_row_index < m_height);
    assert(p_column_index < m_width);
    m_data[p_row_index * m_width + p_column_index] = p_value;
}

//-------------------------------------------------------------------------
template <typename T>
const T &
my_matrix<T>::get_data(unsigned int p_row_index
                      ,unsigned int p_column_index
                      ) const
{
    assert(p_row_index < m_height);
    assert(p_column_index < m_width);
    return m_data[p_row_index * m_width + p_column_index];
}

//-------------------------------------------------------------------------
template <typename T>
unsigned int my_matrix<T>::get_width() const
{
    return m_width;
}

//-------------------------------------------------------------------------
template <typename T>
unsigned int my_matrix<T>::get_height() const
{
    return m_height;
}

//-------------------------------------------------------------------------
template <typename T>
void my_matrix<T>::initialize(const T & a)
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
template <typename T>
my_matrix<T>
my_matrix<T>::extract_matrix(unsigned int p_line_index
                            ,unsigned int p_column_index
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
template <typename T>
std::tuple<T, unsigned int, unsigned int> my_matrix<T>::max() const
{
    return(max_sub_matrix(0, 0));
}

//-------------------------------------------------------------------------
template <typename T>
std::tuple<T, unsigned int, unsigned int>
my_matrix<T>::max_abs() const
{
    return(max_abs_sub_matrix(0, 0));
}

//-------------------------------------------------------------------------
template <typename T>
std::tuple<T, unsigned int, unsigned int>
my_matrix<T>::max_sub_matrix(unsigned int p_min_height
                            ,unsigned int p_min_width
                            ) const
{
    assert(p_min_height < m_height);
    assert(p_min_width < m_width);

    std::tuple<T, unsigned int, unsigned int> l_result(get_data(p_min_height, p_min_width)
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
template <typename T>
std::tuple<T, unsigned int, unsigned int>
my_matrix<T>::max_abs_sub_matrix(unsigned int p_min_height
                                ,unsigned int p_min_width
                                ) const
{
    assert(p_min_height < m_height);
    assert(p_min_width < m_width);

    std::tuple<T, unsigned int, unsigned int> l_result(std::abs(get_data(p_min_height, p_min_width))
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
template <typename T>
std::tuple<T, unsigned int>
my_matrix<T>::max_column(unsigned int p_column_index) const
{
    return(max_sub_column(0, p_column_index));
}

//-------------------------------------------------------------------------
template <typename T>
std::tuple<T, unsigned int>
my_matrix<T>::max_sub_column(unsigned int p_row_index
                            ,unsigned int p_column_index
                            ) const
{
    assert(p_row_index < m_height);
    assert(p_column_index < m_width);

    std::tuple<T, unsigned int> l_result(get_data(p_row_index, p_column_index), p_row_index);

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
template <typename T>
std::tuple<T, unsigned int>
my_matrix<T>::max_abs_column(unsigned int p_column_index) const
{
    return(max_abs_sub_column(0, p_column_index));
}

//-------------------------------------------------------------------------
template <typename T>
std::tuple<T, unsigned int>
my_matrix<T>::max_abs_sub_column(unsigned int p_row_index
                                ,unsigned int p_column_index
                                ) const
{
    assert(p_row_index < m_height);
    assert(p_column_index < m_width);

    std::tuple<T, unsigned int> l_result(std::abs(get_data(p_row_index, p_column_index))
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
template <typename T>
void my_matrix<T>::swap_line(unsigned int p_row_index_1
                            ,unsigned int p_row_index_2
                            )
{
    assert(p_row_index_1 < m_height);
    assert(p_row_index_2 < m_height);

    T l_tmp_data;

    for(unsigned int l_column_index = 0 ; l_column_index < m_width; l_column_index++)
    {
        l_tmp_data= get_data(p_row_index_1, l_column_index);
        set_data(p_row_index_1, l_column_index, get_data(p_row_index_2, l_column_index));
        set_data(p_row_index_2, l_column_index, l_tmp_data);
    }
}

//-------------------------------------------------------------------------
template <typename T>
void my_matrix<T>::swap_column(unsigned int p_column_index_1
                              ,unsigned int p_column_index_2
                              )
{
    assert(p_column_index_1 < m_width);
    assert(p_column_index_2 < m_width);

    T l_tmp_data;

    for(unsigned int l_row_index = 0; l_row_index < m_height; l_row_index++)
    {
        l_tmp_data= get_data(l_row_index, p_column_index_1);
        set_data(l_row_index, p_column_index_1, get_data(l_row_index, p_column_index_2));
        set_data(l_row_index, p_column_index_2, l_tmp_data);
    }
}

//-------------------------------------------------------------------------
template <typename T>
my_matrix<T>
my_matrix<T>::mult(const my_matrix & p_matrix)
{
    if(m_width != p_matrix.get_height())
    {
        throw quicky_exception::quicky_logic_exception("my_matrix.class: fatal error ! you try to multiplicate two matrix with incompatible sizes", __LINE__, __FILE__);
    }
    my_matrix l_result(m_height, p_matrix.get_width());
    T l_total;

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
template <typename T>
std::string
my_matrix<T>::to_string() const
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
template <typename T>
my_matrix<T>::~my_matrix()
{
    delete[] m_data;
}

//-------------------------------------------------------------------------
template <typename T>
my_matrix<T>::my_matrix(const my_matrix & p_matrix)
:m_width(p_matrix.m_width)
,m_height(p_matrix.m_height)
,m_data(new T[m_width * m_height])
{
    memcpy(m_data, p_matrix.m_data, sizeof(T) * m_width * m_height);
}

//-------------------------------------------------------------------------
template <typename T>
bool
my_matrix<T>::operator==(const my_matrix & p_matrix) const
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
template <typename T>
my_matrix<T>::my_matrix(my_matrix && p_matrix)
:m_width(p_matrix.m_width)
,m_height(p_matrix.m_height)
,m_data(p_matrix.m_data)
{
    p_matrix.m_width = 0;
    p_matrix.m_height = 0;
    p_matrix.m_data = NULL;
}

//-----------------------------------------------------------------------------
template <typename T>
my_matrix<T> &
my_matrix<T>::operator=(const my_matrix & p_matrix)
{
    m_width = p_matrix.m_width;
    m_height = p_matrix.m_height;
    delete[] m_data;
    m_data = new T[m_width * m_height];
    memcpy(m_data, p_matrix.m_data, m_width * m_height * sizeof(T));
    return *this;
}

//-----------------------------------------------------------------------------
template <typename T>
my_matrix<T> &
my_matrix<T>::operator=(my_matrix && p_matrix)
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
template <typename T>
my_matrix<T>::my_matrix()
:m_width(0)
,m_height(0)
,m_data(NULL)
{

}

#ifdef SIMPLEX_SELF_TEST
//-------------------------------------------------------------------------
bool test_my_matrix();
#endif // SIMPLEX_SELF_TEST

#endif // _MY_MATRIX_H_
// EOF
