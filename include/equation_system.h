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

#ifndef _MY_EQUATION_SYSTEM_H_
#define _MY_EQUATION_SYSTEM_H_

#include "my_square_matrix.h"

template <typename T>
class SystemEquation
{
  private:
	my_square_matrix<T> m_matrix;
	my_matrix<T> m_coef;

  public:


    SystemEquation(const my_square_matrix<T> & p_matrix
                  ,const my_matrix<T> & p_coef
                  );

    my_matrix<T>
    solve();

    std::string to_string() const;
};

//-----------------------------------------------------------------------------
template <typename T>
SystemEquation<T>::SystemEquation(const my_square_matrix<T> & p_matrix
                                 ,const my_matrix<T> & p_coef
                                 )
:m_matrix(p_matrix)
,m_coef(p_coef)
{
    if(p_matrix.get_height() != p_coef.get_height())
    {
        throw quicky_exception::quicky_logic_exception("Equation p_matrix and coefficient p_matrix have incompatible sizes", __LINE__, __FILE__);
    }
}

//-----------------------------------------------------------------------------
template <typename T>
my_matrix<T>
SystemEquation<T>::solve()
{
    if(m_matrix.get_determ()==0)
    {
        return my_matrix<T>();
    }

    unsigned int l_width = m_matrix.get_width();
    bool l_pivot;
    std::tuple<T, unsigned int> l_max;
    my_matrix<T> l_result(l_width, 1);

    // Matrix triangularisation
    for(unsigned int i = 0; i < l_width ; ++i)
    {
        l_pivot = false;

        if(m_matrix.get_data(i,i) == 0)
        {
            l_pivot=true;
            l_max = m_matrix.max_abs_sub_column(i + 1,i);
            m_matrix.swap_line(i, std::get<1>(l_max));
            m_coef.swap_line(i, std::get<1>(l_max));
        }

        if(l_pivot == false || (l_pivot == true && m_matrix.get_data(i,i) != 0))
        {
            for(unsigned int i2 = i+1; i2 < l_width; ++i2)
            {
                T m = m_matrix.get_data(i2, i) / m_matrix.get_data(i,i);

                m_coef.set_data(i2, 0, m_coef.get_data(i2, 0) - m * m_coef.get_data(i,0));
                for(unsigned int j = 0; j < l_width; ++j)
                {
                    if(j == i)
                    {
                        m_matrix.set_data(i2, j, 0);
                    }
                    else
                    {
                        m_matrix.set_data(i2, j, m_matrix.get_data(i2, j) - m * m_matrix.get_data(i,j));
                    }
                }
            }
        }
    }

    // Triangular system resolution
    for(unsigned int j = l_width - 1; j < l_width; --j)
    {
        l_result.set_data(j, 0, m_coef.get_data(j, 0) / m_matrix.get_data(j,j));
        for(unsigned int i = j + 1; i < l_width; ++i)
        {
            l_result.set_data(j, 0, l_result.get_data(j, 0) - (l_result.get_data(i, 0) * m_matrix.get_data(j,i)) / m_matrix.get_data(j,j));
        }
    }
    return(l_result);
}

//-----------------------------------------------------------------------------
template <typename T>
std::string SystemEquation<T>::to_string() const
{
    std::string l_string("dimension");
    l_string += std::to_string(m_matrix.get_width()) +"\n";
    for(unsigned int i = 0; i  < m_matrix.get_height(); ++i)
    {
        for(unsigned int j = 0; j <m_matrix.get_width(); ++j)
        {
            l_string += std::to_string(m_matrix.get_data(i,j)) +"\t";
        }
        l_string += "=" + std::to_string(m_coef.get_data(i,0)) +"\n";
    }
    return(l_string);
}

#ifdef SIMPLEX_SELF_TEST
bool test_equation_system();
#endif // SIMPLEX_SELF_TEST

#endif // _MY_EQUATION_SYSTEM_H_
// EOF
