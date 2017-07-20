/*    This file is part of simplex
      Copyright (C) 2017  Julien Thevenon ( julien_thevenon at yahoo.fr )

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
#ifndef _SIMPLEX_ARRAY_H_
#define _SIMPLEX_ARRAY_H_

#include "simplex_array_base.h"
#include <cstring>
#include <cassert>

namespace simplex
{
  template <typename COEF_TYPE>
  class simplex_array: public simplex_array_base<COEF_TYPE>
  {
  public:
    simplex_array(const unsigned int & p_nb_equations,
		  const unsigned int & p_nb_variables
		  );

    /**
       Define coefficient for objective function
       @param p_index : the value should be less than number of variables
       @param value : value of coefficient in the formula Z = SUM(Cj * x)
    */
    inline void set_Z_coef(const unsigned int p_index,
			   const COEF_TYPE & p_value
			   );

    /**
       Define coefficient for objective function
       @param p_index : the value should be less than number of variables
       @param value : value of coefficient in the formula Z = SUM(Cj * x)
    */
    inline
      const COEF_TYPE &
      get_Z_coef(const unsigned int p_index
		 )const;

    /**
       Define coefficient Z0 for objective function
       @param value : value of coefficient in the formula Z = SUM(Cj * x)
    */
    inline void set_Z0_coef(const COEF_TYPE & p_value
		       );

    /**
       Return coefficient Z0 for objective function
       @param value : value of coefficient in the formula Z = SUM(Cj * x)
    */
    inline
      const COEF_TYPE &
      get_Z0_coef(void)const;

    /**
       Define coefficient for B coefficients in A x = b
       @param p_index : the value should be less than total number of equations
       @param value : value of coefficient in b
    */
    inline void set_B_coef(const unsigned int p_index,
			   const COEF_TYPE & p_value
			   );

    /**
       Return coefficient for B coefficients in A x = b
       @param p_index : the value should be less than total number of equations
       @param value : value of coefficient in b
    */
    inline const COEF_TYPE & get_B_coef(const unsigned int p_index
					)const;

    /**
       Define coefficient for A coefficients in A x = b
       @param p_equation_index : the value should be less than total number of equations
       @param p_variable_index : the value should be less than number of variables
       @param value : value of coefficient in A
    */
    inline void set_A_coef(const unsigned int p_equation_index,
			   const unsigned int p_variable_index,
			   const COEF_TYPE & p_value
			   );

    /**
       Return coefficient for A coefficients in A x = b
       @param p_equation_index : the value should be less than total number of equations
       @param p_variable_index : the value should be less than number of variables
       @return value : value of coefficient in A
    */
    inline
      const COEF_TYPE &
      get_A_coef(const unsigned int p_equation_index,
		 const unsigned int p_variable_index
		 ) const;

    inline ~simplex_array(void);

  private:
    /**
       Matrix A coefficients in Ax = b
     */
    COEF_TYPE * m_equation_coefs;

    /**
       b coefficients in Ax = b
     */
    COEF_TYPE * m_b_coefs;

    /**
       Coefficients of objective function in the form Z - SUM(Cj * X) = 0
    */
    COEF_TYPE * m_z_coefs;

    /**
       Z0 coefficient
    */
    COEF_TYPE m_z0;
  };

  //----------------------------------------------------------------------------
  template <typename COEF_TYPE>
  simplex_array<COEF_TYPE>::simplex_array(const unsigned int & p_nb_equations,
					  const unsigned int & p_nb_variables
					  ):
    simplex_array_base<COEF_TYPE>(p_nb_equations,p_nb_variables),
    m_equation_coefs(new COEF_TYPE[p_nb_variables * p_nb_equations]),
    m_b_coefs(new COEF_TYPE[p_nb_equations]),
    m_z_coefs(new COEF_TYPE[p_nb_variables]),
    m_z0(0)
    {
      memset(m_equation_coefs,0,p_nb_variables * p_nb_equations * sizeof(COEF_TYPE));
      memset(m_b_coefs,0,p_nb_equations * sizeof(COEF_TYPE));
      memset(m_z_coefs,0,p_nb_variables * sizeof(COEF_TYPE));
    }

  //----------------------------------------------------------------------------
  template <typename COEF_TYPE>
  void simplex_array<COEF_TYPE>::set_Z_coef(const unsigned int p_index,
				      const COEF_TYPE & p_value
				      )
    {
      assert(m_z_coefs);
      assert(p_index < simplex_array_base<COEF_TYPE>::get_nb_variables());
      m_z_coefs[p_index] = p_value;
    }

  //----------------------------------------------------------------------------
  template <typename COEF_TYPE>
  const COEF_TYPE & simplex_array<COEF_TYPE>::get_Z_coef(const unsigned int p_index
							 )const
    {
      assert(m_z_coefs);
      assert(p_index < simplex_array_base<COEF_TYPE>::get_nb_variables());
      return m_z_coefs[p_index];
    }

  //----------------------------------------------------------------------------
  template <typename COEF_TYPE>
  void simplex_array<COEF_TYPE>::set_Z0_coef(const COEF_TYPE & p_value
					)
    {
      m_z0 = p_value;
    }

  //----------------------------------------------------------------------------
  template <typename COEF_TYPE>
  const COEF_TYPE & simplex_array<COEF_TYPE>::get_Z0_coef(void)const
    {
      return m_z0;
    }

  //----------------------------------------------------------------------------
  template <typename COEF_TYPE>
    void simplex_array<COEF_TYPE>::set_B_coef(const unsigned int p_index,
					      const COEF_TYPE & p_value
					      )
    {
      assert(m_b_coefs);
      assert(p_index < simplex_array_base<COEF_TYPE>::get_nb_equations());
      m_b_coefs[p_index] = p_value;
    }

  //----------------------------------------------------------------------------
  template <typename COEF_TYPE>
    const COEF_TYPE & simplex_array<COEF_TYPE>::get_B_coef(const unsigned int p_index
							   )const
    {
      assert(m_b_coefs);
      assert(p_index < simplex_array_base<COEF_TYPE>::get_nb_equations());
      return m_b_coefs[p_index];
    }

  //----------------------------------------------------------------------------
  template <typename COEF_TYPE>
  void simplex_array<COEF_TYPE>::set_A_coef(const unsigned int p_equation_index,
					    const unsigned int p_variable_index,
					    const COEF_TYPE & p_value
					    )
  {
    assert(m_equation_coefs);
    assert(p_equation_index < simplex_array_base<COEF_TYPE>::get_nb_equations());
    assert(p_variable_index < simplex_array_base<COEF_TYPE>::get_nb_variables());
    m_equation_coefs[p_equation_index * simplex_array_base<COEF_TYPE>::get_nb_variables() + p_variable_index] = p_value;
  }

  //----------------------------------------------------------------------------
  template <typename COEF_TYPE>
  const COEF_TYPE &
  simplex_array<COEF_TYPE>::get_A_coef(const unsigned int p_equation_index,
				       const unsigned int p_variable_index
				       )const
    {
      assert(m_equation_coefs);
      assert(p_equation_index < simplex_array_base<COEF_TYPE>::get_nb_equations());
      assert(p_variable_index < simplex_array_base<COEF_TYPE>::get_nb_variables());
      return m_equation_coefs[p_equation_index * simplex_array_base<COEF_TYPE>::get_nb_variables() + p_variable_index];
    }

  //----------------------------------------------------------------------------
  template <typename COEF_TYPE>
  simplex_array<COEF_TYPE>::~simplex_array(void)
    {
      delete[] m_z_coefs;
      m_z_coefs = nullptr;
      delete[] m_b_coefs;
      m_b_coefs = nullptr;
      delete[] m_equation_coefs;
      m_equation_coefs = nullptr;
    }
}

#endif // _SIMPLEX_ARRAY_H_
// EOF
