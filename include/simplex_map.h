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
#ifndef _SIMPLEX_MAP_H_
#define _SIMPLEX_MAP_H_

#include "simplex_array_base.h"
#include <cassert>
#include <map>

namespace simplex
{
  template <typename COEF_TYPE>
    class simplex_map: public simplex_array_base<COEF_TYPE>
  {
  public:
    simplex_map(const unsigned int & p_nb_equations,
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

    inline ~simplex_map(void);

  private:
    inline
    void set_coef(const unsigned int & p_equation_index,
		  const unsigned int & p_variable_index,
		  const COEF_TYPE & p_value
		  );

    inline
    const COEF_TYPE & get_coef(const unsigned int & p_line_index,
			       const unsigned int & p_column_index
			       )const;

    std::map<unsigned int,COEF_TYPE> * m_coefs;
    COEF_TYPE m_zero;
  };


  //----------------------------------------------------------------------------
  template <typename COEF_TYPE>
  simplex_map<COEF_TYPE>::simplex_map(const unsigned int & p_nb_equations,
				       const unsigned int & p_nb_variables
				       ):
    simplex_array_base<COEF_TYPE>(p_nb_equations,p_nb_variables),
    m_coefs(new std::map<unsigned int,COEF_TYPE>[p_nb_equations + 1]),
    m_zero(0)
    {
    }


  //----------------------------------------------------------------------------
  template <typename COEF_TYPE>
  void simplex_map<COEF_TYPE>::set_coef(const unsigned int & p_equation_index,
					const unsigned int & p_variable_index,
					const COEF_TYPE & p_value
					)
  {
    assert(p_equation_index <= simplex_array_base<COEF_TYPE>::get_nb_equations());
    assert(m_coefs);
    typename std::map<unsigned int,COEF_TYPE> & l_line = m_coefs[p_equation_index];
    typename std::map<unsigned int,COEF_TYPE>::iterator l_iter = l_line.find(p_variable_index);
    if(l_iter == l_line.end())
      {
	if(p_value)
	  {
	    l_line.insert(typename std::map<unsigned int,COEF_TYPE>::value_type(p_variable_index,p_value));
	  }
      }
    else
      {
	if(p_value)
	  {
	    l_iter->second = p_value;
	  }
	else
	  {
	    l_line.erase(l_iter);
	  }
      }
  }

  //----------------------------------------------------------------------------
  template <typename COEF_TYPE>
  const COEF_TYPE & simplex_map<COEF_TYPE>::get_coef(const unsigned int & p_line_index,
						     const unsigned int & p_column_index
						     )const
  {
    assert(p_line_index <= simplex_array_base<COEF_TYPE>::get_nb_equations());
    assert(m_coefs);
    typename std::map<unsigned int,COEF_TYPE>::iterator l_iter = m_coefs[p_line_index].find(p_column_index);
    if(l_iter != m_coefs[p_line_index].end())
      {
	return l_iter->second;
      }
    return m_zero;
  }

  //----------------------------------------------------------------------------
  template <typename COEF_TYPE>
  void simplex_map<COEF_TYPE>::set_Z_coef(const unsigned int p_index,
					  const COEF_TYPE & p_value
					  )
    {
      assert(p_index < simplex_array_base<COEF_TYPE>::get_nb_variables());
      set_coef(simplex_array_base<COEF_TYPE>::get_nb_equations(),p_index,p_value);
    }

  //----------------------------------------------------------------------------
  template <typename COEF_TYPE>
  const COEF_TYPE & simplex_map<COEF_TYPE>::get_Z_coef(const unsigned int p_index
							 )const
    {
      assert(p_index < simplex_array_base<COEF_TYPE>::get_nb_variables());
      return get_coef(simplex_array_base<COEF_TYPE>::get_nb_equations(),p_index);
    }

  //----------------------------------------------------------------------------
  template <typename COEF_TYPE>
  void simplex_map<COEF_TYPE>::set_Z0_coef(const COEF_TYPE & p_value
					     )
    {
      set_coef(simplex_array_base<COEF_TYPE>::get_nb_equations(),simplex_array_base<COEF_TYPE>::get_nb_variables(),p_value);
    }

  //----------------------------------------------------------------------------
  template <typename COEF_TYPE>
  const COEF_TYPE & simplex_map<COEF_TYPE>::get_Z0_coef(void)const
    {
      return get_coef(simplex_array_base<COEF_TYPE>::get_nb_equations(),simplex_array_base<COEF_TYPE>::get_nb_variables());
    }

  //----------------------------------------------------------------------------
  template <typename COEF_TYPE>
    void simplex_map<COEF_TYPE>::set_B_coef(const unsigned int p_index,
					      const COEF_TYPE & p_value
					      )
    {
      assert(p_index < simplex_array_base<COEF_TYPE>::get_nb_equations());
      set_coef(p_index,simplex_array_base<COEF_TYPE>::get_nb_variables(),p_value);
    }

  //----------------------------------------------------------------------------
  template <typename COEF_TYPE>
    const COEF_TYPE & simplex_map<COEF_TYPE>::get_B_coef(const unsigned int p_index
							   )const
    {
      assert(p_index < simplex_array_base<COEF_TYPE>::get_nb_equations());
      return get_coef(p_index,simplex_array_base<COEF_TYPE>::get_nb_variables());
    }

  //----------------------------------------------------------------------------
  template <typename COEF_TYPE>
  void simplex_map<COEF_TYPE>::set_A_coef(const unsigned int p_equation_index,
					    const unsigned int p_variable_index,
					    const COEF_TYPE & p_value
					    )
  {
    assert(p_equation_index < simplex_array_base<COEF_TYPE>::get_nb_equations());
    assert(p_variable_index < simplex_array_base<COEF_TYPE>::get_nb_variables());
    set_coef(p_equation_index,p_variable_index, p_value);
  }

  //----------------------------------------------------------------------------
  template <typename COEF_TYPE>
  const COEF_TYPE &
  simplex_map<COEF_TYPE>::get_A_coef(const unsigned int p_equation_index,
				       const unsigned int p_variable_index
				       )const
    {
      assert(p_equation_index < simplex_array_base<COEF_TYPE>::get_nb_equations());
      assert(p_variable_index < simplex_array_base<COEF_TYPE>::get_nb_variables());
      return get_coef(p_equation_index, p_variable_index);
    }

  //----------------------------------------------------------------------------
  template <typename COEF_TYPE>
  simplex_map<COEF_TYPE>::~simplex_map(void)
    {
      delete[] m_coefs;
    }
}
#endif // _SIMPLEX_MAP_H_
// EOF
