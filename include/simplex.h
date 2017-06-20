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
#ifndef _SIMPLEX_H_
#define _SIMPLEX_H_

#include <quicky_exception.h>
#include <cstring>
#include <type_traits>
#include <cassert>
#include <iostream>
#include <sstream>
#include <string>
#include <limits>
		
namespace simplex
{
  typedef enum class equation_type
    {
      UNDEFINED = 0,
      EQUATION,
      INEQUATION_LT,
      INEQUATION_GT
    } t_equation_type;


  template <typename COEF_TYPE>
  class simplex
  {
  public:

    /**
       This constructor prepare a simplex in Canonical form
       Max z = cx
       a1 * x1 + a2 * x2 + ... + an * xn <= bn
       xi >= 0 for all xi

       Or

       Min z = cx
       a1 * x1 + a2 * x2 + ... + an * xn >= bn
       xi >= 0 for all xi
    */
    inline simplex(unsigned int p_nb_variables,
		   unsigned int p_nb_inequations_lt,
		   unsigned int p_nb_equations,
		   unsigned int p_nb_inequations_gt
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
       Define coefficient for B coefficients in A x = b
       @param p_index : the value should be less than total number of equations
       @param value : value of coefficient in b
    */
    inline void set_B_coef(const unsigned int p_index,
			   const COEF_TYPE & p_value
			   );

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

    inline ~simplex(void);

    /**
       Display simplex array representation
       @param p_stream stream where the display should be done
       @return the modified stream
     */
    inline std::ostream & display_array(std::ostream & p_stream);

    /**
       Define equation type
       @param Equation index
       @param Equation type
     */
    inline void define_equation_type(const unsigned int & p_equation_index,
				     const t_equation_type & p_equation_type
				     );
    /**
       Method performing pivot to change the base
       The A coefficient A[row,column] should be !0
       @param p_row_index Row index
       @param p_column_index Column index
     */
    inline void pivot(const unsigned int p_row_index,
		      const unsigned int p_column_index
		      );

    inline bool find_max(COEF_TYPE & p_max, bool & p_infinite);

  private:
    /**
       Define coefficient for A coefficients in A x = b
       @param p_equation_index : the value should be less than total number of equations
       @param p_variable_index : the value should be less than total number of variables
       including adjustment variables
       @param value : value of coefficient in A
    */
    inline void set_internal_coef(const unsigned int p_equation_index,
				  const unsigned int p_variable_index,
				  const COEF_TYPE & p_value
				  );

    /**
       Return coefficient for A coefficients in A x = b
       @param p_equation_index : the value should be less than total number of equations
       @param p_variable_index : the value should be less than total number of variables
       including adjustment variables
       @return value : value of coefficient in A
    */
    inline
      const COEF_TYPE &
      get_internal_coef(const unsigned int p_equation_index,
			const unsigned int p_variable_index
			) const;


    simplex(void) = delete;

    /**
       Method to set adjustement variable in the right place
       @param index of inequation
       @param value to set as coefficient
    */
    inline void set_adjustement_variable(const unsigned int & p_equation_index,
					 const COEF_TYPE & p_value
					 );
    /**
       Method to determine the next input variable for pivot operation
       when searching max optimum
       @param reference on variable where to store the input variable index if any
       @return boolean indicating if an input variable was found: TRUE = found
     */
    inline bool get_max_input_variable_index(unsigned int & p_variable_index) const;

    /**
       Method to determine the next output variable for pivot operation
       @param index of input variable
       @param reference on variable where to store the input variable index if any
       @return boolean indicating if an input variable was found
     */
    inline bool get_output_variable_index(unsigned int p_input_variable_index,
					      unsigned int & p_variable_index
					      )const;

    /**
       Method to determine the next input variable for pivot operation
       when searching min optimum
       @param reference on variable where to store the input variable index if any
       @return boolean indicating if an input variable was found: TRUE = found
     */
    inline bool get_min_input_variable_index(unsigned int & p_variable_index) const;

    /**
       Variable numbers without adjustment variables
    */
    unsigned int m_nb_variables;

    /**
       Number of inequations of type x1 + x2 + ... + xn <= value
     */
    unsigned int m_nb_inequations_lt;

    /**
       Number of equations of type x1 + x2 + ... + xn = value
     */
    unsigned int m_nb_equations;

    /**
       Number of inequations of type x1 + x2 + ... + xn >= value
     */
    unsigned int m_nb_inequations_gt;

    /**
       Variable counting the number of base variables specified
     */
    unsigned int m_nb_base_variables;

    /**
       Variable counting the number of adjustement variable specified
     */
    unsigned int m_nb_defined_adjustment_variables;

    /**
       Total number of variables including adjustment variables
    */
    unsigned int m_nb_all_variables;

    /**
       Total numer of equations after adjuments variables have been added
    */
    unsigned int m_nb_total_equations;

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

    /**
       Equation types
     */
    t_equation_type * m_equation_types;
  };

  //----------------------------------------------------------------------------
  std::ostream & operator<<(std::ostream & p_stream,
			    const t_equation_type & p_equation_type
			    )
    {
      switch(p_equation_type)
	{
	case t_equation_type::UNDEFINED:
	  p_stream << "\"undef\"";
	  break;
	case t_equation_type::EQUATION:
	  p_stream << "\"=\"";
	  break;
	case t_equation_type::INEQUATION_LT:
	  p_stream << "\"<=\"";
	break;
	case t_equation_type::INEQUATION_GT:
	  p_stream << "\">=\"";
	break;
	default:
	  throw quicky_exception::quicky_logic_exception("Unknown equation_type value : "+ std::to_string((unsigned int)p_equation_type),__LINE__,__FILE__);
	}
      return p_stream;
    }

  //----------------------------------------------------------------------------
  template <typename COEF_TYPE>
  simplex<COEF_TYPE>::simplex(unsigned int p_nb_variables,
			      unsigned int p_nb_inequations_lt,
			      unsigned int p_nb_equations,
			      unsigned int p_nb_inequations_gt
		   ):
    m_nb_variables(p_nb_variables),
    m_nb_inequations_lt(p_nb_inequations_lt),
    m_nb_equations(p_nb_equations),
    m_nb_inequations_gt(p_nb_inequations_gt),
    m_nb_base_variables(p_nb_inequations_lt + p_nb_inequations_gt),
    m_nb_defined_adjustment_variables(0),
    m_nb_all_variables(p_nb_variables + m_nb_base_variables),
    m_nb_total_equations(p_nb_inequations_lt + p_nb_equations + p_nb_inequations_gt),
    m_equation_coefs(new COEF_TYPE[m_nb_all_variables * m_nb_total_equations]),
    m_b_coefs(new COEF_TYPE[m_nb_total_equations]),
    m_z_coefs(new COEF_TYPE[m_nb_all_variables]),
    m_z0(0),
    m_equation_types(new t_equation_type[m_nb_total_equations])
      {
	static_assert(std::is_signed<COEF_TYPE>::value,"Simplex template parameter shoudl be signed");
	memset(m_equation_coefs,0,m_nb_all_variables * m_nb_total_equations * sizeof(COEF_TYPE));
	memset(m_b_coefs,0,m_nb_total_equations * sizeof(COEF_TYPE));
	memset(m_z_coefs,0,m_nb_all_variables * sizeof(COEF_TYPE));
	memset(m_equation_types, 0, m_nb_total_equations * sizeof(t_equation_type));
      }

  //----------------------------------------------------------------------------
  template <typename COEF_TYPE>
  void simplex<COEF_TYPE>::set_Z_coef(const unsigned int p_index,
				      const COEF_TYPE & p_value
				      )
    {
      assert(m_z_coefs);
      assert(p_index < m_nb_variables);
      // Negate version of value due to the change of Z form between function call and storage variable
      m_z_coefs[p_index] = - p_value;
    }

  //----------------------------------------------------------------------------
  template <typename COEF_TYPE>
  void simplex<COEF_TYPE>::set_B_coef(const unsigned int p_index,
				      const COEF_TYPE & p_value
				      )
    {
      assert(m_b_coefs);
      assert(p_index < m_nb_total_equations);
      m_b_coefs[p_index] = p_value;
    }

  //----------------------------------------------------------------------------
  template <typename COEF_TYPE>
  void simplex<COEF_TYPE>::set_A_coef(const unsigned int p_equation_index,
				      const unsigned int p_variable_index,
				      const COEF_TYPE & p_value
				      )
  {
    assert(m_equation_coefs);
    assert(p_equation_index < m_nb_total_equations);
    assert(p_variable_index < m_nb_variables);
    m_equation_coefs[p_equation_index * m_nb_all_variables + p_variable_index] = p_value;
  }

  //----------------------------------------------------------------------------
  template <typename COEF_TYPE>
    const COEF_TYPE &
    simplex<COEF_TYPE>::get_A_coef(const unsigned int p_equation_index,
				   const unsigned int p_variable_index
				   )const
    {
      assert(m_equation_coefs);
      assert(p_equation_index < m_nb_total_equations);
      assert(p_variable_index < m_nb_variables);
      return m_equation_coefs[p_equation_index * m_nb_all_variables + p_variable_index];
    }

  //----------------------------------------------------------------------------
  template <typename COEF_TYPE>
  void simplex<COEF_TYPE>::set_internal_coef(const unsigned int p_equation_index,
					     const unsigned int p_variable_index,
					     const COEF_TYPE & p_value
					     )
  {
    assert(m_equation_coefs);
    assert(p_equation_index < m_nb_total_equations);
    assert(p_variable_index < m_nb_all_variables);
    m_equation_coefs[p_equation_index * m_nb_all_variables + p_variable_index] = p_value;
  }

  //----------------------------------------------------------------------------
  template <typename COEF_TYPE>
    const COEF_TYPE &
    simplex<COEF_TYPE>::get_internal_coef(const unsigned int p_equation_index,
					  const unsigned int p_variable_index
					  )const
    {
      assert(m_equation_coefs);
      assert(p_equation_index < m_nb_total_equations);
      assert(p_variable_index < m_nb_all_variables);
      return m_equation_coefs[p_equation_index * m_nb_all_variables + p_variable_index];
    }

  //----------------------------------------------------------------------------
  template <typename COEF_TYPE>
  void simplex<COEF_TYPE>::pivot(const unsigned int p_row_index,
				 const unsigned int p_column_index
				 )
  {
    assert(p_row_index < m_nb_total_equations);
    assert(p_column_index < m_nb_all_variables);
    COEF_TYPE l_pivot = get_internal_coef(p_row_index,p_column_index);
    assert(l_pivot);

    // Pivoting Z
     COEF_TYPE l_q = m_z_coefs[p_column_index];
     for(unsigned int l_index = 0;
	 l_index < m_nb_all_variables;
	 ++l_index
	 )
       {
	 COEF_TYPE l_u = get_internal_coef(p_row_index,l_index);
	 m_z_coefs[l_index] = m_z_coefs[l_index] - (l_q * l_u) / l_pivot;
       }
     m_z0 = m_z0 - (l_q * m_b_coefs[p_row_index]) / l_pivot;

    // Pivoting other rows
    for(unsigned int l_row_index = 0;
	l_row_index < m_nb_total_equations;
	++l_row_index
	)
      {
	if(l_row_index != p_row_index)
	  {
	    COEF_TYPE l_q = get_internal_coef(l_row_index,p_column_index);
	    m_b_coefs[l_row_index] = m_b_coefs[l_row_index] - (l_q * m_b_coefs[p_row_index]) / l_pivot;
	    for(unsigned int l_index = 0;
		l_index < m_nb_all_variables;
		++l_index
		)
	      {
		COEF_TYPE l_u = get_internal_coef(p_row_index,l_index);
		set_internal_coef(l_row_index, l_index, get_internal_coef(l_row_index,l_index) - (l_q * l_u) / l_pivot);
	      }
	  }
      }

    // Particular case of pivot row
    for(unsigned int l_index = 0;
	l_index < m_nb_variables;
	++l_index
	)
      {
	set_internal_coef(p_row_index,l_index,get_internal_coef(p_row_index,l_index) / l_pivot);
      }
    m_b_coefs[p_row_index] = m_b_coefs[p_row_index] / l_pivot;
  }

  //----------------------------------------------------------------------------
  template <typename COEF_TYPE>
  std::ostream & simplex<COEF_TYPE>::display_array(std::ostream & p_stream)
    {
      p_stream << "Z\t";
     for(unsigned int l_index = 0;
	 l_index < m_nb_all_variables;
	 ++l_index
	 )
       {
	 p_stream << m_z_coefs[l_index] << "\t";
       }
     p_stream << "|\t" << m_z0 << std::endl;
     for(unsigned int l_row_index = 0;
	 l_row_index < m_nb_total_equations;
	 ++l_row_index
	 )
       {
	 p_stream << " \t";
	 for(unsigned int l_index = 0;
	     l_index < m_nb_all_variables;
	     ++l_index
	     )
	   {
	     p_stream << get_internal_coef(l_row_index, l_index) << "\t";
	   }
	 p_stream << "|\t" << m_b_coefs[l_row_index] << std::endl;
      }
     return p_stream;
    }

  //----------------------------------------------------------------------------
  template <typename COEF_TYPE>
  void simplex<COEF_TYPE>::set_adjustement_variable(const unsigned int & p_equation_index,
						    const COEF_TYPE & p_value
						    )
    {
      unsigned int l_column_index = m_nb_variables + m_nb_defined_adjustment_variables;
      set_internal_coef(p_equation_index, l_column_index, p_value);
      ++m_nb_defined_adjustment_variables;
    }

  //----------------------------------------------------------------------------
  template <typename COEF_TYPE>
  bool simplex<COEF_TYPE>::get_max_input_variable_index(unsigned int & p_variable_index) const
    {
      for(unsigned int l_index = 0;
	  l_index < m_nb_all_variables;
	  ++l_index
	  )
	{
	  if(m_z_coefs[l_index] < 0)
	    {
	      p_variable_index = l_index;
	      return true;
	    }
	}
      return false;
    }

  //----------------------------------------------------------------------------
  template <typename COEF_TYPE>
  bool simplex<COEF_TYPE>::get_min_input_variable_index(unsigned int & p_variable_index) const
    {
      for(unsigned int l_index = 0;
	  l_index < m_nb_all_variables;
	  ++l_index
	  )
	{
	  if(m_z_coefs[l_index] > 0)
	    {
	      p_variable_index = l_index;
	      return true;
	    }
	}
      return false;
    }

  //----------------------------------------------------------------------------
  template <typename COEF_TYPE>
  bool simplex<COEF_TYPE>::get_output_variable_index(unsigned int p_input_variable_index,
							 unsigned int & p_variable_index
							 )const
    {
      assert(p_input_variable_index < m_nb_all_variables);
      bool l_found = false;
      COEF_TYPE l_min = std::numeric_limits<COEF_TYPE>::max();
      for(unsigned int l_index = 0;
	  l_index < m_nb_total_equations;
	  ++l_index
	  )
	{
	  COEF_TYPE l_divider = get_internal_coef(l_index,p_input_variable_index);
	  if(l_divider > 0)
	    {
	      COEF_TYPE l_result = m_b_coefs[l_index] / l_divider;
	      if(l_result < l_min)
		{
		  l_min = l_result;
		  p_variable_index = l_index;
		  l_found = true;
		}
	    }
	}
      return l_found;
    }

  //----------------------------------------------------------------------------
  template <typename COEF_TYPE>
  void simplex<COEF_TYPE>::define_equation_type(const unsigned int & p_equation_index,
						const t_equation_type & p_equation_type
						)
    {
      assert(p_equation_index < m_nb_total_equations);
      m_equation_types[p_equation_index] = p_equation_type;
      switch(p_equation_type)
	{
	case t_equation_type::UNDEFINED:
	  throw quicky_exception::quicky_logic_exception("Try to set undefined equation type for equation" + std::to_string(p_equation_index),__LINE__,__FILE__);
	  break;
	case t_equation_type::EQUATION:
	  // Nothing to do
	  break;
	case t_equation_type::INEQUATION_LT:
	  set_adjustement_variable(p_equation_index,1);
	  break;
	case t_equation_type::INEQUATION_GT:
	  set_adjustement_variable(p_equation_index,-1);
	  break;
	default:
	  throw quicky_exception::quicky_logic_exception("Unknown equation_type value : "+ std::to_string((unsigned int)p_equation_type),__LINE__,__FILE__);
	}
    }

  //----------------------------------------------------------------------------
  template <typename COEF_TYPE>
  bool simplex<COEF_TYPE>::find_max(COEF_TYPE & p_max, bool & p_infinite)
    {
      std::cout << "---------------------------------" << std::endl;
      display_array(std::cout);
      p_infinite = false;
      bool l_input_found = true;
      unsigned int l_input_variable_index = 0;
      while(true == (l_input_found = get_max_input_variable_index(l_input_variable_index)))
	{
	  unsigned int l_ouput_variable_index = 0;
	  if(get_output_variable_index(l_input_variable_index,l_ouput_variable_index))
	    {
	      pivot(l_ouput_variable_index,l_input_variable_index);
	      std::cout << "---------------------------------" << std::endl;
	      display_array(std::cout);
	    }
	  else
	    {
	      p_infinite = true;
	      return false;
	    }
	}
      p_max = m_z0;
      return true;
    }

  //----------------------------------------------------------------------------
  template <typename COEF_TYPE>
  simplex<COEF_TYPE>::~simplex(void)
    {
      delete[] m_equation_types;
      m_equation_types = nullptr;
      delete[] m_z_coefs;
      m_z_coefs = nullptr;
      delete[] m_b_coefs;
      m_b_coefs = nullptr;
      delete[] m_equation_coefs;
      m_equation_coefs = nullptr;
    }
}

#endif // _SIMPLEX_H_
// EOF
