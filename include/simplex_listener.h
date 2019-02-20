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
#ifndef _SIMPLEX_LISTENER_H_
#define _SIMPLEX_LISTENER_H_


#include "simplex_listener_target_if.h"
#include <iostream>

namespace simplex
{
  template <typename COEF_TYPE>
  class simplex_listener
  {
    public:
    inline simplex_listener(const simplex_listener_target_if<COEF_TYPE> & p_simplex,
			                std::ostream & p_ostream = std::cout
			               );
    inline void start_iteration(const unsigned int & p_nb_iteration);
    inline void new_input_var_event(const unsigned int & p_input_variable_index);
    inline void new_output_var_event(const unsigned int & p_input_variable_index);
    inline void new_Z0(COEF_TYPE p_z0);
    private:
    unsigned int m_nb_iteration;
    const simplex_listener_target_if<COEF_TYPE> & m_simplex;
    std::ostream & m_ostream;
    bool m_start;
  };

  //----------------------------------------------------------------------------
  template <typename COEF_TYPE>
  simplex_listener<COEF_TYPE>::simplex_listener(const simplex_listener_target_if<COEF_TYPE> & p_simplex,
							                    std::ostream & p_ostream
							                   ):
    m_nb_iteration(0),
    m_simplex(p_simplex),
    m_ostream(p_ostream),
    m_start(true)
  {
  }

  //----------------------------------------------------------------------------
  template <typename COEF_TYPE>
  void simplex_listener<COEF_TYPE>::start_iteration(const unsigned int & p_nb_iteration)
  {
    m_nb_iteration = p_nb_iteration;
    if(m_start)
      {
        m_simplex.display_array(m_ostream);
        m_start = false;
      }
  }
 
  //----------------------------------------------------------------------------
  template <typename COEF_TYPE>
  void simplex_listener<COEF_TYPE>::new_input_var_event(const unsigned int & p_input_variable_index
					     )
  {
      m_ostream << "Iteration[" << m_nb_iteration << "] : New input variable selected : " << p_input_variable_index << std::endl;
  }
  
  //----------------------------------------------------------------------------
  template <typename COEF_TYPE>
  void simplex_listener<COEF_TYPE>::new_output_var_event(const unsigned int & p_output_variable_index
					     )
  {
    m_ostream << "Iteration[" << m_nb_iteration << "] : New output variable selected : " << p_output_variable_index << std::endl;
  }

  //----------------------------------------------------------------------------
  template<typename COEF_TYPE>
  void simplex_listener<COEF_TYPE>::new_Z0(const COEF_TYPE p_z0)
  {
    m_ostream << "Iteration[" << m_nb_iteration << "] : New Z0 : " << p_z0 << std::endl;
    m_simplex.display_array(m_ostream);
    std::vector<COEF_TYPE> l_values = m_simplex.get_variable_values();
    for(unsigned int l_index = 0; l_index < l_values.size(); ++l_index)
    {
        m_ostream << "Iteration[" << m_nb_iteration << "] : Variable value[" << l_index << "] = " << l_values[l_index] << std::endl;
    }
  }

}
#endif // _SIMPLEX_LISTENER_H_
//EOF
