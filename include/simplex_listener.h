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

#include <iostream>

namespace simplex
{
  class simplex_listener
  {
    public:
    inline simplex_listener(void);
    inline void start_iteration(const unsigned int & p_nb_iteration);
    inline void new_input_var_event(const unsigned int & p_input_variable_index);
    inline void new_output_var_event(const unsigned int & p_input_variable_index);
    template<typename COEF_TYPE>
    inline void new_Z0(const COEF_TYPE p_z0);
    private:
    unsigned int m_nb_iteration;
  };

  //----------------------------------------------------------------------------
  simplex_listener::simplex_listener(void):
    m_nb_iteration(0)
  {
  }

  //----------------------------------------------------------------------------
  void simplex_listener::start_iteration(const unsigned int & p_nb_iteration)
  {
    m_nb_iteration = p_nb_iteration;
  }
 
  //----------------------------------------------------------------------------
  void simplex_listener::new_input_var_event(const unsigned int & p_input_variable_index
					     )
  {
    std::cout << "Iteration[" << m_nb_iteration << "] : New input variable selected : " << p_input_variable_index << std::endl;
  }
  
  //----------------------------------------------------------------------------
  void simplex_listener::new_output_var_event(const unsigned int & p_output_variable_index
					     )
  {
    std::cout << "Iteration[" << m_nb_iteration << "] : New output variable selected : " << p_output_variable_index << std::endl;
  }

  //----------------------------------------------------------------------------
  template<typename COEF_TYPE>
  void simplex_listener::new_Z0(const COEF_TYPE p_z0)
  {
    std::cout << "Iteration[" << m_nb_iteration << "] : New Z0 : " << p_z0 << std::endl;
  }
  
}
#endif // _SIMPLEX_LISTENER_H_
//EOF
