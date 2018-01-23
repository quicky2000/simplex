/*    This file is part of simplex
      Copyright (C) 2018  Julien Thevenon ( julien_thevenon at yahoo.fr )

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

#ifndef SIMPLEX_SOLVER_INTEGER_BASE_H
#define SIMPLEX_SOLVER_INTEGER_BASE_H

#include "simplex_solver_base.h"

namespace simplex
{
    template <typename COEF_TYPE,typename ARRAY_TYPE>
    class simplex_solver_integer_base: public simplex_solver_base<COEF_TYPE, ARRAY_TYPE>
    {

      public:
        simplex_solver_integer_base() = delete;

        simplex_solver_integer_base(unsigned int p_nb_variables,
                                    unsigned int p_nb_inequations_lt,
                                    unsigned int p_nb_equations,
                                    unsigned int p_nb_inequations_gt
                                   );

        ~simplex_solver_integer_base();

        void set_Z_coef(const unsigned int p_index,
                        const COEF_TYPE & p_value
                       );

        /**
                 Method implementing simplex algorithm to find max optimum solution
                 The problem must be in solved form
                 @param reference on variable where result will be stored
                 @param reference on a boolean value that will receive true if max
                 is infinite
                 @return value indicating if a max was found
                 */
        template <class LISTENER=simplex_listener<COEF_TYPE,ARRAY_TYPE>>
        bool
        find_max(COEF_TYPE & p_max,
                 bool & p_infinite,
                 LISTENER * p_listener = NULL
                );

      private:
        COEF_TYPE * m_original_Z_coefs;
    };

    //-------------------------------------------------------------------------
    template <typename COEF_TYPE, typename ARRAY_TYPE>
    simplex::simplex_solver_integer_base<COEF_TYPE, ARRAY_TYPE>::simplex_solver_integer_base(unsigned int p_nb_variables,
                                                                                             unsigned int p_nb_inequations_lt,
                                                                                             unsigned int p_nb_equations,
                                                                                             unsigned int p_nb_inequations_gt
                                                                                            ):
            simplex_solver_base<COEF_TYPE,ARRAY_TYPE>(p_nb_variables,
                                                      p_nb_inequations_lt,
                                                      p_nb_equations,
                                                      p_nb_inequations_gt
                                                     )
    {
        m_original_Z_coefs = new COEF_TYPE[this->get_nb_all_variables()];
        memset(m_original_Z_coefs, 0, this->get_nb_all_variables() * sizeof(COEF_TYPE));
        static_assert(std::is_integral<COEF_TYPE>::value, "Simplex solver acccept only integer types");
        static_assert(std::is_signed<COEF_TYPE>::value, "Simplex solver acccept only signed types");
    }

    //-------------------------------------------------------------------------
    template <typename COEF_TYPE, typename ARRAY_TYPE>
    simplex::simplex_solver_integer_base<COEF_TYPE, ARRAY_TYPE>::~simplex_solver_integer_base()
    {
        delete m_original_Z_coefs;
    }

    //-------------------------------------------------------------------------
    template <typename COEF_TYPE, typename ARRAY_TYPE>
    void
    simplex::simplex_solver_integer_base<COEF_TYPE, ARRAY_TYPE>::set_Z_coef(const unsigned int p_index,
                                                                            const COEF_TYPE & p_value
                                                                           )
    {
        simplex_solver_base<COEF_TYPE, ARRAY_TYPE>::set_Z_coef(p_index,
                                                               p_value
                                                              );
        m_original_Z_coefs[p_index] = -p_value;
    }

    //-------------------------------------------------------------------------
    template <typename COEF_TYPE, typename ARRAY_TYPE>
    template <class LISTENER>
    bool
    simplex::simplex_solver_integer_base<COEF_TYPE, ARRAY_TYPE>::find_max(COEF_TYPE & p_max,
                                                                          bool & p_infinite,
                                                                          LISTENER *p_listener
                                                                         )
    {
        bool l_result = simplex_solver_base<COEF_TYPE, ARRAY_TYPE>::find_max(p_max,
                                                                             p_infinite,
                                                                             p_listener
                                                                            );
        COEF_TYPE l_computed_max = 0;
        for(unsigned int l_index = 0;
            l_index < this->get_nb_variables();
            ++l_index
                )
        {
            std::cout << "Check variable " << l_index << std::endl;
            unsigned int l_row_index = this->get_base_variables_position(l_index);
            if(std::numeric_limits<unsigned int>::max() != l_row_index)
            {
                std::cout << "Corresponding row index : " << l_row_index << std::endl;
                COEF_TYPE l_coef = this->get_array().get_A_coef(l_row_index, l_index);
                std::cout << "Corresponding coef : " << l_coef << std::endl;
                COEF_TYPE l_B_coef = this->get_array().get_B_coef(l_row_index);
                std::cout << "Corresponding B coef : " << l_B_coef << std::endl;
                COEF_TYPE l_var_value = l_B_coef / l_coef;
                std::cout << "Variable value : " << l_var_value << std::endl;
                COEF_TYPE l_z_coef = m_original_Z_coefs[l_index];
                std::cout << "Corresponding Z coef : " << l_z_coef << std::endl;
                l_computed_max -= l_z_coef * l_var_value;
            }
        }
        p_max = l_computed_max;
        return l_result;
    }
}

#endif //SIMPLEX_SOLVER_INTEGER_BASE_H
// EOF
