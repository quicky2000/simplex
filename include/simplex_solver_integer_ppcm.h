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

#ifndef SIMPLEX_SOLVER_INTEGER_PPCM_H
#define SIMPLEX_SOLVER_INTEGER_PPCM_H

#include "simplex_array.h"
#include "simplex_solver_integer_base.h"

namespace simplex
{
    template <typename COEF_TYPE, typename ARRAY_TYPE=simplex_array<COEF_TYPE>>
    class simplex_solver_integer_ppcm: public simplex_solver_integer_base<COEF_TYPE,ARRAY_TYPE>
    {
      public:
        simplex_solver_integer_ppcm() = delete;

        simplex_solver_integer_ppcm(unsigned int p_nb_variables,
                               unsigned int p_nb_inequations_lt,
                               unsigned int p_nb_equations,
                               unsigned int p_nb_inequations_gt
                              );

      private:
        /**
         * Method performing pivot to change the base
         * The A coefficient A[row,column] should be !0
         * @param p_row_index Row index
         * @param p_column_index Column index
         * */
        inline void pivot(const unsigned int p_row_index,
                          const unsigned int p_column_index
                         ) override ;

    };

    //-------------------------------------------------------------------------
    template <typename COEF_TYPE, typename ARRAY_TYPE>
    simplex_solver_integer_ppcm<COEF_TYPE,ARRAY_TYPE>::simplex_solver_integer_ppcm(unsigned int p_nb_variables,
                                                                                   unsigned int p_nb_inequations_lt,
                                                                                   unsigned int p_nb_equations,
                                                                                   unsigned int p_nb_inequations_gt
                                                                                  ):
            simplex_solver_integer_base<COEF_TYPE,ARRAY_TYPE>(p_nb_variables,
                                                              p_nb_inequations_lt,
                                                              p_nb_equations,
                                                              p_nb_inequations_gt
                                                             )
    {
        static_assert(std::is_integral<COEF_TYPE>::value, "Simplex solver acccept only integer types");
        static_assert(std::is_signed<COEF_TYPE>::value, "Simplex solver acccept only signed types");
    }

    //-------------------------------------------------------------------------
    template <typename COEF_TYPE, typename ARRAY_TYPE>
    void
    simplex_solver_integer_ppcm<COEF_TYPE,ARRAY_TYPE>::pivot(const unsigned int p_row_index,
                                                        const unsigned int p_column_index
                                                       )
    {
        assert(p_row_index < this->get_nb_total_equations());
        assert(p_column_index < this->get_nb_all_variables());
        COEF_TYPE l_pivot = this->get_internal_coef(p_row_index,p_column_index);
        assert(l_pivot);

        // Pivoting Z
        COEF_TYPE l_q = this->get_array().get_Z_coef(p_column_index);

        // Determine PPCM and pivoting to be able to set pivoting Z to 0
        COEF_TYPE l_ppcm = quicky_utils::fract<COEF_TYPE>::PPCM(l_q,l_pivot);
        COEF_TYPE l_mult_z = -l_ppcm / l_q;
        COEF_TYPE l_mult_array = -l_ppcm / l_pivot;

        std::cout << "Z line <= (Z * " << l_mult_z << ") - (R[" << p_row_index << "] * " << l_mult_array << ")" << std::endl;
        COEF_TYPE l_pgcd = 0;
        COEF_TYPE l_new_coef = 0;
        for(unsigned int l_index = 0;
            l_index < this->get_nb_all_variables();
            ++l_index
                )
        {
            COEF_TYPE l_u = this->get_internal_coef(p_row_index,l_index);
            l_new_coef = this->get_array().get_Z_coef(l_index) * l_mult_z - l_u * l_mult_array;
            this->accumulate_PGCD(l_pgcd, l_new_coef);
            this->get_array().set_Z_coef(l_index, l_new_coef);
        }
        assert(!this->get_array().get_Z_coef(p_column_index));
        l_new_coef = this->get_array().get_Z0_coef() * l_mult_z - this->get_array().get_B_coef(p_row_index) * l_mult_array;
        this->accumulate_PGCD(l_pgcd, l_new_coef);
        this->get_array().set_Z0_coef(l_new_coef);

        // Divide Z row by PGCD if necessary
        if(l_pgcd > 1)
        {
            std::cout << "Z line <= Z / " << l_pgcd << std::endl;
            for (unsigned int l_index = 0;
                 l_index < this->get_nb_all_variables();
                 ++l_index
                    )
            {
                this->get_array().set_Z_coef(l_index, this->get_array().get_Z_coef(l_index) / l_pgcd);
            }
            this->get_array().set_Z0_coef(this->get_array().get_Z0_coef() / l_pgcd);
        }
        // Pivoting other rows
        for(unsigned int l_row_index = 0;
            l_row_index < this->get_nb_total_equations();
            ++l_row_index
                )
        {
            if(l_row_index != p_row_index)
            {
                COEF_TYPE l_q = this->get_internal_coef(l_row_index,p_column_index);
                if(l_q)
                {
                    // Determine PPCM and pivoting to be able to set pivoting Z to 0
                    l_ppcm = quicky_utils::fract<COEF_TYPE>::PPCM(l_q,
                                                                  l_pivot
                                                                 );
                    if(l_q > 0)
                    {
                        l_mult_z = l_ppcm / l_q;
                        l_mult_array = l_ppcm / l_pivot;
                    } else
                    {
                        l_mult_z = -l_ppcm / l_q;
                        l_mult_array = -l_ppcm / l_pivot;
                    }

                    std::cout << "R[" << l_row_index << "] <= (R[" << l_row_index << "] * " << l_mult_z << ") - (R["
                              << p_row_index << "] * " << l_mult_array << ")" << std::endl;
                    l_pgcd = 0;
                    l_new_coef = this->get_array().get_B_coef(l_row_index) * l_mult_z - this->get_array().get_B_coef(p_row_index) * l_mult_array;
                    this->accumulate_PGCD(l_pgcd, l_new_coef);
                    this->get_array().set_B_coef(l_row_index, l_new_coef);
                    for (unsigned int l_index = 0;
                         l_index < this->get_nb_all_variables();
                         ++l_index
                        )
                    {
                        COEF_TYPE l_u = this->get_internal_coef(p_row_index, l_index);
                        l_new_coef = this->get_internal_coef(l_row_index, l_index) * l_mult_z - l_u * l_mult_array;
                        this->accumulate_PGCD(l_pgcd, l_new_coef);
                        this->set_internal_coef(l_row_index, l_index, l_new_coef);
                    }

                    if (l_pgcd > 1)
                    {
                        std::cout << "R[" << l_row_index << "] <= R[" << l_row_index << "] / " << l_pgcd << std::endl;
                        for (unsigned int l_index = 0;
                             l_index < this->get_nb_all_variables();
                             ++l_index
                            )
                        {
                            this->set_internal_coef(l_row_index, l_index, this->get_internal_coef(l_row_index, l_index) / l_pgcd);
                        }
                        this->get_array().set_B_coef(l_row_index, this->get_array().get_B_coef(l_row_index) / l_pgcd);
                    }
                }
            }
        }

        // Particular case of pivot row
        l_pgcd = 0;
        for(unsigned int l_index = 0;
            l_index < this->get_nb_variables();
            ++l_index
           )
        {
            this->accumulate_PGCD(l_pgcd,this->get_internal_coef(p_row_index, l_index));
        }
        this->accumulate_PGCD(l_pgcd,this->get_array().get_B_coef(p_row_index));
        if(l_pgcd > 1)
        {
            std::cout << "R[" << p_row_index << "] <= R[" << p_row_index << "] / " << l_pgcd << std::endl;
            for(unsigned int l_index = 0;
                l_index < this->get_nb_variables();
                ++l_index
               )
            {
                this->set_internal_coef(p_row_index,l_index,this->get_internal_coef(p_row_index,l_index) / l_pgcd);
            }
            this->get_array().set_B_coef(p_row_index, this->get_array().get_B_coef(p_row_index) / l_pgcd);
        }
    }
}


#endif //SIMPLEX_SOLVER_INTEGER_PPCM_H
