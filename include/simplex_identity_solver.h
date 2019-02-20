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

#ifndef SIMPLEX_SIMPLEX_IDENTITY_SOLVER_H
#define SIMPLEX_SIMPLEX_IDENTITY_SOLVER_H

#include "simplex_listener.h"
#include "simplex_array.h"
#include <cassert>
#include <type_traits>
#include <string>
#include <iostream>

namespace simplex
{
    /**
     * Simplex solver to solve simplex problem having an identity matrix
     * This is the case when there are only inequations of form <
     * General form of simplex is Ax = b
     * Here we assume that A = (A' I) with I identity matrix
     * @tparam COEF_TYPE Type of array coef
     */
    template <typename COEF_TYPE>
    class simplex_identity_solver: public simplex_listener_target_if<COEF_TYPE>
    {
      public:
        /**
         * Constructor
         * @param p_nb_variables Number of variable
         * @param p_nb_equations Number of equations
         */
        simplex_identity_solver(unsigned int p_nb_variables
                               ,unsigned int p_nb_equations
                               );

        ~simplex_identity_solver();

        /**
         * Define coefficient for A' coefficients in A'I x = b
         * @param p_equation_index : the value should be less than number of equations
         * @param p_variable_index : the value should be less than number of variables
         * @param value : value of coefficient in A
         */
        void set_A_coef(const unsigned int p_equation_index
                       ,const unsigned int p_variable_index
                       ,const COEF_TYPE & p_value
                       );

        /**
         * Define coefficient for B coefficients in A x = b
         * @param p_equation_index : the value should be less than total number of equations
         * @param value : value of coefficient in b
         */
        void set_B_coef(const unsigned int p_equation_index
                       ,const COEF_TYPE & p_value
                       );

        /**
         * Define coefficient for objective function
         * @param p_variable_index : the value should be less than number of variables
         * @param value : value of coefficient in the formula Z = SUM(Cj * x)
         */
        void set_Z_coef(const unsigned int p_variable_index
                       ,const COEF_TYPE & p_value
                       );

        /**
         * Define coefficient Z0 for objective function
         * @param p_value : value of coefficient in the formula Z = SUM(Cj * x)
         */
        void set_Z0_coef(const COEF_TYPE & p_value);

        /**
         * Display simplex array representation
         * @param p_stream stream where the display should be done
         * @return the modified stream
         */
        std::ostream & display_array(std::ostream & p_stream)const override;

        /**
         * Return value of variables
         * @return value of variables
         */
        std::vector<COEF_TYPE> get_variable_values() const override;

        /**
         * Method implementing simplex algorithm to find max optimum solution
         * The problem must be in solved form
         * @param p_max reference on variable where result will be stored
         * @param p_infinite reference on a boolean value that will receive true if max
         * is infinite
         * @return value indicating if a max was found
         */
        template <class LISTENER=simplex_listener<COEF_TYPE>>
        bool
        find_max(COEF_TYPE & p_max
                ,bool & p_infinite
                ,LISTENER * p_listener = NULL
                );

      private:

        /**
         * Method performing pivot to change the base
         * The A' coefficient A'[row,column] should be !0
         * @param p_row_index Row index
         * @param p_column_index Column index
         */
        void pivot(const unsigned int p_row_index
                  ,const unsigned int p_column_index
                  ,const unsigned int p_base_column
                  );

        /**
         * Pivot Z coefficient at column index
         * @param p_pivot pivot coefficient
         * @param l_q
         * @param p_row_index
         * @param p_column_index
         */
        void pivot_Z(const COEF_TYPE & p_pivot
                    ,const COEF_TYPE & l_q
                    ,unsigned int p_row_index
                    ,unsigned int p_column_index
                    );

        /**
         * Pivot A coefficient at column index
         * @param p_pivot pivot coefficient
         * @param l_q
         * @param p_row_index
         * @param p_column_index
         */
        void
        pivot_A_rows(unsigned int p_start_row_index
                    ,unsigned int p_end_row_index
                    ,const COEF_TYPE & p_pivot
                    ,unsigned int p_pivot_row_index
                    ,unsigned int p_pivot_column_index
                    ,unsigned int p_base_column
                    );

        void
        pivot_A_coef(unsigned int p_start_column
                    ,unsigned int p_end_column
                    ,COEF_TYPE p_pivot
                    ,COEF_TYPE p_q
                    ,unsigned int p_row_index
                    ,unsigned int p_pivot_row_index
                    );

        /**
         * Return name of variable, if p_index >= variable numbers this is an
         * adjustment variable
         * @param p_index Variable index
         * @return variable name
         */
        const std::string & get_variable_name(unsigned int p_index) const;

        /**
         * Return variable index corresponding to base_variable index
         * @param p_index base variable index
         * @return variable index
         */
        unsigned int get_base_variable_index(unsigned int p_index) const;

        /**
         * Record variable index in base variable table
         * @param p_index index in base
         * @param p_variable_index variable index
         */
        void set_base_variable_index(unsigned int p_index
                                    ,unsigned int p_variable_index
                                    );
        /**
         * Return variable index corresponding to array column
         * @param p_column_index column index
         * @return variable index
         */
        unsigned int get_array_variable_index(unsigned int p_column_index) const;

        /**
         * Record variable index corresponding to A' column
         * @param p_index column index
         * @param p_variable_index variable index
         */
        void set_array_variable_index(unsigned int p_index
                                     ,unsigned int p_variable_index
                                     );

        /**
         * Return base variable index corresponding to equation/ array line
         * @param p_equation_index equation/line index
         * @return variable index
         */
        unsigned int get_equation_base_column_index(unsigned int p_equation_index) const;

        /**
         * Method to determine the next input variable for pivot operation
         * when searching max optimum
         * @param p_variable_index on variable where to store the input variable index if any
         * * @return boolean indicating if an input variable was found: TRUE = found
         */
        bool get_max_input_column_index(unsigned int & p_variable_index) const;

        /**
         * Method to determine the equation index corresponding to next output
         * variable for pivot operation
         * @param p_column_index index of input variable
         * @param p_equation_index reference on variable where to store the output equation index if any
         * @return boolean indicating if an input variable was found
         */
        bool
        get_output_equation_index(unsigned int p_column_index,
                                  unsigned int & p_equation_index
                                 ) const;

        /**
         * Total number of variables: declared by user + base variables used to
         * transform < inequations to equations
         */
        unsigned int m_nb_total_variables;

        /**
         * Coefficient array A' in ( Ax =b <=> (A'I x) = B
         */
        simplex_array<COEF_TYPE> m_array;

        /**
         Store variable index corresponding to each base variable
         */
         unsigned int * m_base_variables_index;

         /**
          * Store variable index corresponding to each A' column
          */
        unsigned int * m_array_variables_index;

        /**
         * Store base variable index corresponding to an equation
         */
        unsigned int * m_equation_base_variable_index;

        /**
         * Variable names
         */
         std::string * m_variable_names;
    };

    //-------------------------------------------------------------------------
    template <typename COEF_TYPE>
    simplex_identity_solver<COEF_TYPE>::simplex_identity_solver(unsigned int p_nb_variables
                                                               ,unsigned int p_nb_equations
                                                               )
                                                               :m_nb_total_variables(p_nb_variables + p_nb_equations)
                                                               ,m_array(p_nb_equations, p_nb_variables)
                                                               ,m_base_variables_index(new unsigned int[p_nb_equations])
                                                               ,m_array_variables_index(new unsigned int[p_nb_variables])
                                                               ,m_equation_base_variable_index(new unsigned int[p_nb_equations])
                                                               ,m_variable_names( new std::string[m_nb_total_variables])
    {
        static_assert(std::is_signed<COEF_TYPE>::value, "Simplex template parameter should be signed");

        // At the beginning array variables are normal variables
        for(unsigned int l_index = 0; l_index < p_nb_variables; ++l_index)
        {
            m_array_variables_index[l_index] = l_index;
            m_variable_names[l_index] = "X" + std::to_string(l_index);
        }

        // At the beginning base variables are adjustment variables
        for(unsigned int l_index = 0; l_index < p_nb_equations; ++l_index)
        {
            m_base_variables_index[l_index] = l_index + p_nb_variables;
            m_equation_base_variable_index[l_index] = l_index;
            m_variable_names[l_index + p_nb_variables] = "E" + std::to_string(l_index);
        }
    }

    //-------------------------------------------------------------------------
    template <typename COEF_TYPE>
    simplex_identity_solver<COEF_TYPE>::~simplex_identity_solver()
    {
        delete[] m_array_variables_index;
        delete[] m_base_variables_index;
        delete[] m_equation_base_variable_index;
        delete[] m_variable_names;
    }

    //-------------------------------------------------------------------------
    template <typename COEF_TYPE>
    void
    simplex_identity_solver<COEF_TYPE>::set_A_coef(const unsigned int p_equation_index
                                                  ,const unsigned int p_variable_index
                                                  ,const COEF_TYPE & p_value
                                                  )
    {
        m_array.set_A_coef(p_equation_index, p_variable_index, p_value);
    }

    //-------------------------------------------------------------------------
    template <typename COEF_TYPE>
    void
    simplex_identity_solver<COEF_TYPE>::set_B_coef(const unsigned int p_equation_index
                                                  ,const COEF_TYPE & p_value
                                                  )
    {
        m_array.set_B_coef(p_equation_index, p_value);
    }

    //-------------------------------------------------------------------------
    template <typename COEF_TYPE>
    void
    simplex_identity_solver<COEF_TYPE>::set_Z_coef(const unsigned int p_variable_index
                                                  ,const COEF_TYPE & p_value
                                                  )
    {
        m_array.set_Z_coef(p_variable_index, -p_value);
    }

    //-------------------------------------------------------------------------
    template <typename COEF_TYPE>
    void
    simplex_identity_solver<COEF_TYPE>::set_Z0_coef(const COEF_TYPE & p_value)
    {
        m_array.set_Z0_coef(p_value);
    }

    //-------------------------------------------------------------------------
    template <typename COEF_TYPE>
    std::ostream &
    simplex_identity_solver<COEF_TYPE>::display_array(std::ostream & p_stream) const
    {
        // Display variable names
        p_stream << "Z";
        for(unsigned int l_index = 0; l_index < m_array.get_nb_variables(); ++l_index)
        {
            p_stream << "\t" << get_variable_name(get_array_variable_index(l_index));
        }
        p_stream << std::endl;

        // Display Z coefs
        p_stream << "Z\t";
        for(unsigned int l_index = 0;
            l_index < m_array.get_nb_variables();
            ++l_index
           )
        {
            p_stream << m_array.get_Z_coef(l_index) << "\t";
        }
        p_stream << "|\t" << m_array.get_Z0_coef() << std::endl;

        // Display coef
        for(unsigned int l_row_index = 0;
            l_row_index < m_array.get_nb_equations();
            ++l_row_index
           )
        {
            p_stream << get_variable_name(get_base_variable_index(get_equation_base_column_index(l_row_index))) << "\t";
            for(unsigned int l_index = 0;
                l_index < m_array.get_nb_variables();
                ++l_index
               )
            {
                p_stream << m_array.get_A_coef(l_row_index, l_index) << "\t";
            }
            p_stream << "|\t" << m_array.get_B_coef(l_row_index) << std::endl;
        }
        p_stream << "Base variables:" << std::endl;
        for(unsigned int l_index = 0; l_index < m_array.get_nb_equations(); ++l_index)
        {
            p_stream << get_variable_name(get_base_variable_index(l_index)) << "\t";
        }
        p_stream << std::endl;
        return p_stream;
    }

    //-------------------------------------------------------------------------
    template <typename COEF_TYPE>
    const std::string &
    simplex_identity_solver<COEF_TYPE>::get_variable_name(unsigned int p_index) const
    {
        assert(p_index < m_nb_total_variables);
        return m_variable_names[p_index];
    }

    //-------------------------------------------------------------------------
    template <typename COEF_TYPE>
    unsigned int
    simplex_identity_solver<COEF_TYPE>::get_array_variable_index(unsigned int p_column_index) const
    {
        assert(p_column_index < m_array.get_nb_variables());
        return m_array_variables_index[p_column_index];
    }

    //-------------------------------------------------------------------------
    template <typename COEF_TYPE>
    void
    simplex_identity_solver<COEF_TYPE>::set_array_variable_index(unsigned int p_index
                                                                ,unsigned int p_variable_index
                                                                )
    {
        assert(p_index < m_array.get_nb_variables());
        m_array_variables_index[p_index] = p_variable_index;
    }

    //-------------------------------------------------------------------------
    template <typename COEF_TYPE>
    unsigned int
    simplex_identity_solver<COEF_TYPE>::get_equation_base_column_index(unsigned int p_equation_index) const
    {
        assert(p_equation_index < m_array.get_nb_equations());
        return m_equation_base_variable_index[p_equation_index];
    }

    //-------------------------------------------------------------------------
    template <typename COEF_TYPE>
    unsigned int
    simplex_identity_solver<COEF_TYPE>::get_base_variable_index(unsigned int p_index) const
    {
        assert(p_index < m_array.get_nb_equations());
        return m_base_variables_index[p_index];
    }

    //-------------------------------------------------------------------------
    template <typename COEF_TYPE>
    void
    simplex_identity_solver<COEF_TYPE>::set_base_variable_index(unsigned int p_index
                                                               ,unsigned int p_variable_index
                                                               )
    {
        assert(p_index < m_array.get_nb_equations());
        assert(p_variable_index < m_nb_total_variables);
        m_base_variables_index[p_index] = p_variable_index;
    }

    //-------------------------------------------------------------------------
    template <typename COEF_TYPE>
    bool
    simplex_identity_solver<COEF_TYPE>::get_max_input_column_index(unsigned int & p_variable_index) const
    {
        for(unsigned int l_index = 0; l_index < m_array.get_nb_variables(); ++l_index)
        {
            if(m_array.get_Z_coef(l_index) < 0)
            {
                p_variable_index = l_index;
                return true;
            }
        }
        return false;
    }

    //-------------------------------------------------------------------------
    template <typename COEF_TYPE>
    template <class LISTENER>
    bool
    simplex_identity_solver<COEF_TYPE>::find_max(COEF_TYPE & p_max
                                                ,bool & p_infinite
                                                ,LISTENER *p_listener
                                                )
    {
        p_infinite = false;
        // A' Column index that will be used for Pivot
        unsigned int l_input_column_index = 0;
        unsigned int l_nb_iteration = 0;
        while(get_max_input_column_index(l_input_column_index))
        {
            // Index of variable corresponding to this column
            unsigned int l_input_variable_index = get_array_variable_index(l_input_column_index);
            if(p_listener)
            {
                p_listener->start_iteration(l_nb_iteration);
                p_listener->new_input_var_event(l_input_variable_index);
            }
            // A' line/Equation index that will be used for Pivot
            unsigned int l_output_equation_index = 0;
            if(get_output_equation_index(l_input_column_index, l_output_equation_index))
            {
                // I Column index corresponding to output equation
                unsigned int l_output_base_column_index = get_equation_base_column_index(l_output_equation_index);

                // Output variable index corresponding to Output I column
                unsigned int l_output_variable_index = get_base_variable_index(l_output_base_column_index);
                if(p_listener)
                {
                    p_listener->new_output_var_event(l_output_variable_index);
                }
                set_base_variable_index(l_output_base_column_index, l_input_variable_index);
                set_array_variable_index(l_input_column_index, l_output_variable_index);

                // Pivot
                pivot(l_output_equation_index, l_input_column_index, l_output_base_column_index);

                // Update Z0
                if(p_listener)
                {
                    p_listener->new_Z0(m_array.get_Z0_coef());
                }
            }
            else
            {
                p_infinite = true;
                return false;
            }
            ++l_nb_iteration;
        }
        p_max = m_array.get_Z0_coef();
        return true;
    }

    //-------------------------------------------------------------------------
    template <typename COEF_TYPE>
    bool
    simplex_identity_solver<COEF_TYPE>::get_output_equation_index(unsigned int p_column_index
                                                                 ,unsigned int & p_equation_index
                                                                 ) const
    {
        assert(p_column_index < m_array.get_nb_variables());
        unsigned int l_index = 0;
        while(m_array.get_A_coef(l_index,p_column_index) <= 0 && l_index < m_array.get_nb_equations())
        {
            ++l_index;
        }
        if(l_index == m_array.get_nb_equations())
        {
            return false;
        }
        COEF_TYPE l_min = m_array.get_B_coef(l_index) / m_array.get_A_coef(l_index,p_column_index);
        p_equation_index = l_index;
        ++l_index;
        while(l_index < m_array.get_nb_equations())
        {
            COEF_TYPE l_divider = m_array.get_A_coef(l_index,p_column_index);
            if(l_divider > 0)
            {
                COEF_TYPE l_result = m_array.get_B_coef(l_index) / l_divider;
                if(l_result < l_min)
                {
                    l_min = l_result;
                    p_equation_index = l_index;
                }
            }
            ++l_index;
        }
        return true;
    }

    //-------------------------------------------------------------------------
    template <typename COEF_TYPE>
    void
    simplex_identity_solver<COEF_TYPE>::pivot(const unsigned int p_pivot_row_index
                                             ,const unsigned int p_pivot_column_index
                                             ,const unsigned int p_base_column
                                             )
    {
        assert(p_pivot_row_index < m_array.get_nb_equations());
        assert(p_pivot_column_index < m_array.get_nb_variables());
        COEF_TYPE l_pivot = m_array.get_A_coef(p_pivot_row_index, p_pivot_column_index);
        assert(l_pivot);

        // Pivoting Z row
        //-----------------------
        COEF_TYPE l_q = m_array.get_Z_coef(p_pivot_column_index);
#ifdef DEBUG_SIMPLEX
        if(l_pivot != 1 || l_q)
        {
            std::cout << "Z line <= Z - " << l_q << " * R[" << p_row_index << "] / " << l_pivot << std::endl;
        }
#endif // DEBUG_SIMPLEX
        for(unsigned int l_index = 0;
            l_index < p_pivot_column_index;
            ++l_index
                )
        {
            pivot_Z(l_pivot, l_q, p_pivot_row_index, l_index);
        }
        for(unsigned int l_index = p_pivot_column_index + 1;
            l_index < m_array.get_nb_variables();
            ++l_index
                )
        {
            pivot_Z(l_pivot, l_q, p_pivot_row_index, l_index);
        }
        COEF_TYPE l_u = p_base_column == p_pivot_row_index ? 1 : 0;
        // Z coef is 0 as this column correspond to output base variable
        m_array.set_Z_coef(p_pivot_column_index, - l_q * (l_u / l_pivot));

        m_array.set_Z0_coef(m_array.get_Z0_coef() - (l_q * m_array.get_B_coef(p_pivot_row_index)) / l_pivot);

        // Pivoting other rows except pivot row
        //----------------------
        pivot_A_rows(0,
                     p_pivot_row_index,
                     l_pivot,
                     p_pivot_row_index,
                     p_pivot_column_index,
                     0
                    );
        pivot_A_rows(p_pivot_row_index + 1,
                     m_array.get_nb_equations(),
                     l_pivot,
                     p_pivot_row_index,
                     p_pivot_column_index,
                     0
                    );

        // Particular case of pivot row
        for(unsigned int l_index = 0;
            l_index < p_pivot_column_index;
            ++l_index
                )
        {
            m_array.set_A_coef(p_pivot_row_index,l_index,m_array.get_A_coef(p_pivot_row_index,l_index) / l_pivot);
        }
        for(unsigned int l_index = p_pivot_column_index + 1;
            l_index < m_array.get_nb_variables();
            ++l_index
                )
        {
            m_array.set_A_coef(p_pivot_row_index,l_index,m_array.get_A_coef(p_pivot_row_index,l_index) / l_pivot);
        }
        // Particular case of pivot row/column
        m_array.set_A_coef(p_pivot_row_index, p_pivot_column_index,(p_base_column == p_pivot_row_index ? 1 : 0) / l_pivot);

        // Particulare case of pivot row B coef
        m_array.set_B_coef(p_pivot_row_index, m_array.get_B_coef(p_pivot_row_index) / l_pivot);

    }

    //-------------------------------------------------------------------------
    template <typename COEF_TYPE>
    void
    simplex_identity_solver<COEF_TYPE>::pivot_Z(const COEF_TYPE & p_pivot
                                               ,const COEF_TYPE & p_q
                                               ,unsigned int p_row_index
                                               ,unsigned int p_column_index
                                               )
    {
        COEF_TYPE l_u = m_array.get_A_coef(p_row_index, p_column_index);
        m_array.set_Z_coef(p_column_index, m_array.get_Z_coef(p_column_index) - p_q * (l_u / p_pivot));
    }

    //-------------------------------------------------------------------------
    template <typename COEF_TYPE>
    void
    simplex_identity_solver<COEF_TYPE>::pivot_A_rows(unsigned int p_start_row_index
                                                    ,unsigned int p_end_row_index
                                                    ,const COEF_TYPE & p_pivot
                                                    ,unsigned int p_pivot_row_index
                                                    ,unsigned int p_pivot_column_index
                                                    ,unsigned int p_base_column
                                                    )
    {
        assert(p_end_row_index <= m_array.get_nb_equations());
        for(unsigned int l_row_index = p_start_row_index; l_row_index < p_end_row_index; ++l_row_index)
        {
            COEF_TYPE l_q = m_array.get_A_coef(l_row_index,p_pivot_column_index);
            m_array.set_B_coef(l_row_index, m_array.get_B_coef(l_row_index) - (l_q * m_array.get_B_coef(p_pivot_row_index)) / p_pivot);
            pivot_A_coef(0, p_pivot_column_index, p_pivot, l_q, l_row_index, p_pivot_row_index);
            pivot_A_coef(p_pivot_column_index + 1, m_array.get_nb_variables(), p_pivot, l_q, l_row_index, p_pivot_row_index);
            // Case of pivot column that will be computed using I output variable colum
            COEF_TYPE l_u = p_pivot_row_index == p_base_column ? 1 : 0;
            m_array.set_A_coef(l_row_index, p_pivot_column_index, (l_row_index == p_base_column ?  1 : 0) - (l_q * l_u) / p_pivot);

        }
    }

    //-------------------------------------------------------------------------
    template <typename COEF_TYPE>
    void
    simplex_identity_solver<COEF_TYPE>::pivot_A_coef(unsigned int p_start_column
                                                    ,unsigned int p_end_column
                                                    ,COEF_TYPE p_pivot
                                                    ,COEF_TYPE p_q
                                                    ,unsigned int p_row_index
                                                    ,unsigned int p_pivot_row_index
                                                    )
    {
        assert(p_end_column <= m_array.get_nb_variables());
        for(unsigned int l_column_index = p_start_column;
            l_column_index < p_end_column;
            ++l_column_index
           )
        {
            COEF_TYPE l_u = m_array.get_A_coef(p_pivot_row_index, l_column_index);
            m_array.set_A_coef(p_row_index, l_column_index, m_array.get_A_coef(p_row_index, l_column_index) - (p_q * l_u) / p_pivot);
        }
    }

    //-------------------------------------------------------------------------
    template <typename COEF_TYPE>
    std::vector<COEF_TYPE> simplex_identity_solver<COEF_TYPE>::get_variable_values() const
    {
        std::vector<COEF_TYPE> l_result(m_array.get_nb_variables(), 0);
        for(unsigned int l_index = 0; l_index < m_array.get_nb_equations(); ++l_index)
        {
            if(m_base_variables_index[l_index] < m_array.get_nb_variables())
            {
                l_result[m_base_variables_index[l_index]] = m_array.get_B_coef(l_index);
            }
        }
        return l_result;
    }

}
#endif //SIMPLEX_SIMPLEX_IDENTITY_SOLVER_H
