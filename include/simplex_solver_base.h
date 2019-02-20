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

#ifndef SIMPLEX_SOLVER_BASE_H
#define SIMPLEX_SOLVER_BASE_H

#include "simplex_listener.h"
#include "quicky_exception.h"
#include "equation_system.h"
#include <sstream>
#include <limits>
#include <vector>
#include <memory.h>

namespace simplex
{

    typedef enum class equation_type
    {
        UNDEFINED = 0,
        EQUATION,
        INEQUATION_LT,
        INEQUATION_GT
    } t_equation_type;

    template <typename COEF_TYPE,typename ARRAY_TYPE>
    class simplex_solver_base: public simplex_listener_target_if<COEF_TYPE>
    {
      public:
        typedef COEF_TYPE t_coef_type;
        typedef ARRAY_TYPE t_array_type;

        /**
         * This constructor prepare a simplex in Canonical form
         * Max z = cx
         * a1 * x1 + a2 * x2 + ... + an * xn <= bn
         * xi >= 0 for all xi
         *
         * Or
         * Min z = cx
         * a1 * x1 + a2 * x2 + ... + an * xn >= bn
         * xi >= 0 for all xi
         */
        inline simplex_solver_base(unsigned int p_nb_variables,
                                   unsigned int p_nb_inequations_lt,
                                   unsigned int p_nb_equations,
                                   unsigned int p_nb_inequations_gt
                                  );

        /**
         * Define coefficient for objective function
         * @param p_index : the value should be less than number of variables
         * @param value : value of coefficient in the formula Z = SUM(Cj * x)
         */
        inline void set_Z_coef(const unsigned int p_index,
                               const COEF_TYPE & p_value
                              );

        /**
         * Define coefficient for B coefficients in A x = b
         * @param p_index : the value should be less than total number of equations
         * @param value : value of coefficient in b
         */
        inline void set_B_coef(const unsigned int p_index,
                               const COEF_TYPE & p_value
                              );

        /**
         * Define coefficient for A coefficients in A x = b
         * @param p_equation_index : the value should be less than total number of equations
         * @param p_variable_index : the value should be less than number of variables
         * @param value : value of coefficient in A
         */
        inline void set_A_coef(const unsigned int p_equation_index,
                               const unsigned int p_variable_index,
                               const COEF_TYPE & p_value
                              );

        /**
         * Return coefficient for A coefficients in A x = b
         * @param p_equation_index : the value should be less than total number of equations
         * @param p_variable_index : the value should be less than number of variables
         * @return value : value of coefficient in A
         */
        inline
        const COEF_TYPE &
        get_A_coef(const unsigned int p_equation_index,
                   const unsigned int p_variable_index
                  ) const;

        /**
         * Display simplex array representation
         * @param p_stream stream where the display should be done
         * @return the modified stream
         */
        inline std::ostream & display_array(std::ostream & p_stream)const override;

        /**
         * Return value of variables
         * @return value of variables
         */
        std::vector<COEF_TYPE> get_variable_values() const override;

        /**
         * Define equation type
         * @param Equation index
         * @param Equation type
         */
        inline
        void define_equation_type(const unsigned int & p_equation_index,
                                  const simplex::equation_type & p_equation_type
                                 );

        /**
         * Declare that a variable is a base variable
         * @param variable index in simplex array
         */
        inline
        void define_base_variable(const unsigned int & p_variable_index);

        /**
         * Return variable index stored as base variable defined by parameter
         * @param index of base variable
         * @return index of variable used as base variable
         */
        inline
        const unsigned int &
        get_base_variable(const unsigned int & p_index)const;

        /**
         * Return total number of equations ( <= + >= + = )
         */
        inline
        const unsigned int
        get_total_nb_equation(void)const;

        /**
         Method implementing simplex algorithm to find max optimum solution
         The problem must be in solved form
         @param reference on variable where result will be stored
         @param reference on a boolean value that will receive true if max
         is infinite
         @return value indicating if a max was found
         */
        template <class LISTENER=simplex_listener<COEF_TYPE>>
        bool
        find_max(COEF_TYPE & p_max,
                 bool & p_infinite,
                 LISTENER * p_listener = NULL
                );

        /**
         *  Method checking if variable values passed as parameters respect
         *  simplex problem constraints
         * @param p_values variable values
         * @return true if constraints are respected
         */
        inline bool check_variables(const std::vector<COEF_TYPE> & p_values);

        inline virtual ~simplex_solver_base(void);

      protected:

        bool is_base_variable(const unsigned int & p_index)const;

        unsigned int
        get_nb_all_variables() const;

        unsigned int
        get_nb_total_equations() const;

        const ARRAY_TYPE &
        get_array() const;

        ARRAY_TYPE &
        get_array();

        unsigned int
        get_nb_variables() const;

        unsigned int
        get_nb_base_variables_defined() const;

        unsigned int
        get_base_variables_position(unsigned int p_index) const;

        void
        set_base_variables_position(unsigned int p_index,
                                    unsigned int p_value
                                   );

        /**
         * Set index of base variable in equation defined by parameter
         * @param p_index of base variable
         * @param p_value index of variable used as base variable
         */
        inline
        void
        set_base_variable(const unsigned int & p_index, const unsigned int & p_value)const;

        /**
            * Define coefficient for A coefficients in A x = b
            * @param p_equation_index : the value should be less than total number of equations
            * @param p_variable_index : the value should be less than total number of variables
            * including adjustment variables
            * @param value : value of coefficient in A
            */
        inline
        void
        set_internal_coef(const unsigned int p_equation_index,
                          const unsigned int p_variable_index,
                          const COEF_TYPE & p_value
                         );

        /**
              * Return coefficient for A coefficients in A x = b
              * @param p_equation_index : the value should be less than total number of equations
              * @param p_variable_index : the value should be less than total number of variables
              * including adjustment variables
              * @return value : value of coefficient in A
              */
        inline
        const COEF_TYPE &
        get_internal_coef(const unsigned int p_equation_index,
                          const unsigned int p_variable_index
                         ) const;

        /**
           * Method to determine the next input variable for pivot operation
           * when searching max optimum
           * @param reference on variable where to store the input variable index if any
           * @return boolean indicating if an input variable was found: TRUE = found
           */
        inline
        bool
        get_max_input_variable_index(unsigned int & p_variable_index) const;

      private:

        /**
           Method performing pivot to change the base
           The A coefficient A[row,column] should be !0
           @param p_row_index Row index
           @param p_column_index Column index
        */
        virtual void pivot(const unsigned int p_row_index, 
                           const unsigned int p_column_index
                          ) = 0;

        /**
           Method to determine the equation index corresponding to next output
           variable for pivot operation
           @param index of input variable
           @param reference on variable where to store the output equation index if any
           @return boolean indicating if an input variable was found
         */
        virtual
        bool
        get_output_equation_index(unsigned int p_input_variable_index,
                                  unsigned int & p_equation_index
                                 )const = 0;

        /**
         * Method to set adjustement variable in the right place
         * @param index of inequation
         * @param value to set as coefficient
         */
        inline
        void
        set_adjustement_variable(const unsigned int & p_equation_index,
                                 const COEF_TYPE & p_value
                                );

        /**
         * Method to determine the next input variable for pivot operation
         * when searching min optimum
         * @param reference on variable where to store the input variable index if any
         * @return boolean indicating if an input variable was found: TRUE = found
         */
        inline
        bool
        get_min_input_variable_index(unsigned int & p_variable_index) const;

        simplex_solver_base(void) = delete;

        /**
         * Variable numbers without adjustment variables
         */
        unsigned int m_nb_variables;

        /**
         * Number of inequations of type x1 + x2 + ... + xn <= value
         */
        unsigned int m_nb_inequations_lt;

        /**
         * Number of equations of type x1 + x2 + ... + xn = value
         */
        unsigned int m_nb_equations;

        /**
         * Number of inequations of type x1 + x2 + ... + xn >= value
         */
        unsigned int m_nb_inequations_gt;

        /**
         * Variable counting the number of base variables
         */
        unsigned int m_nb_adjustment_variable;

        /**
         * Variable counting the number of adjustement variable specified
         */
        unsigned int m_nb_defined_adjustment_variables;

        /**
         * Total number of variables including adjustment variables
         */
        unsigned int m_nb_all_variables;

        /**
         * Total numer of equations after adjuments variables have been added
         */
        unsigned int m_nb_total_equations;

        /**
         * Array storing equations
         */
        ARRAY_TYPE m_array;

        /**
         * Equation types
         */
        t_equation_type * m_equation_types;

        /**
         * Base variables, array to store base variables index
         */
        unsigned int * m_base_variables;

        /**
         * Array storing if variable is a base variable or not
         */
        unsigned int * m_base_variables_position;

        /**
         * Nb base variables defined
         */
        unsigned int m_nb_base_variables_defined;
    };

    //----------------------------------------------------------------------------
    std::ostream &
    operator<<(std::ostream & p_stream,
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
    template <typename COEF_TYPE,typename ARRAY_TYPE>
    simplex::simplex_solver_base<COEF_TYPE, ARRAY_TYPE>::simplex_solver_base(unsigned int p_nb_variables,
                                                                             unsigned int p_nb_inequations_lt,
                                                                             unsigned int p_nb_equations,
                                                                             unsigned int p_nb_inequations_gt
                                                                            ):
            m_nb_variables(p_nb_variables),
            m_nb_inequations_lt(p_nb_inequations_lt),
            m_nb_equations(p_nb_equations),
            m_nb_inequations_gt(p_nb_inequations_gt),
            m_nb_adjustment_variable(p_nb_inequations_lt + p_nb_inequations_gt),
            m_nb_defined_adjustment_variables(0),
            m_nb_all_variables(p_nb_variables + m_nb_adjustment_variable),
            m_nb_total_equations(p_nb_inequations_lt + p_nb_equations + p_nb_inequations_gt),
            m_array(m_nb_total_equations,m_nb_all_variables),
            m_equation_types(new simplex::equation_type[m_nb_total_equations]),
            m_base_variables(new unsigned int[m_nb_total_equations]),
            m_base_variables_position(new unsigned int[m_nb_all_variables]),
            m_nb_base_variables_defined(0)
    {
        static_assert(std::is_signed<COEF_TYPE>::value, "Simplex template parameter should be signed");
        for(unsigned int l_index = 0;
            l_index < m_nb_total_equations;
            ++l_index
                )
        {
            m_base_variables[l_index] = ::std::numeric_limits<unsigned int>::max();
        }
        memset(m_base_variables_position, 0xFF, m_nb_all_variables * sizeof(unsigned int));
        memset(m_equation_types, 0, m_nb_total_equations * sizeof(simplex::equation_type));
    }

    //----------------------------------------------------------------------------
    template <typename COEF_TYPE,typename ARRAY_TYPE>
    void simplex::simplex_solver_base<COEF_TYPE, ARRAY_TYPE>::set_Z_coef(const unsigned int p_index,
                                                                         const COEF_TYPE & p_value
                                                                        )
    {
        assert(p_index < m_nb_variables);
        m_array.set_Z_coef(p_index,-p_value);
    }

    //----------------------------------------------------------------------------
    template <typename COEF_TYPE,typename ARRAY_TYPE>
    void simplex::simplex_solver_base<COEF_TYPE, ARRAY_TYPE>::set_B_coef(const unsigned int p_index,
                                                                         const COEF_TYPE & p_value
                                                                        )
    {
        m_array.set_B_coef(p_index, p_value);
    }

    //----------------------------------------------------------------------------
    template <typename COEF_TYPE,typename ARRAY_TYPE>
    void simplex::simplex_solver_base<COEF_TYPE, ARRAY_TYPE>::set_A_coef(const unsigned int p_equation_index,
                                                                         const unsigned int p_variable_index,
                                                                         const COEF_TYPE & p_value
                                                                        )
    {
        // Keep assert because we differentiate access to internal coefs
        assert(p_variable_index < m_nb_variables);
        m_array.set_A_coef(p_equation_index,p_variable_index,p_value);
    }

    //----------------------------------------------------------------------------
    template <typename COEF_TYPE,typename ARRAY_TYPE>
    const COEF_TYPE &
    simplex::simplex_solver_base<COEF_TYPE, ARRAY_TYPE>::get_A_coef(const unsigned int p_equation_index,
                                                                    const unsigned int p_variable_index
                                                                   )const
    {
        // Keep assert because we differentiate access to internal coefs
        assert(p_variable_index < m_nb_variables);
        return m_array.get_A_coef(p_equation_index,p_variable_index);
    }

    //----------------------------------------------------------------------------
    template <typename COEF_TYPE,typename ARRAY_TYPE>
    void simplex::simplex_solver_base<COEF_TYPE, ARRAY_TYPE>::set_internal_coef(const unsigned int p_equation_index,
                                                                                const unsigned int p_variable_index,
                                                                                const COEF_TYPE & p_value
                                                                               )
    {
        // Keep assert because we differentiate access to internal coefs
        assert(p_variable_index < m_nb_all_variables);
        //    assert(p_variable_index >= m_nb_variables);
        m_array.set_A_coef(p_equation_index, p_variable_index, p_value);
    }

    //----------------------------------------------------------------------------
    template <typename COEF_TYPE,typename ARRAY_TYPE>
    const COEF_TYPE &
    simplex::simplex_solver_base<COEF_TYPE, ARRAY_TYPE>::get_internal_coef(const unsigned int p_equation_index,
                                                                           const unsigned int p_variable_index
                                                                          )const
    {
        // Keep assert because we differentiate access to internal coefs
        assert(p_variable_index < m_nb_all_variables);
        //      assert(p_variable_index >= m_nb_variables);
        return m_array.get_A_coef(p_equation_index, p_variable_index);
    }

    //----------------------------------------------------------------------------
    template <typename COEF_TYPE,typename ARRAY_TYPE>
    std::ostream & simplex::simplex_solver_base<COEF_TYPE, ARRAY_TYPE>::display_array(std::ostream & p_stream)const
    {
        p_stream << "Z\t";
        for(unsigned int l_index = 0;
            l_index < m_nb_all_variables;
            ++l_index
                )
        {
            if(l_index < m_nb_variables)
            {
                p_stream << "X" << 1 + l_index;
            }
            else
            {
                p_stream << "E" << 1 + l_index - m_nb_variables;
            }
            p_stream << "\t" ;
        }
        p_stream << std::endl;
        p_stream << "Z\t";
        for(unsigned int l_index = 0;
            l_index < m_nb_all_variables;
            ++l_index
                )
        {
            p_stream << m_array.get_Z_coef(l_index) << "\t";
        }
        p_stream << "|\t" << m_array.get_Z0_coef() << std::endl;
        for(unsigned int l_row_index = 0;
            l_row_index < m_nb_total_equations;
            ++l_row_index
                )
        {
            unsigned int l_var_index = m_base_variables[l_row_index];
            if(::std::numeric_limits<unsigned int>::max() != l_var_index)
            {
                if(l_var_index < m_nb_variables)
                {
                    p_stream << "X" << 1 + l_var_index;
                }
                else
                {
                    p_stream << "E" << 1 + l_var_index - m_nb_variables;
                }
            }
            else
            {
                p_stream << " ";
            }
            p_stream << "\t";
            for(unsigned int l_index = 0;
                l_index < m_nb_all_variables;
                ++l_index
                    )
            {
                p_stream << get_internal_coef(l_row_index, l_index) << "\t";
            }
            p_stream << "|\t" << m_array.get_B_coef(l_row_index) << std::endl;
        }
        return p_stream;
    }

    //----------------------------------------------------------------------------
    template <typename COEF_TYPE,typename ARRAY_TYPE>
    void simplex::simplex_solver_base<COEF_TYPE, ARRAY_TYPE>::set_adjustement_variable(const unsigned int & p_equation_index,
                                                                                       const COEF_TYPE & p_value
                                                                                      )
    {
        unsigned int l_column_index = m_nb_variables + m_nb_defined_adjustment_variables;
        set_internal_coef(p_equation_index, l_column_index, p_value);
        if(!m_nb_equations)
        {
            define_base_variable(l_column_index);
        }
        ++m_nb_defined_adjustment_variables;
    }

    //----------------------------------------------------------------------------
    template <typename COEF_TYPE,typename ARRAY_TYPE>
    bool simplex::simplex_solver_base<COEF_TYPE, ARRAY_TYPE>::get_max_input_variable_index(unsigned int & p_variable_index) const
    {
        for(unsigned int l_index = 0;
            l_index < m_nb_all_variables;
            ++l_index
                )
        {
            if(m_array.get_Z_coef(l_index) < 0)
            {
                p_variable_index = l_index;
                return true;
            }
        }
        return false;
    }

    //----------------------------------------------------------------------------
    template <typename COEF_TYPE,typename ARRAY_TYPE>
    bool simplex::simplex_solver_base<COEF_TYPE, ARRAY_TYPE>::get_min_input_variable_index(unsigned int & p_variable_index) const
    {
        for(unsigned int l_index = 0;
            l_index < m_nb_all_variables;
            ++l_index
                )
        {
            if(m_array.get_Z_coef(l_index) > 0)
            {
                p_variable_index = l_index;
                return true;
            }
        }
        return false;
    }

    //----------------------------------------------------------------------------
    template <typename COEF_TYPE,typename ARRAY_TYPE>
    void simplex_solver_base<COEF_TYPE,ARRAY_TYPE>::define_equation_type(const unsigned int & p_equation_index,
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
                set_adjustement_variable(p_equation_index,(COEF_TYPE)1);
                break;
            case t_equation_type::INEQUATION_GT:
                set_adjustement_variable(p_equation_index,(COEF_TYPE)-1);
                break;
            default:
                throw quicky_exception::quicky_logic_exception("Unknown equation_type value : "+ std::to_string((unsigned int)p_equation_type),__LINE__,__FILE__);
        }
    }

    //----------------------------------------------------------------------------
    template <typename COEF_TYPE,typename ARRAY_TYPE>
    void simplex_solver_base<COEF_TYPE,ARRAY_TYPE>::define_base_variable(const unsigned int & p_variable_index)
    {
        assert(p_variable_index < m_nb_all_variables);
        assert(m_nb_base_variables_defined < m_nb_total_equations);
        m_base_variables[m_nb_base_variables_defined] = p_variable_index;
        m_base_variables_position[p_variable_index] = m_nb_base_variables_defined;
        ++m_nb_base_variables_defined;
    }

    //----------------------------------------------------------------------------
    template <typename COEF_TYPE,typename ARRAY_TYPE>
    const unsigned int & simplex_solver_base<COEF_TYPE,ARRAY_TYPE>::get_base_variable(const unsigned int & p_index)const
    {
        assert(p_index < m_nb_total_equations);
        return m_base_variables[p_index];
    }

    //----------------------------------------------------------------------------
    template <typename COEF_TYPE,typename ARRAY_TYPE>
    void
    simplex_solver_base<COEF_TYPE,ARRAY_TYPE>::set_base_variable(const unsigned int & p_index,
                                                                 const unsigned int & p_value
                                                                ) const
    {
        assert(p_index < m_nb_total_equations);
        m_base_variables[p_index] = p_value;
    }

    //----------------------------------------------------------------------------
    template <typename COEF_TYPE,typename ARRAY_TYPE>
    const unsigned int simplex_solver_base<COEF_TYPE,ARRAY_TYPE>::get_total_nb_equation(void)const
    {
        return m_nb_total_equations;
    }

    //----------------------------------------------------------------------------
    template <typename COEF_TYPE,typename ARRAY_TYPE>
    simplex_solver_base<COEF_TYPE,ARRAY_TYPE>::~simplex_solver_base(void)
    {
        delete[] m_base_variables_position;
        delete[] m_base_variables;
        delete[] m_equation_types;
        m_equation_types = nullptr;
    }

    //----------------------------------------------------------------------------
    template <typename COEF_TYPE,typename ARRAY_TYPE>
    unsigned int
    simplex_solver_base<COEF_TYPE,ARRAY_TYPE>::get_nb_total_equations() const
    {
        return m_nb_total_equations;
    }

    //----------------------------------------------------------------------------
    template <typename COEF_TYPE,typename ARRAY_TYPE>
    unsigned int
    simplex_solver_base<COEF_TYPE,ARRAY_TYPE>::get_nb_all_variables() const
    {
        return m_nb_all_variables;
    }

    //----------------------------------------------------------------------------
    template <typename COEF_TYPE,typename ARRAY_TYPE>
    const ARRAY_TYPE &
    simplex_solver_base<COEF_TYPE,ARRAY_TYPE>::get_array() const
    {
        return m_array;
    }

    //----------------------------------------------------------------------------
    template <typename COEF_TYPE,typename ARRAY_TYPE>
    ARRAY_TYPE &
    simplex_solver_base<COEF_TYPE,ARRAY_TYPE>::get_array()
    {
        return m_array;
    }

    //----------------------------------------------------------------------------
    template <typename COEF_TYPE,typename ARRAY_TYPE>
    unsigned int
    simplex_solver_base<COEF_TYPE,ARRAY_TYPE>::get_nb_variables() const
    {
        return m_nb_variables;
    }

    //----------------------------------------------------------------------------
    template <typename COEF_TYPE,typename ARRAY_TYPE>
    unsigned int
    simplex_solver_base<COEF_TYPE,ARRAY_TYPE>::get_nb_base_variables_defined() const
    {
        return m_nb_base_variables_defined;
    }

    //----------------------------------------------------------------------------
    template <typename COEF_TYPE,typename ARRAY_TYPE>
    unsigned int
    simplex_solver_base<COEF_TYPE,ARRAY_TYPE>::get_base_variables_position(unsigned int p_index) const
    {
        assert(p_index < m_nb_all_variables);
        return m_base_variables_position[p_index];
    }

    //----------------------------------------------------------------------------
    template <typename COEF_TYPE,typename ARRAY_TYPE>
    void
    simplex_solver_base<COEF_TYPE,ARRAY_TYPE>::set_base_variables_position(unsigned int p_index,
                                                                           unsigned int p_value
                                                                          )
    {
        assert(p_index < m_nb_all_variables);
        m_base_variables_position[p_index] = p_value;
    }

    //----------------------------------------------------------------------------
    template <typename COEF_TYPE,typename ARRAY_TYPE>
    template <class LISTENER>
    bool
    simplex_solver_base<COEF_TYPE,ARRAY_TYPE>::find_max(COEF_TYPE & p_max,
                                                        bool & p_infinite,
                                                        LISTENER * p_listener
                                                       )
    {
        if(m_nb_base_variables_defined != m_nb_total_equations)
        {
            throw quicky_exception::quicky_runtime_exception("Not enough base variables defined : " + std::to_string(m_nb_base_variables_defined) + " < " + std::to_string(m_nb_total_equations) , __LINE__, __FILE__);
        }
        for(unsigned int l_index = 0;
            l_index < m_nb_total_equations;
            ++l_index
                )
        {
            unsigned int l_var_index = m_base_variables[l_index];
            assert(l_var_index != ::std::numeric_limits<unsigned int>::max());
            if(m_array.get_Z_coef(l_var_index))
            {
                throw quicky_exception::quicky_runtime_exception("Z coef of base variable in column " + std::to_string(l_var_index) + " should be 0 or this is not a base variable", __LINE__, __FILE__);
            }
        }
        p_infinite = false;
        bool l_input_found = true;
        unsigned int l_input_variable_index = 0;
        unsigned int l_nb_iteration = 0;
        while(true == (l_input_found = get_max_input_variable_index(l_input_variable_index)))
        {
            if(p_listener)
            {
                p_listener->start_iteration(l_nb_iteration);
                p_listener->new_input_var_event(l_input_variable_index);
            }
            assert(::std::numeric_limits<unsigned int>::max() == m_base_variables_position[l_input_variable_index]);
            assert(m_array.get_Z_coef(l_input_variable_index));
            unsigned int l_output_equation_index = 0;
            if(get_output_equation_index(l_input_variable_index, l_output_equation_index))
            {
                unsigned int l_output_variable_index = m_base_variables[l_output_equation_index];
                if(p_listener)
                {
                    p_listener->new_output_var_event(l_output_variable_index);
                }
                assert(l_output_equation_index == m_base_variables_position[l_output_variable_index]);
                assert(!m_array.get_Z_coef(l_output_variable_index));
                pivot(l_output_equation_index,l_input_variable_index);
                assert(m_array.get_Z_coef(l_output_variable_index));
                assert(!m_array.get_Z_coef(l_input_variable_index));
                m_base_variables_position[l_output_variable_index] = std::numeric_limits<unsigned int>::max();
                m_base_variables_position[l_input_variable_index] = l_output_equation_index;
                m_base_variables[l_output_equation_index] = l_input_variable_index;

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
    template <typename COEF_TYPE,
              typename ARRAY_TYPE
             >
    bool
    simplex_solver_base<COEF_TYPE, ARRAY_TYPE>::is_base_variable(const unsigned int & p_index) const
    {
        return std::numeric_limits<unsigned int>::max() != m_base_variables_position[p_index];
    }

    //-------------------------------------------------------------------------
    template <typename COEF_TYPE,
              typename ARRAY_TYPE
             >
    bool
    simplex_solver_base<COEF_TYPE, ARRAY_TYPE>::check_variables(const std::vector<COEF_TYPE> & p_values)
    {
        assert(p_values.size() <= m_nb_variables);
        for(unsigned int l_equation_index = 0;
            l_equation_index < m_nb_equations;
            ++l_equation_index
           )
        {
            COEF_TYPE l_sum(0);
            for(unsigned int l_variable_index = 0;
                l_variable_index < p_values.size();
                ++l_variable_index
               )
            {
                l_sum += m_array.get_A_coef(l_equation_index, l_variable_index) * p_values[l_variable_index];
            }
            COEF_TYPE l_b_coef = m_array.get_B_coef(l_equation_index);
            switch(m_equation_types[l_equation_index])
            {
                case t_equation_type::UNDEFINED:
                    throw quicky_exception::quicky_logic_exception("Undefined equation when checking simplex variables", __LINE__, __FILE__);
                    break;
                case t_equation_type::EQUATION:
                    if(l_sum != l_b_coef)
                    {
                        return false;
                    }
                    break;
                case t_equation_type::INEQUATION_LT:
                    if(l_sum > l_b_coef)
                    {
                        return false;
                    }
                    break;
                case t_equation_type::INEQUATION_GT:
                    if(l_sum < l_b_coef)
                    {
                        return false;
                    }
            }
        }
        return true;
    }

    //-------------------------------------------------------------------------
    template <typename COEF_TYPE, typename ARRAY_TYPE>
    std::vector<COEF_TYPE>
    simplex_solver_base<COEF_TYPE, ARRAY_TYPE>::get_variable_values() const
    {
        std::vector<COEF_TYPE> l_result(m_nb_variables, (COEF_TYPE)0);
        my_square_matrix<COEF_TYPE> l_matrix(m_nb_equations);
        my_matrix<COEF_TYPE> l_coefs(m_nb_equations, 1);
        // Create Equation system
        for(unsigned int l_equation_index = 0; l_equation_index < m_nb_equations; ++l_equation_index)
        {
            for(unsigned int l_variable_index = 0; l_variable_index < m_nb_equations; ++l_variable_index)
            {
                l_matrix.set_data(l_equation_index, l_variable_index, m_array.get_A_coef(l_equation_index, get_base_variable(l_variable_index)));
            }
            l_coefs.set_data(l_equation_index, 0, m_array.get_B_coef(l_equation_index));
        }
        my_equation_system<COEF_TYPE> l_equation_system(l_matrix, l_coefs);

        my_matrix<COEF_TYPE> l_result_matrix = l_equation_system.solve();
        assert(l_result_matrix.get_height() == m_nb_equations);

        // Fill result
        for(unsigned int l_index = 0; l_index < m_nb_equations; ++l_index)
        {
            l_result[get_base_variable(l_index)] = l_result_matrix.get_data(l_index, 0);
        }
        return l_result;
    }
}
#endif //SIMPLEX_SOLVER_BASE_H
// EOF
