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

#ifndef SIMPLEX_SOLVER_GLPK_H
#define SIMPLEX_SOLVER_GLPK_H

#include "simplex_solver_base.h"
#include "glpk.h"
#include <string>
#include <map>
#include <cstring>

namespace simplex
{
    class simplex_solver_glpk
    {
      public:
        inline simplex_solver_glpk(unsigned int p_nb_variables
                                  ,unsigned int p_nb_equations
                                  );
        inline ~simplex_solver_glpk();

        /**
         * Define coefficient for objective function
         * @param p_index : the value should be less than number of variables
         * @param value : value of coefficient in the formula Z = SUM(Cj * x)
         */
        inline void set_Z_coef(const unsigned int p_index
                              ,const double & p_value
                              );

        /**
         * Define coefficient for B coefficients in A x = b
         * @param p_index : the value should be less than total number of equations
         * @param value : value of coefficient in b
         */
        inline void set_B_coef(const unsigned int p_index
                              ,const double & p_value
                              );

        /**
         * Define coefficient for A coefficients in A x = b
         * @param p_equation_index : the value should be less than total number of equations
         * @param p_variable_index : the value should be less than number of variables
         * @param value : value of coefficient in A
         */
        inline void set_A_coef(const unsigned int p_equation_index
                              ,const unsigned int p_variable_index
                              ,const double & p_value
                              );

        /**
         * Define equation type
         * @param Equation index
         * @param Equation type
         */
        inline
        void define_equation_type(const unsigned int & p_equation_index
                                 ,const simplex::equation_type & p_equation_type
                                 );

        /**
         * Method implementing simplex algorithm to find max optimum solution
         * @param reference on variable where result will be stored
         * @param reference on a boolean value that will receive true if max
         * is infinite
         * @return value indicating if a max was found
         */
        bool
        find_max(double & p_max
                ,bool & p_infinite
                );

      private:
        static inline
        int convert(const simplex::equation_type & p_equation_type);

        glp_prob * m_problem;
        double * m_B_coefs;
        simplex::equation_type * m_equation_types;
        unsigned int m_nb_equations;
        unsigned int m_nb_variables;
        std::map<std::pair<unsigned int, unsigned int>,double> m_A_coefs;
        bool m_prepared;
    };

    //-------------------------------------------------------------------------
    simplex_solver_glpk::simplex_solver_glpk(unsigned int p_nb_variables
                                            ,unsigned int p_nb_equations
                                            )
    : m_problem(glp_create_prob())
    , m_B_coefs(new double[p_nb_equations])
    , m_equation_types(new simplex::equation_type[p_nb_equations])
    , m_nb_equations(p_nb_equations)
    , m_nb_variables(p_nb_variables)
    , m_prepared(false)
    {
        glp_set_prob_name(m_problem, "Problem");
        glp_add_rows(m_problem, p_nb_equations);
        glp_add_cols(m_problem, p_nb_variables);
        memset(m_equation_types, 0, m_nb_equations * sizeof(simplex::equation_type));
        memset(m_B_coefs, 0, m_nb_equations * sizeof(double));
    }

    //-------------------------------------------------------------------------
    simplex_solver_glpk::~simplex_solver_glpk()
    {
        delete[] m_equation_types;
        delete[] m_B_coefs;
        glp_delete_prob(m_problem);
    }

    //-------------------------------------------------------------------------
    void
    simplex_solver_glpk::set_Z_coef(const unsigned int p_index
                                   ,const double & p_value
                                   )
    {
        assert(p_index < m_nb_variables);
        assert(!m_prepared);
        glp_set_obj_coef(m_problem, 1 + p_index, p_value);
    }

    //-------------------------------------------------------------------------
    void
    simplex_solver_glpk::set_B_coef(const unsigned int p_index
                                   ,const double & p_value
                                   )
    {
        assert(p_index < m_nb_equations);
        assert(!m_prepared);
        m_B_coefs[p_index] = p_value;
    }

    //-------------------------------------------------------------------------
    void
    simplex_solver_glpk::set_A_coef(const unsigned int p_equation_index
                                   ,const unsigned int p_variable_index
                                   ,const double & p_value
                                   )
    {
        assert(p_variable_index < m_nb_variables);
        assert(p_equation_index < m_nb_equations);
        assert(!m_prepared);
        m_A_coefs[std::pair<unsigned int, unsigned int>(p_equation_index, p_variable_index)] = p_value;
    }

    //-------------------------------------------------------------------------
    void
    simplex_solver_glpk::define_equation_type(const unsigned int & p_equation_index
                                             ,const simplex::equation_type & p_equation_type
                                             )
    {
        assert(p_equation_index < m_nb_equations);
        assert(!m_prepared);
        m_equation_types[p_equation_index] = p_equation_type;
    }

    //-------------------------------------------------------------------------
    bool
    simplex_solver_glpk::find_max(double & p_max,
                                  bool & p_infinite
                                 )
    {
        if(!m_prepared)
        {
            // Complete problem description
            // Equation descriptions
            for (unsigned int l_index = 0;
                 l_index < m_nb_equations;
                 ++l_index
                )
            {
                double l_lower_bound = simplex::equation_type::EQUATION == m_equation_types[l_index] ||
                                       simplex::equation_type::INEQUATION_GT == m_equation_types[l_index]
                                       ? m_B_coefs[l_index] : 0;
                double l_upper_bound = simplex::equation_type::INEQUATION_LT == m_equation_types[l_index] ? m_B_coefs[l_index] : 0;
                glp_set_row_bnds(m_problem,
                                 l_index + 1,
                                 convert(m_equation_types[l_index]),
                                 l_lower_bound,
                                 l_upper_bound
                                );
                glp_set_row_name(m_problem,
                                 l_index + 1,
                                 ("E" + std::to_string(l_index)).c_str()
                                );
            }
            // Variable descriptions
            for (unsigned int l_index = 0;
                 l_index < m_nb_variables;
                 ++l_index
                )
            {
                glp_set_col_name(m_problem,
                                 l_index + 1,
                                 ("X" + std::to_string(l_index)).c_str()
                                );
                glp_set_col_bnds(m_problem, l_index + 1, GLP_LO, 0.0, 0.0);
            }
            // Coefficient arrays
            int *l_equation_index_list = new int[m_A_coefs.size() + 1];
            int *l_variable_index_list = new int[m_A_coefs.size() + 1];
            double *l_coef_list = new double[m_A_coefs.size() + 1];
            unsigned int l_index = 1;
            for (auto l_iter: m_A_coefs
                    )
            {
                l_equation_index_list[l_index] = 1 + l_iter.first.first;
                l_variable_index_list[l_index] = 1 + l_iter.first.second;
                l_coef_list[l_index] = l_iter.second;
                ++l_index;
            }
            glp_load_matrix(m_problem, m_A_coefs.size(), l_equation_index_list, l_variable_index_list, l_coef_list);
            //delete[] l_equation_index_list;
            //delete[] l_variable_index_list;
            //delete[] l_coef_list;
        }
        glp_set_obj_dir(m_problem, GLP_MAX);
        glp_simplex(m_problem, NULL);
        std::cout << "STATUS= " << glp_get_status(m_problem) << std::endl;
        for (unsigned int l_index = 0;
             l_index < m_nb_variables;
             ++l_index
                )
        {
            std::cout << glp_get_col_prim(m_problem, 1 + l_index) << std::endl;
        }

        p_max = glp_get_obj_val(m_problem);
        return true;
    }

    //-------------------------------------------------------------------------
    int
    simplex_solver_glpk::convert(const simplex::equation_type & p_equation_type)
    {
        switch(p_equation_type)
        {
            case simplex::equation_type::UNDEFINED:
                return GLP_FR;
            case simplex::equation_type::EQUATION:
                return GLP_FX;
            case simplex::equation_type::INEQUATION_LT:
                return GLP_UP;
            case simplex::equation_type::INEQUATION_GT:
                return GLP_LO;
            default:
                return GLP_FR;

        }
    }

}
#endif //SIMPLEX_SOLVER_GLPK_H
