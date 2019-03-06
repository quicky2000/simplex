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

#include "simplex_listener.h"
#include "simplex_solver_base.h"
#include "glpk.h"
#include <string>
#include <map>
#include <cstring>
#include <cassert>

namespace simplex
{
    class simplex_solver_glpk: public simplex_listener_target_if<double>
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
         * @param p_max reference on variable where result will be stored
         * @param p_infinite reference on a boolean value that will receive true if max
         * is infinite
         * @param p_listener optional listener to treat iterations information
         * @return value indicating if a max was found
         */
        bool
        find_max(double & p_max
                ,bool & p_infinite
                ,simplex_listener<double> *p_listener = NULL
                );

        /**
         * Return values of variable at current iteration
         * @return
         */
        std::vector<double> get_variable_values() const override;

        /**
         * Supposed to display array content
         * @param p_stream
         * @return
         */
        std::ostream & display_array(std::ostream & p_stream) const override;
      private:
        /**
         * Method to intercept terminal output from GLPK
         * @param p_info additional info provided when registering hook ( this)
         * @param p_msg message send to terminal by GLPK
         * @return 0 to let GPLK display message, something else to not display
         */
        static inline
        int terminal_hook(void *p_info
                         ,const char *p_msg
                         );

        /**
         * Method treating messages sent by GLPK to terminal output
         * @param p_msg message to be analyzed
         */
        void treat_message(const std::string & p_msg);

        static inline
        int convert(const simplex::equation_type & p_equation_type);

        static inline
        std::string status_to_string(int p_status);

        glp_prob * m_problem;
        double * m_B_coefs;
        simplex::equation_type * m_equation_types;
        unsigned int m_nb_equations;
        unsigned int m_nb_variables;
        std::map<std::pair<unsigned int, unsigned int>,double> m_A_coefs;
        bool m_prepared;

        /**
         * Optional listener when calling find_max
         * Member value to be accessible from treat_message method
         */
        simplex_listener<double> * m_listener;

        /**
         * Store iteration number
         */
        uint64_t m_iteration;
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
    , m_listener(nullptr)
    , m_iteration(0)
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
    simplex_solver_glpk::find_max(double & p_max
                                 ,bool & p_infinite
                                 ,simplex_listener<double> * p_listener
                                 )
    {
        m_listener = p_listener;
        m_iteration = 0;
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
            delete[] l_equation_index_list;
            delete[] l_variable_index_list;
            delete[] l_coef_list;
        }
        glp_set_obj_dir(m_problem, GLP_MAX);
        glp_smcp l_solver_parameters;
        glp_init_smcp(&l_solver_parameters);

        if(m_listener)
        {
            // Intercept terminal outputs to get iteration number
            glp_term_hook(simplex_solver_glpk::terminal_hook, this);
            l_solver_parameters.out_frq = 1;
            l_solver_parameters.msg_lev = GLP_MSG_ALL;

            // Limit number of iteration because we saw in treat_message that
            // variable values are not accessible during call of GLPK solver
            // By this way we will got out of solver when iteration limit is
            // reached an be able to access to variable values
            l_solver_parameters.it_lim = 1;
        }

        do
        {
            glp_simplex(m_problem, &l_solver_parameters);
            if(m_listener)
            {
                m_listener->new_Z0(glp_get_obj_val(m_problem));
            }
#ifdef DEBUG_SIMPLEX_SOLVER_GLPK
            std::cout << "STATUS= " << status_to_string(glp_get_status(m_problem)) << std::endl;
#endif // DEBUG_SIMPLEX_SOLVER_GLPK
        }
        while(GLP_FEAS == glp_get_status(m_problem));
        p_max = glp_get_obj_val(m_problem);
        m_listener = NULL;
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

    //-------------------------------------------------------------------------
    std::string
    simplex_solver_glpk::status_to_string(int p_status)
    {
        switch(p_status)
        {
            case GLP_OPT:
                return "Optimal solution";
            case GLP_FEAS:
                return "Feasible solution";
            case GLP_INFEAS:
                return "Ineasible solution";
            case GLP_NOFEAS:
                return "No feasible solution";
            case GLP_UNBND:
                return "Unbounded solution";
            case GLP_UNDEF:
                return "Undefined solution";
            default:
                throw quicky_exception::quicky_logic_exception("Uknown GLPK solver status: " + std::to_string(p_status), __LINE__, __FILE__);
        }
    }

    //-------------------------------------------------------------------------
    int
    simplex_solver_glpk::terminal_hook(void *p_info
                                      ,const char *p_msg
                                      )
    {
        static_cast<simplex_solver_glpk*>(p_info)->treat_message(p_msg);
        return 1;
    }

    //-------------------------------------------------------------------------
    void
    simplex_solver_glpk::treat_message(const std::string & p_msg)
    {
#ifdef DEBUG_SIMPLEX_SOLVER_GLPK
        std::cout << "Message : " << p_msg;
#endif // DEBUG_SIMPLEX_SOLVER_GLPK
        size_t l_pos = p_msg.find(':');

        // Check if we are in step message
        if(std::string::npos == l_pos)
        {
            return;
        }
        // We start at second character as first one indicate the phase of simplex method
        unsigned long l_step = std::stoul(p_msg.substr(1, l_pos - 1));

        assert(m_listener);
        if(m_iteration == l_step)
        {
            m_listener->start_iteration(l_step);
        }
        else
        {
            m_iteration = l_step;
        }

#ifdef DEBUG_SIMPLEX_SOLVER_GLPK
        l_pos = p_msg.find('=', l_pos);
        assert(std::string::npos != l_pos);
        // Skip spaces
        l_pos = p_msg.find_first_not_of(" \t", l_pos + 1);
        assert(std::string::npos != l_pos);
        std::string l_substr(p_msg.substr(l_pos));

        l_pos = l_substr.find("inf");
        assert(std::string::npos != l_pos);

        std::string l_objective_str = l_substr.substr(0, l_pos);
        double l_objective = std::stod(l_objective_str);
        std::cout << "Step: " << l_step << "\tObjective : " << l_objective << std::endl;
#endif // DEBUG_SIMPLEX_SOLVER_GLPK
    }

    //-------------------------------------------------------------------------
    std::vector<double>
    simplex_solver_glpk::get_variable_values() const
    {
        assert(m_problem);
        assert(GLP_FEAS == glp_get_status(m_problem) || GLP_OPT == glp_get_status(m_problem));
        std::vector<double> l_result;
        for (unsigned int l_index = 0;
             l_index < m_nb_variables;
             ++l_index
                )
        {
            l_result.push_back(glp_get_col_prim(m_problem, 1 + l_index));
        }
        return l_result;
    }

    //-------------------------------------------------------------------------
    std::ostream &
    simplex_solver_glpk::display_array(std::ostream & p_stream) const
    {
        // Not implemented because don't know how to access coefficients at
        // this time
        return p_stream;
    }

}
#endif //SIMPLEX_SOLVER_GLPK_H
