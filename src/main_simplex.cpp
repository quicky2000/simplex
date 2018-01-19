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
#ifdef SIMPLEX_SELF_TEST
#include <iostream>
#include <fstream>
#include "safe_types.h"
#include "fract.h"
#include "simplex_listener.h"
#include "simplex_map.h"
#include "simplex_solver.h"
#include "simplex_solver_integer.h"

template <typename SIMPLEX_TYPE>
void test_case1(void);

template <typename SIMPLEX_TYPE>
void test_case2(void);

template <typename SIMPLEX_TYPE>
void test_case3(void);

//------------------------------------------------------------------------------
int main(int argc,char ** argv)
{
    try
    {
        std::cout << "============ TEST CASE 1 ==============" << std::endl;
        test_case1<simplex::simplex_solver<double>>();
        std::cout << "============ TEST CASE 1 ==============" << std::endl;
        test_case1<simplex::simplex_solver_integer<int32_t>>();
        std::cout << "============ TEST CASE 1 ==============" << std::endl;
        test_case1<simplex::simplex_solver<double,simplex::simplex_map<double>>>();
        std::cout << "============ TEST CASE 1 bis==============" << std::endl;
        test_case1<simplex::simplex_solver<quicky_utils::fract<uint32_t>,simplex::simplex_map<quicky_utils::fract<uint32_t>>>>();
        std::cout << "============ TEST CASE 1 ter ==========" << std::endl;
        test_case1<simplex::simplex_solver<quicky_utils::fract<quicky_utils::safe_uint32_t>,simplex::simplex_map<quicky_utils::fract<quicky_utils::safe_uint32_t>>>>();
        std::cout << "============ TEST CASE 2 ==============" << std::endl;
        test_case2<simplex::simplex_solver<double>>();
        std::cout << "============ TEST CASE 3 ==============" << std::endl;
        test_case3<simplex::simplex_solver<double>>();
    }
    catch(quicky_exception::quicky_runtime_exception & e)
    {
        std::cout << "ERROR : " << e.what() << std::endl ;
        return(-1);
    }
    catch(quicky_exception::quicky_logic_exception & e)
    {
        std::cout << "ERROR : " << e.what() << std::endl ;
        return(-1);
    }
    return 0;

}

//-----------------------------------------------------------------------------
template <typename SIMPLEX_TYPE>
void test_case1(void)
{
    // Example
    // Max Z -1 * X1 - 4 * X2 - 3 * X3           = 0
    //        2 * X1 + 2 * X2 + 1 * X3 + X4      = 4
    //            X1 + 2 * X2 + 2 * X3 +    + X5 = 6
    SIMPLEX_TYPE l_simplex(5, // Number of variables : x1 and x2
                           0, // Number of inequations with the form A x <= b
                           2, // Number of equations with the form A x = b
                           0  // Number of inequations with the form A x >= b
                          );
    l_simplex.set_Z_coef(0,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_Z_coef(1,(typename SIMPLEX_TYPE::t_coef_type)4);
    l_simplex.set_Z_coef(2,(typename SIMPLEX_TYPE::t_coef_type)3);
    l_simplex.set_B_coef(0,(typename SIMPLEX_TYPE::t_coef_type)4);
    l_simplex.set_B_coef(1,(typename SIMPLEX_TYPE::t_coef_type)6);
    l_simplex.set_A_coef(0,0,(typename SIMPLEX_TYPE::t_coef_type)2);
    l_simplex.set_A_coef(0,1,(typename SIMPLEX_TYPE::t_coef_type)2);
    l_simplex.set_A_coef(0,2,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef(0,3,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef(1,0,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef(1,1,(typename SIMPLEX_TYPE::t_coef_type)2);
    l_simplex.set_A_coef(1,2,(typename SIMPLEX_TYPE::t_coef_type)2);
    l_simplex.set_A_coef(1,4,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.define_equation_type(0,simplex::t_equation_type::EQUATION);
    l_simplex.define_equation_type(1,simplex::t_equation_type::EQUATION);
    l_simplex.define_base_variable(3);
    l_simplex.define_base_variable(4);

    typename SIMPLEX_TYPE::t_coef_type l_max = (typename SIMPLEX_TYPE::t_coef_type)0;
    bool l_infinite = false;
    simplex::simplex_listener<typename SIMPLEX_TYPE::t_coef_type,typename SIMPLEX_TYPE::t_array_type> l_listener(l_simplex);
    if(l_simplex.find_max(l_max,l_infinite,&l_listener))
    {
        std::cout << "Max = " << l_max << std::endl ;
        for(unsigned int l_index = 0;
            l_index < l_simplex.get_total_nb_equation();
            ++l_index
           )
        {
            std::cout << "Base variable[" << l_index << "] is X" << l_simplex.get_base_variable(l_index) + 1 << std::endl;
        }
    }
    else if(l_infinite)
    {
        std::cout << "Inifinite Max" << std::endl;
    }
    else
    {
        std::cout << "No Max found !?" << std::endl;
    }
    assert(l_max == 10);
}

//-----------------------------------------------------------------------------
template <typename SIMPLEX_TYPE>
void test_case2(void)
{
    // Example
    // Max z = 1000 x1 + 1200 x2
    // 10 x1 + 5 x2 <= 200
    // 2 x1 + 3 x2 <= 60
    // x1 <= 34
    // x2 <= 14
    // x1 x2 >= 0

    // Manual search
    {
        unsigned int l_max = 0;
        for(unsigned int l_x1 = 0;
            l_x1 <= 34;
            ++l_x1
           )
        {
            for(unsigned int l_x2 = 0;
                l_x2 < 14;
                ++l_x2
               )
            {
                if(2 * l_x1 + 3 * l_x2 <= 60 && 10 * l_x1 + 5 * l_x2 <=200)
                {
                    unsigned int l_result = 1000 * l_x1 + 1200 * l_x2;
                    if(l_result > l_max)
                    {
                        l_max = l_result;
                        std::cout << "(" << l_x1 << "," << l_x2 << ") = " << l_max << std::endl;
                    }
                }
            }
        }
    }

    SIMPLEX_TYPE l_simplex(2, // Number of variables : x1 and x2
                           4, // Number of inequations with the form A x <= b
                           0, // Number of equations with the form A x = b
                           0  // Number of inequations with the form A x >= b
                          );
    l_simplex.set_Z_coef(0,(typename SIMPLEX_TYPE::t_coef_type)1000);
    l_simplex.set_Z_coef(1,(typename SIMPLEX_TYPE::t_coef_type)1200);
    l_simplex.set_B_coef(0,(typename SIMPLEX_TYPE::t_coef_type)200);
    l_simplex.set_B_coef(1,(typename SIMPLEX_TYPE::t_coef_type)60);
    l_simplex.set_B_coef(2,(typename SIMPLEX_TYPE::t_coef_type)34);
    l_simplex.set_B_coef(3,(typename SIMPLEX_TYPE::t_coef_type)14);
    l_simplex.set_A_coef(0,0,(typename SIMPLEX_TYPE::t_coef_type)10);
    l_simplex.set_A_coef(0,1,(typename SIMPLEX_TYPE::t_coef_type)5);
    l_simplex.set_A_coef(1,0,(typename SIMPLEX_TYPE::t_coef_type)2);
    l_simplex.set_A_coef(1,1,(typename SIMPLEX_TYPE::t_coef_type)3);
    l_simplex.set_A_coef(2,0,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef(3,1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.define_equation_type(0,simplex::t_equation_type::INEQUATION_LT);
    l_simplex.define_equation_type(1,simplex::t_equation_type::INEQUATION_LT);
    l_simplex.define_equation_type(2,simplex::t_equation_type::INEQUATION_LT);
    l_simplex.define_equation_type(3,simplex::t_equation_type::INEQUATION_LT);

    typename SIMPLEX_TYPE::t_coef_type l_max = 0;
    bool l_infinite = false;
    if(l_simplex.find_max(l_max,l_infinite))
    {
        std::cout << "Max = " << l_max << std::endl ;
    }
    else if(l_infinite)
    {
        std::cout << "Inifinite Max" << std::endl;
    }
    else
    {
        std::cout << "No Max found !?" << std::endl;
    }
    assert(27000 == l_max );
}

// 

//   2      3     1
//  A+B  | B+C | C+D

//   X2     X6   X7

// 1B      3B    2B
// D+C |  C+B  | B+A

// X10     X15   X17

typedef enum class variable
{
    P1_0=0, // 1
    P2_0,   // 2
    P3_0,   // 3
    P1_1,   // 4
    P2_1,   // 5
    P3_1,   // 6
    P1_2,   // 7
    P2_2,   // 8
    P3_2,   // 9
    P1b_0,  // 10
    P2b_0,  // 11
    P3b_0,  // 12
    P1b_1,  // 13
    P2b_1,  // 14
    P3b_1,  // 15
    P1b_2,  // 16
    P2b_2,  // 17
    P3b_2   // 18
} t_variable;

const std::ostream & display(std::ostream & p_stream,const t_variable & p_variable)
{
    switch(p_variable)
    {
        case t_variable::P1_0:
            p_stream << "P1_0" ;
            break;
        case t_variable::P2_0:
            p_stream << "P2_0" ;
            break;
        case t_variable::P3_0:
            p_stream << "P3_0" ;
            break;
        case t_variable::P1_1:
            p_stream << "P1_1" ;
            break;
        case t_variable::P2_1:
            p_stream << "P2_1" ;
            break;
        case t_variable::P3_1:
            p_stream << "P3_1" ;
            break;
        case t_variable::P1_2:
            p_stream << "P1_2" ;
            break;
        case t_variable::P2_2:
            p_stream << "P2_2" ;
            break;
        case t_variable::P3_2:
            p_stream << "P3_2" ;
            break;
        case t_variable::P1b_0:
            p_stream << "P1b_0" ;
            break;
        case t_variable::P2b_0:
            p_stream << "P2b_0" ;
            break;
        case t_variable::P3b_0:
            p_stream << "P3b_0" ;
            break;
        case t_variable::P1b_1:
            p_stream << "P1b_1" ;
            break;
        case t_variable::P2b_1:
            p_stream << "P2b_1" ;
            break;
        case t_variable::P3b_1:
            p_stream << "P3b_1" ;
            break;
        case t_variable::P1b_2:
            p_stream << "P1b_2" ;
            break;
        case t_variable::P2b_2:
            p_stream << "P2b_2" ;
            break;
        case t_variable::P3b_2:
            p_stream << "P3b_2" ;
            break;
        default:
            throw quicky_exception::quicky_logic_exception("Unknown t_variable value : "+std::to_string((unsigned int)p_variable),__LINE__,__FILE__);
    }
    return p_stream;
}

typedef enum class equation
{
    Pos0=0,
    Pos1,
    Pos2,
    P1,
    P2,
    P3,
    P1_0_P2_1,
    P1_0_P3_1,
    P1_0_P2b_1,
    P1_0_P3b_1,
    P1b_0_P2_1,
    P1b_0_P3_1,
    P1b_0_P2b_1,
    P1b_0_P3b_1,
    P2_0_P1_1,
    P2_0_P3_1,
    P2_0_P1b_1,
    P2_0_P3b_1,
    P2b_0_P1_1,
    P2b_0_P3_1,
    P2b_0_P1b_1,
    P2b_0_P3b_1,
    P3_0_P1_1,
    P3_0_P2_1,
    P3_0_P1b_1,
    P3_0_P2b_1,
    P3b_0_P1_1,
    P3b_0_P2_1,
    P3b_0_P1b_1,
    P3b_0_P2b_1,
    P1_1_P2_2,
    P1_1_P3_2,
    P1_1_P2b_2,
    P1_1_P3b_2,
    P1b_1_P2_2,
    P1b_1_P3_2,
    P1b_1_P2b_2,
    P1b_1_P3b_2,
    P2_1_P1_2,
    P2_1_P3_2,
    P2_1_P1b_2,
    P2_1_P3b_2,
    P2b_1_P1_2,
    P2b_1_P3_2,
    P2b_1_P1b_2,
    P2b_1_P3b_2,
    P3_1_P1_2,
    P3_1_P2_2,
    P3_1_P1b_2,
    P3_1_P2b_2,
    P3b_1_P1_2,
    P3b_1_P2_2,
    P3b_1_P1b_2,
    P3b_1_P2b_2
} t_equation;

typedef enum class yes_no
{
    NO=1,
    YES
} t_yes_no;

template <typename SIMPLEX_TYPE>
void test_case3(void)
{
    // 3 position equations
    // 3 pieces equations

    // 18 variables : 3 pieces at 3 positions with 2 orientations : 3*3*2
    // 54 equations : 48 : (3pieces * 2orientations) * (2pieces * 2orientations) * 2 segments
    SIMPLEX_TYPE l_simplex(18, // Number of variables : x1 and x2
                           54, // Number of inequations with the form A x <= b
                           0, // Number of equations with the form A x = b
                           0  // Number of inequations with the form A x >= b
                          );
    for(unsigned int l_index = 0;
        l_index < 18;
        ++l_index
            )
    {
        l_simplex.set_Z_coef(l_index,(typename SIMPLEX_TYPE::t_coef_type)1);
    }

    // Position equations
    l_simplex.set_A_coef((unsigned int)t_equation::Pos0,
                         (unsigned int)t_variable::P1_0,
                         (typename SIMPLEX_TYPE::t_coef_type)1
                        );
    l_simplex.set_A_coef((unsigned int)t_equation::Pos0,
                         (unsigned int)t_variable::P2_0,
                         (typename SIMPLEX_TYPE::t_coef_type)1
                        );
    l_simplex.set_A_coef((unsigned int)t_equation::Pos0,
                         (unsigned int)t_variable::P3_0,
                         (typename SIMPLEX_TYPE::t_coef_type)1
                        );
    l_simplex.set_A_coef((unsigned int)t_equation::Pos0,
                         (unsigned int)t_variable::P1b_0,
                         (typename SIMPLEX_TYPE::t_coef_type)1
                        );
    l_simplex.set_A_coef((unsigned int)t_equation::Pos0,
                         (unsigned int)t_variable::P2b_0,
                         (typename SIMPLEX_TYPE::t_coef_type)1
                        );
    l_simplex.set_A_coef((unsigned int)t_equation::Pos0,
                         (unsigned int)t_variable::P3b_0,
                         (typename SIMPLEX_TYPE::t_coef_type)1
                        );

    l_simplex.set_B_coef((unsigned int)t_equation::Pos0,
                         (typename SIMPLEX_TYPE::t_coef_type)1
                        );
    l_simplex.define_equation_type((unsigned int)t_equation::Pos0,
                                   simplex::t_equation_type::INEQUATION_LT
                                  );


    l_simplex.set_A_coef((unsigned int)t_equation::Pos1,
                         (unsigned int)t_variable::P1_1,
                         (typename SIMPLEX_TYPE::t_coef_type)1
                        );
    l_simplex.set_A_coef((unsigned int)t_equation::Pos1,
                         (unsigned int)t_variable::P2_1,
                         (typename SIMPLEX_TYPE::t_coef_type)1
                        );
    l_simplex.set_A_coef((unsigned int)t_equation::Pos1,
                         (unsigned int)t_variable::P3_1,
                         (typename SIMPLEX_TYPE::t_coef_type)1
                        );
    l_simplex.set_A_coef((unsigned int)t_equation::Pos1,
                         (unsigned int)t_variable::P1b_1,
                         (typename SIMPLEX_TYPE::t_coef_type)1
                        );
    l_simplex.set_A_coef((unsigned int)t_equation::Pos1,
                         (unsigned int)t_variable::P2b_1,
                         (typename SIMPLEX_TYPE::t_coef_type)1
                        );
    l_simplex.set_A_coef((unsigned int)t_equation::Pos1,
                         (unsigned int)t_variable::P3b_1,
                         (typename SIMPLEX_TYPE::t_coef_type)1
                        );

    l_simplex.set_B_coef((unsigned int)t_equation::Pos1,
                         (typename SIMPLEX_TYPE::t_coef_type)1
                        );
    l_simplex.define_equation_type((unsigned int)t_equation::Pos1,
                                   simplex::t_equation_type::INEQUATION_LT
                                  );

    l_simplex.set_A_coef((unsigned int)t_equation::Pos2,
                         (unsigned int)t_variable::P1_2,
                         (typename SIMPLEX_TYPE::t_coef_type)1
                        );
    l_simplex.set_A_coef((unsigned int)t_equation::Pos2,
                         (unsigned int)t_variable::P2_2,
                         (typename SIMPLEX_TYPE::t_coef_type)1
                        );
    l_simplex.set_A_coef((unsigned int)t_equation::Pos2,
                         (unsigned int)t_variable::P3_2,
                         (typename SIMPLEX_TYPE::t_coef_type)1
                        );
    l_simplex.set_A_coef((unsigned int)t_equation::Pos2,
                         (unsigned int)t_variable::P1b_2,
                         (typename SIMPLEX_TYPE::t_coef_type)1
                        );
    l_simplex.set_A_coef((unsigned int)t_equation::Pos2,
                         (unsigned int)t_variable::P2b_2,
                         (typename SIMPLEX_TYPE::t_coef_type)1
                        );
    l_simplex.set_A_coef((unsigned int)t_equation::Pos2,
                         (unsigned int)t_variable::P3b_2,
                         (typename SIMPLEX_TYPE::t_coef_type)1
                        );

    l_simplex.set_B_coef((unsigned int)t_equation::Pos2,
                         (typename SIMPLEX_TYPE::t_coef_type)1
                        );
    l_simplex.define_equation_type((unsigned int)t_equation::Pos2,
                                   simplex::t_equation_type::INEQUATION_LT
                                  );

    // Pieces equation
    l_simplex.set_A_coef((unsigned int)t_equation::P1,
                         (unsigned int)t_variable::P1_0,
                         (typename SIMPLEX_TYPE::t_coef_type)1
                        );
    l_simplex.set_A_coef((unsigned int)t_equation::P1,
                         (unsigned int)t_variable::P1_1,
                         (typename SIMPLEX_TYPE::t_coef_type)1
                        );
    l_simplex.set_A_coef((unsigned int)t_equation::P1,
                         (unsigned int)t_variable::P1_2,
                         (typename SIMPLEX_TYPE::t_coef_type)1
                        );
    l_simplex.set_A_coef((unsigned int)t_equation::P1,
                         (unsigned int)t_variable::P1b_0,
                         (typename SIMPLEX_TYPE::t_coef_type)1
                        );
    l_simplex.set_A_coef((unsigned int)t_equation::P1,
                         (unsigned int)t_variable::P1b_1,
                         (typename SIMPLEX_TYPE::t_coef_type)1
                        );
    l_simplex.set_A_coef((unsigned int)t_equation::P1,
                         (unsigned int)t_variable::P1b_2,
                         (typename SIMPLEX_TYPE::t_coef_type)1
                        );

    l_simplex.set_B_coef((unsigned int)t_equation::P1,
                         (typename SIMPLEX_TYPE::t_coef_type)1
                        );
    l_simplex.define_equation_type((unsigned int)t_equation::P1,
                                   simplex::t_equation_type::INEQUATION_LT
                                  );

    l_simplex.set_A_coef((unsigned int)t_equation::P2,
                         (unsigned int)t_variable::P2_0,
                         (typename SIMPLEX_TYPE::t_coef_type)1
                        );
    l_simplex.set_A_coef((unsigned int)t_equation::P2,
                         (unsigned int)t_variable::P2_1,
                         (typename SIMPLEX_TYPE::t_coef_type)1
                        );
    l_simplex.set_A_coef((unsigned int)t_equation::P2,
                         (unsigned int)t_variable::P2_2,
                         (typename SIMPLEX_TYPE::t_coef_type)1
                        );
    l_simplex.set_A_coef((unsigned int)t_equation::P2,
                         (unsigned int)t_variable::P2b_0,
                         (typename SIMPLEX_TYPE::t_coef_type)1
                        );
    l_simplex.set_A_coef((unsigned int)t_equation::P2,
                         (unsigned int)t_variable::P2b_1,
                         (typename SIMPLEX_TYPE::t_coef_type)1
                        );
    l_simplex.set_A_coef((unsigned int)t_equation::P2,
                         (unsigned int)t_variable::P2b_2,
                         (typename SIMPLEX_TYPE::t_coef_type)1
                        );

    l_simplex.set_B_coef((unsigned int)t_equation::P2,
                         (typename SIMPLEX_TYPE::t_coef_type)1
                        );
    l_simplex.define_equation_type((unsigned int)t_equation::P2,
                                   simplex::t_equation_type::INEQUATION_LT
                                  );

    l_simplex.set_A_coef((unsigned int)t_equation::P3,
                         (unsigned int)t_variable::P3_0,
                         (typename SIMPLEX_TYPE::t_coef_type)1
                        );
    l_simplex.set_A_coef((unsigned int)t_equation::P3,
                         (unsigned int)t_variable::P3_1,
                         (typename SIMPLEX_TYPE::t_coef_type)1
                        );
    l_simplex.set_A_coef((unsigned int)t_equation::P3,
                         (unsigned int)t_variable::P3_2,
                         (typename SIMPLEX_TYPE::t_coef_type)1
                        );
    l_simplex.set_A_coef((unsigned int)t_equation::P3,
                         (unsigned int)t_variable::P3b_0,
                         (typename SIMPLEX_TYPE::t_coef_type)1
                        );
    l_simplex.set_A_coef((unsigned int)t_equation::P3,
                         (unsigned int)t_variable::P3b_1,
                         (typename SIMPLEX_TYPE::t_coef_type)1
                        );
    l_simplex.set_A_coef((unsigned int)t_equation::P3,
                         (unsigned int)t_variable::P3b_2,
                         (typename SIMPLEX_TYPE::t_coef_type)1
                        );

    l_simplex.set_B_coef((unsigned int)t_equation::P3,
                         (typename SIMPLEX_TYPE::t_coef_type)1
                        );
    l_simplex.define_equation_type((unsigned int)t_equation::P3,
                                   simplex::t_equation_type::INEQUATION_LT
                                  );

    // Combination equations

    l_simplex.set_A_coef((unsigned int)t_equation::P1_0_P2_1,(unsigned int)t_variable::P1_0,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef((unsigned int)t_equation::P1_0_P2_1,(unsigned int)t_variable::P2_1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_B_coef((unsigned int)t_equation::P1_0_P2_1,(typename SIMPLEX_TYPE::t_coef_type)t_yes_no::NO);
    l_simplex.define_equation_type((unsigned int)t_equation::P1_0_P2_1,simplex::t_equation_type::INEQUATION_LT);

    l_simplex.set_A_coef((unsigned int)t_equation::P1_0_P3_1,(unsigned int)t_variable::P1_0,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef((unsigned int)t_equation::P1_0_P3_1,(unsigned int)t_variable::P3_1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_B_coef((unsigned int)t_equation::P1_0_P3_1,(typename SIMPLEX_TYPE::t_coef_type)t_yes_no::NO);
    l_simplex.define_equation_type((unsigned int)t_equation::P1_0_P3_1,simplex::t_equation_type::INEQUATION_LT);

    l_simplex.set_A_coef((unsigned int)t_equation::P1_0_P2b_1,(unsigned int)t_variable::P1_0,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef((unsigned int)t_equation::P1_0_P2b_1,(unsigned int)t_variable::P2b_1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_B_coef((unsigned int)t_equation::P1_0_P2b_1,(typename SIMPLEX_TYPE::t_coef_type)t_yes_no::NO);
    l_simplex.define_equation_type((unsigned int)t_equation::P1_0_P2b_1,simplex::t_equation_type::INEQUATION_LT);

    l_simplex.set_A_coef((unsigned int)t_equation::P1_0_P3b_1,(unsigned int)t_variable::P1_0,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef((unsigned int)t_equation::P1_0_P3b_1,(unsigned int)t_variable::P3b_1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_B_coef((unsigned int)t_equation::P1_0_P3b_1,(typename SIMPLEX_TYPE::t_coef_type)t_yes_no::NO);
    l_simplex.define_equation_type((unsigned int)t_equation::P1_0_P3b_1,simplex::t_equation_type::INEQUATION_LT);

    l_simplex.set_A_coef((unsigned int)t_equation::P1b_0_P2_1,(unsigned int)t_variable::P1b_0,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef((unsigned int)t_equation::P1b_0_P2_1,(unsigned int)t_variable::P2_1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_B_coef((unsigned int)t_equation::P1b_0_P2_1,(typename SIMPLEX_TYPE::t_coef_type)t_yes_no::NO);
    l_simplex.define_equation_type((unsigned int)t_equation::P1b_0_P2_1,simplex::t_equation_type::INEQUATION_LT);

    l_simplex.set_A_coef((unsigned int)t_equation::P1b_0_P3_1,(unsigned int)t_variable::P1b_0,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef((unsigned int)t_equation::P1b_0_P3_1,(unsigned int)t_variable::P3_1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_B_coef((unsigned int)t_equation::P1b_0_P3_1,(typename SIMPLEX_TYPE::t_coef_type)t_yes_no::NO);
    l_simplex.define_equation_type((unsigned int)t_equation::P1b_0_P3_1,simplex::t_equation_type::INEQUATION_LT);

    l_simplex.set_A_coef((unsigned int)t_equation::P1b_0_P2b_1,(unsigned int)t_variable::P1b_0,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef((unsigned int)t_equation::P1b_0_P2b_1,(unsigned int)t_variable::P2b_1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_B_coef((unsigned int)t_equation::P1b_0_P2b_1,(typename SIMPLEX_TYPE::t_coef_type)t_yes_no::NO);
    l_simplex.define_equation_type((unsigned int)t_equation::P1b_0_P2b_1,simplex::t_equation_type::INEQUATION_LT);

    l_simplex.set_A_coef((unsigned int)t_equation::P1b_0_P3b_1,(unsigned int)t_variable::P1b_0,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef((unsigned int)t_equation::P1b_0_P3b_1,(unsigned int)t_variable::P3b_1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_B_coef((unsigned int)t_equation::P1b_0_P3b_1,(typename SIMPLEX_TYPE::t_coef_type)t_yes_no::YES);
    l_simplex.define_equation_type((unsigned int)t_equation::P1b_0_P3b_1,simplex::t_equation_type::INEQUATION_LT);

    l_simplex.set_A_coef((unsigned int)t_equation::P2_0_P1_1,(unsigned int)t_variable::P2_0,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef((unsigned int)t_equation::P2_0_P1_1,(unsigned int)t_variable::P1_1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_B_coef((unsigned int)t_equation::P2_0_P1_1,(typename SIMPLEX_TYPE::t_coef_type)t_yes_no::NO);
    l_simplex.define_equation_type((unsigned int)t_equation::P2_0_P1_1,simplex::t_equation_type::INEQUATION_LT);

    l_simplex.set_A_coef((unsigned int)t_equation::P2_0_P3_1,(unsigned int)t_variable::P2_0,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef((unsigned int)t_equation::P2_0_P3_1,(unsigned int)t_variable::P3_1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_B_coef((unsigned int)t_equation::P2_0_P3_1,(typename SIMPLEX_TYPE::t_coef_type)t_yes_no::YES);
    l_simplex.define_equation_type((unsigned int)t_equation::P2_0_P3_1,simplex::t_equation_type::INEQUATION_LT);

    l_simplex.set_A_coef((unsigned int)t_equation::P2_0_P1b_1,(unsigned int)t_variable::P2_0,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef((unsigned int)t_equation::P2_0_P1b_1,(unsigned int)t_variable::P1b_1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_B_coef((unsigned int)t_equation::P2_0_P1b_1,(typename SIMPLEX_TYPE::t_coef_type)t_yes_no::NO);
    l_simplex.define_equation_type((unsigned int)t_equation::P2_0_P1b_1,simplex::t_equation_type::INEQUATION_LT);

    l_simplex.set_A_coef((unsigned int)t_equation::P2_0_P3b_1,(unsigned int)t_variable::P2_0,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef((unsigned int)t_equation::P2_0_P3b_1,(unsigned int)t_variable::P3b_1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_B_coef((unsigned int)t_equation::P2_0_P3b_1,(typename SIMPLEX_TYPE::t_coef_type)t_yes_no::NO);
    l_simplex.define_equation_type((unsigned int)t_equation::P2_0_P3b_1,simplex::t_equation_type::INEQUATION_LT);

    l_simplex.set_A_coef((unsigned int)t_equation::P2b_0_P1_1,(unsigned int)t_variable::P2b_0,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef((unsigned int)t_equation::P2b_0_P1_1,(unsigned int)t_variable::P1_1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_B_coef((unsigned int)t_equation::P2b_0_P1_1,(typename SIMPLEX_TYPE::t_coef_type)t_yes_no::NO);
    l_simplex.define_equation_type((unsigned int)t_equation::P2b_0_P1_1,simplex::t_equation_type::INEQUATION_LT);

    l_simplex.set_A_coef((unsigned int)t_equation::P2b_0_P3_1,(unsigned int)t_variable::P2b_0,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef((unsigned int)t_equation::P2b_0_P3_1,(unsigned int)t_variable::P3_1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_B_coef((unsigned int)t_equation::P2b_0_P3_1,(typename SIMPLEX_TYPE::t_coef_type)t_yes_no::NO);
    l_simplex.define_equation_type((unsigned int)t_equation::P2b_0_P3_1,simplex::t_equation_type::INEQUATION_LT);

    l_simplex.set_A_coef((unsigned int)t_equation::P2b_0_P1b_1,(unsigned int)t_variable::P2b_0,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef((unsigned int)t_equation::P2b_0_P1b_1,(unsigned int)t_variable::P1b_1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_B_coef((unsigned int)t_equation::P2b_0_P1b_1,(typename SIMPLEX_TYPE::t_coef_type)t_yes_no::NO);
    l_simplex.define_equation_type((unsigned int)t_equation::P2b_0_P1b_1,simplex::t_equation_type::INEQUATION_LT);

    l_simplex.set_A_coef((unsigned int)t_equation::P2b_0_P3b_1,(unsigned int)t_variable::P2b_0,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef((unsigned int)t_equation::P2b_0_P3b_1,(unsigned int)t_variable::P3b_1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_B_coef((unsigned int)t_equation::P2b_0_P3b_1,(typename SIMPLEX_TYPE::t_coef_type)t_yes_no::NO);
    l_simplex.define_equation_type((unsigned int)t_equation::P2b_0_P3b_1,simplex::t_equation_type::INEQUATION_LT);

    l_simplex.set_A_coef((unsigned int)t_equation::P3_0_P1_1,(unsigned int)t_variable::P3_0,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef((unsigned int)t_equation::P3_0_P1_1,(unsigned int)t_variable::P1_1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_B_coef((unsigned int)t_equation::P3_0_P1_1,(typename SIMPLEX_TYPE::t_coef_type)t_yes_no::YES);
    l_simplex.define_equation_type((unsigned int)t_equation::P3_0_P1_1,simplex::t_equation_type::INEQUATION_LT);

    l_simplex.set_A_coef((unsigned int)t_equation::P3_0_P2_1,(unsigned int)t_variable::P3_0,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef((unsigned int)t_equation::P3_0_P2_1,(unsigned int)t_variable::P2_1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_B_coef((unsigned int)t_equation::P3_0_P2_1,(typename SIMPLEX_TYPE::t_coef_type)t_yes_no::NO);
    l_simplex.define_equation_type((unsigned int)t_equation::P3_0_P2_1,simplex::t_equation_type::INEQUATION_LT);

    l_simplex.set_A_coef((unsigned int)t_equation::P3_0_P1b_1,(unsigned int)t_variable::P3_0,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef((unsigned int)t_equation::P3_0_P1b_1,(unsigned int)t_variable::P1b_1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_B_coef((unsigned int)t_equation::P3_0_P1b_1,(typename SIMPLEX_TYPE::t_coef_type)t_yes_no::NO);
    l_simplex.define_equation_type((unsigned int)t_equation::P3_0_P1b_1,simplex::t_equation_type::INEQUATION_LT);

    l_simplex.set_A_coef((unsigned int)t_equation::P3_0_P2b_1,(unsigned int)t_variable::P3_0,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef((unsigned int)t_equation::P3_0_P2b_1,(unsigned int)t_variable::P2b_1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_B_coef((unsigned int)t_equation::P3_0_P2b_1,(typename SIMPLEX_TYPE::t_coef_type)t_yes_no::NO);
    l_simplex.define_equation_type((unsigned int)t_equation::P3_0_P2b_1,simplex::t_equation_type::INEQUATION_LT);

    l_simplex.set_A_coef((unsigned int)t_equation::P3b_0_P1_1,(unsigned int)t_variable::P3b_0,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef((unsigned int)t_equation::P3b_0_P1_1,(unsigned int)t_variable::P1_1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_B_coef((unsigned int)t_equation::P3b_0_P1_1,(typename SIMPLEX_TYPE::t_coef_type)t_yes_no::NO);
    l_simplex.define_equation_type((unsigned int)t_equation::P3b_0_P1_1,simplex::t_equation_type::INEQUATION_LT);

    l_simplex.set_A_coef((unsigned int)t_equation::P3b_0_P2_1,(unsigned int)t_variable::P3b_0,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef((unsigned int)t_equation::P3b_0_P2_1,(unsigned int)t_variable::P2_1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_B_coef((unsigned int)t_equation::P3b_0_P2_1,(typename SIMPLEX_TYPE::t_coef_type)t_yes_no::NO);
    l_simplex.define_equation_type((unsigned int)t_equation::P3b_0_P2_1,simplex::t_equation_type::INEQUATION_LT);

    l_simplex.set_A_coef((unsigned int)t_equation::P3b_0_P1b_1,(unsigned int)t_variable::P3b_0,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef((unsigned int)t_equation::P3b_0_P1b_1,(unsigned int)t_variable::P1b_1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_B_coef((unsigned int)t_equation::P3b_0_P1b_1,(typename SIMPLEX_TYPE::t_coef_type)t_yes_no::NO);
    l_simplex.define_equation_type((unsigned int)t_equation::P3b_0_P1b_1,simplex::t_equation_type::INEQUATION_LT);

    l_simplex.set_A_coef((unsigned int)t_equation::P3b_0_P2b_1,(unsigned int)t_variable::P3b_0,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef((unsigned int)t_equation::P3b_0_P2b_1,(unsigned int)t_variable::P2b_1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_B_coef((unsigned int)t_equation::P3b_0_P2b_1,(typename SIMPLEX_TYPE::t_coef_type)t_yes_no::YES);
    l_simplex.define_equation_type((unsigned int)t_equation::P3b_0_P2b_1,simplex::t_equation_type::INEQUATION_LT);

    l_simplex.set_A_coef((unsigned int)t_equation::P1_1_P2_2,(unsigned int)t_variable::P1_1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef((unsigned int)t_equation::P1_1_P2_2,(unsigned int)t_variable::P2_2,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_B_coef((unsigned int)t_equation::P1_1_P2_2,(typename SIMPLEX_TYPE::t_coef_type)t_yes_no::NO);
    l_simplex.define_equation_type((unsigned int)t_equation::P1_1_P2_2,simplex::t_equation_type::INEQUATION_LT);

    l_simplex.set_A_coef((unsigned int)t_equation::P1_1_P3_2,(unsigned int)t_variable::P1_1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef((unsigned int)t_equation::P1_1_P3_2,(unsigned int)t_variable::P3_2,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_B_coef((unsigned int)t_equation::P1_1_P3_2,(typename SIMPLEX_TYPE::t_coef_type)t_yes_no::NO);
    l_simplex.define_equation_type((unsigned int)t_equation::P1_1_P3_2,simplex::t_equation_type::INEQUATION_LT);

    l_simplex.set_A_coef((unsigned int)t_equation::P1_1_P2b_2,(unsigned int)t_variable::P1_1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef((unsigned int)t_equation::P1_1_P2b_2,(unsigned int)t_variable::P2b_2,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_B_coef((unsigned int)t_equation::P1_1_P2b_2,(typename SIMPLEX_TYPE::t_coef_type)t_yes_no::NO);
    l_simplex.define_equation_type((unsigned int)t_equation::P1_1_P2b_2,simplex::t_equation_type::INEQUATION_LT);

    l_simplex.set_A_coef((unsigned int)t_equation::P1_1_P3b_2,(unsigned int)t_variable::P1_1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef((unsigned int)t_equation::P1_1_P3b_2,(unsigned int)t_variable::P3b_2,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_B_coef((unsigned int)t_equation::P1_1_P3b_2,(typename SIMPLEX_TYPE::t_coef_type)t_yes_no::NO);
    l_simplex.define_equation_type((unsigned int)t_equation::P1_1_P3b_2,simplex::t_equation_type::INEQUATION_LT);

    l_simplex.set_A_coef((unsigned int)t_equation::P1b_1_P2_2,(unsigned int)t_variable::P1b_1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef((unsigned int)t_equation::P1b_1_P2_2,(unsigned int)t_variable::P2_2,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_B_coef((unsigned int)t_equation::P1b_1_P2_2,(typename SIMPLEX_TYPE::t_coef_type)t_yes_no::NO);
    l_simplex.define_equation_type((unsigned int)t_equation::P1b_1_P2_2,simplex::t_equation_type::INEQUATION_LT);

    l_simplex.set_A_coef((unsigned int)t_equation::P1b_1_P3_2,(unsigned int)t_variable::P1b_1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef((unsigned int)t_equation::P1b_1_P3_2,(unsigned int)t_variable::P3_2,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_B_coef((unsigned int)t_equation::P1b_1_P3_2,(typename SIMPLEX_TYPE::t_coef_type)t_yes_no::NO);
    l_simplex.define_equation_type((unsigned int)t_equation::P1b_1_P3_2,simplex::t_equation_type::INEQUATION_LT);

    l_simplex.set_A_coef((unsigned int)t_equation::P1b_1_P2b_2,(unsigned int)t_variable::P1b_1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef((unsigned int)t_equation::P1b_1_P2b_2,(unsigned int)t_variable::P2b_2,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_B_coef((unsigned int)t_equation::P1b_1_P2b_2,(typename SIMPLEX_TYPE::t_coef_type)t_yes_no::NO);
    l_simplex.define_equation_type((unsigned int)t_equation::P1b_1_P2b_2,simplex::t_equation_type::INEQUATION_LT);

    l_simplex.set_A_coef((unsigned int)t_equation::P1b_1_P3b_2,(unsigned int)t_variable::P1b_1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef((unsigned int)t_equation::P1b_1_P3b_2,(unsigned int)t_variable::P3b_2,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_B_coef((unsigned int)t_equation::P1b_1_P3b_2,(typename SIMPLEX_TYPE::t_coef_type)t_yes_no::YES);
    l_simplex.define_equation_type((unsigned int)t_equation::P1b_1_P3b_2,simplex::t_equation_type::INEQUATION_LT);

    l_simplex.set_A_coef((unsigned int)t_equation::P2_1_P1_2,(unsigned int)t_variable::P2_1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef((unsigned int)t_equation::P2_1_P1_2,(unsigned int)t_variable::P1_2,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_B_coef((unsigned int)t_equation::P2_1_P1_2,(typename SIMPLEX_TYPE::t_coef_type)t_yes_no::NO);
    l_simplex.define_equation_type((unsigned int)t_equation::P2_1_P1_2,simplex::t_equation_type::INEQUATION_LT);

    l_simplex.set_A_coef((unsigned int)t_equation::P2_1_P3_2,(unsigned int)t_variable::P2_1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef((unsigned int)t_equation::P2_1_P3_2,(unsigned int)t_variable::P3_2,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_B_coef((unsigned int)t_equation::P2_1_P3_2,(typename SIMPLEX_TYPE::t_coef_type)t_yes_no::YES);
    l_simplex.define_equation_type((unsigned int)t_equation::P2_1_P3_2,simplex::t_equation_type::INEQUATION_LT);

    l_simplex.set_A_coef((unsigned int)t_equation::P2_1_P1b_2,(unsigned int)t_variable::P2_1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef((unsigned int)t_equation::P2_1_P1b_2,(unsigned int)t_variable::P1b_2,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_B_coef((unsigned int)t_equation::P2_1_P1b_2,(typename SIMPLEX_TYPE::t_coef_type)t_yes_no::NO);
    l_simplex.define_equation_type((unsigned int)t_equation::P2_1_P1b_2,simplex::t_equation_type::INEQUATION_LT);

    l_simplex.set_A_coef((unsigned int)t_equation::P2_1_P3b_2,(unsigned int)t_variable::P2_1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef((unsigned int)t_equation::P2_1_P3b_2,(unsigned int)t_variable::P3b_2,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_B_coef((unsigned int)t_equation::P2_1_P3b_2,(typename SIMPLEX_TYPE::t_coef_type)t_yes_no::NO);
    l_simplex.define_equation_type((unsigned int)t_equation::P2_1_P3b_2,simplex::t_equation_type::INEQUATION_LT);

    l_simplex.set_A_coef((unsigned int)t_equation::P2b_1_P1_2,(unsigned int)t_variable::P2b_1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef((unsigned int)t_equation::P2b_1_P1_2,(unsigned int)t_variable::P1_2,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_B_coef((unsigned int)t_equation::P2b_1_P1_2,(typename SIMPLEX_TYPE::t_coef_type)t_yes_no::NO);
    l_simplex.define_equation_type((unsigned int)t_equation::P2b_1_P1_2,simplex::t_equation_type::INEQUATION_LT);

    l_simplex.set_A_coef((unsigned int)t_equation::P2b_1_P3_2,(unsigned int)t_variable::P2b_1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef((unsigned int)t_equation::P2b_1_P3_2,(unsigned int)t_variable::P3_2,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_B_coef((unsigned int)t_equation::P2b_1_P3_2,(typename SIMPLEX_TYPE::t_coef_type)t_yes_no::NO);
    l_simplex.define_equation_type((unsigned int)t_equation::P2b_1_P3_2,simplex::t_equation_type::INEQUATION_LT);

    l_simplex.set_A_coef((unsigned int)t_equation::P2b_1_P1b_2,(unsigned int)t_variable::P2b_1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef((unsigned int)t_equation::P2b_1_P1b_2,(unsigned int)t_variable::P1b_2,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_B_coef((unsigned int)t_equation::P2b_1_P1b_2,(typename SIMPLEX_TYPE::t_coef_type)t_yes_no::NO);
    l_simplex.define_equation_type((unsigned int)t_equation::P2b_1_P1b_2,simplex::t_equation_type::INEQUATION_LT);

    l_simplex.set_A_coef((unsigned int)t_equation::P2b_1_P3b_2,(unsigned int)t_variable::P2b_1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef((unsigned int)t_equation::P2b_1_P3b_2,(unsigned int)t_variable::P3b_2,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_B_coef((unsigned int)t_equation::P2b_1_P3b_2,(typename SIMPLEX_TYPE::t_coef_type)t_yes_no::NO);
    l_simplex.define_equation_type((unsigned int)t_equation::P2b_1_P3b_2,simplex::t_equation_type::INEQUATION_LT);

    l_simplex.set_A_coef((unsigned int)t_equation::P3_1_P1_2,(unsigned int)t_variable::P3_1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef((unsigned int)t_equation::P3_1_P1_2,(unsigned int)t_variable::P1_2,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_B_coef((unsigned int)t_equation::P3_1_P1_2,(typename SIMPLEX_TYPE::t_coef_type)t_yes_no::YES);
    l_simplex.define_equation_type((unsigned int)t_equation::P3_1_P1_2,simplex::t_equation_type::INEQUATION_LT);

    l_simplex.set_A_coef((unsigned int)t_equation::P3_1_P2_2,(unsigned int)t_variable::P3_1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef((unsigned int)t_equation::P3_1_P2_2,(unsigned int)t_variable::P2_2,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_B_coef((unsigned int)t_equation::P3_1_P2_2,(typename SIMPLEX_TYPE::t_coef_type)t_yes_no::NO);
    l_simplex.define_equation_type((unsigned int)t_equation::P3_1_P2_2,simplex::t_equation_type::INEQUATION_LT);

    l_simplex.set_A_coef((unsigned int)t_equation::P3_1_P1b_2,(unsigned int)t_variable::P3_1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef((unsigned int)t_equation::P3_1_P1b_2,(unsigned int)t_variable::P1b_2,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_B_coef((unsigned int)t_equation::P3_1_P1b_2,(typename SIMPLEX_TYPE::t_coef_type)t_yes_no::NO);
    l_simplex.define_equation_type((unsigned int)t_equation::P3_1_P1b_2,simplex::t_equation_type::INEQUATION_LT);

    l_simplex.set_A_coef((unsigned int)t_equation::P3_1_P2b_2,(unsigned int)t_variable::P3_1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef((unsigned int)t_equation::P3_1_P2b_2,(unsigned int)t_variable::P2b_2,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_B_coef((unsigned int)t_equation::P3_1_P2b_2,(typename SIMPLEX_TYPE::t_coef_type)t_yes_no::NO);
    l_simplex.define_equation_type((unsigned int)t_equation::P3_1_P2b_2,simplex::t_equation_type::INEQUATION_LT);

    l_simplex.set_A_coef((unsigned int)t_equation::P3b_1_P1_2,(unsigned int)t_variable::P3b_1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef((unsigned int)t_equation::P3b_1_P1_2,(unsigned int)t_variable::P1_2,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_B_coef((unsigned int)t_equation::P3b_1_P1_2,(typename SIMPLEX_TYPE::t_coef_type)t_yes_no::NO);
    l_simplex.define_equation_type((unsigned int)t_equation::P3b_1_P1_2,simplex::t_equation_type::INEQUATION_LT);

    l_simplex.set_A_coef((unsigned int)t_equation::P3b_1_P2_2,(unsigned int)t_variable::P3b_1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef((unsigned int)t_equation::P3b_1_P2_2,(unsigned int)t_variable::P2_2,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_B_coef((unsigned int)t_equation::P3b_1_P2_2,(typename SIMPLEX_TYPE::t_coef_type)t_yes_no::NO);
    l_simplex.define_equation_type((unsigned int)t_equation::P3b_1_P2_2,simplex::t_equation_type::INEQUATION_LT);

    l_simplex.set_A_coef((unsigned int)t_equation::P3b_1_P1b_2,(unsigned int)t_variable::P3b_1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef((unsigned int)t_equation::P3b_1_P1b_2,(unsigned int)t_variable::P1b_2,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_B_coef((unsigned int)t_equation::P3b_1_P1b_2,(typename SIMPLEX_TYPE::t_coef_type)t_yes_no::NO);
    l_simplex.define_equation_type((unsigned int)t_equation::P3b_1_P1b_2,simplex::t_equation_type::INEQUATION_LT);

    l_simplex.set_A_coef((unsigned int)t_equation::P3b_1_P2b_2,(unsigned int)t_variable::P3b_1,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_A_coef((unsigned int)t_equation::P3b_1_P2b_2,(unsigned int)t_variable::P2b_2,(typename SIMPLEX_TYPE::t_coef_type)1);
    l_simplex.set_B_coef((unsigned int)t_equation::P3b_1_P2b_2,(typename SIMPLEX_TYPE::t_coef_type)t_yes_no::YES);
    l_simplex.define_equation_type((unsigned int)t_equation::P3b_1_P2b_2,simplex::t_equation_type::INEQUATION_LT);


    //  l_simplex.set_B_coef(0,200);
    //  l_simplex.set_B_coef(1,60);
    //  l_simplex.set_B_coef(2,34);
    //  l_simplex.set_B_coef(3,14);
    //  l_simplex.set_A_coef(0,1,5);
    //  l_simplex.set_A_coef(1,0,2);
    //  l_simplex.set_A_coef(1,1,3);
    //  l_simplex.set_A_coef(2,0,1);
    //  l_simplex.set_A_coef(3,1,1);

    //  l_simplex.define_equation_type(1,simplex::t_equation_type::INEQUATION_LT);
    //  l_simplex.define_equation_type(2,simplex::t_equation_type::INEQUATION_LT);
    //  l_simplex.define_equation_type(3,simplex::t_equation_type::INEQUATION_LT);

    std::ofstream l_output_file;
    l_output_file.open("test_case3.log");
    if(!l_output_file.is_open())
    {
        throw quicky_exception::quicky_runtime_exception("Unable to open file test_case3.log",__LINE__,__FILE__);
    }
    l_simplex.display_array(l_output_file);

    typename SIMPLEX_TYPE::t_coef_type l_max = 0;
    bool l_infinite = false;
    simplex::simplex_listener<typename SIMPLEX_TYPE::t_coef_type> l_listener(l_simplex,l_output_file);
    if(l_simplex.find_max(l_max,l_infinite,&l_listener))
    {
        std::cout << "Max = " << l_max << std::endl ;
    }
    else if(l_infinite)
    {
        std::cout << "Inifinite Max" << std::endl;
    }
    else
    {
        std::cout << "No Max found !?" << std::endl;
    }
    assert(3 == l_max);
    l_output_file.close();
}

#endif // SIMPLEX_SELF_TEST
//EOF
