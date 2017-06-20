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
#include "quicky_exception.h"
#include "simplex.h"
#include <iostream>


void test_case1(void)
{
  // Example
  // Max Z -1 * X1 - 4 * X2 - 3 * X3           = 0
  //        2 * X1 + 2 * X2 + 1 * X3 + X4      = 4
  //            X1 + 2 * X2 + 2 * X3 +    + X5 = 6
  simplex::simplex<double> l_simplex(5, // Number of variables : x1 and x2
				  0, // Number of inequations with the form A x <= b
				  2, // Number of equations with the form A x = b
				  0  // Number of inequations with the form A x >= b
				  );
  l_simplex.set_Z_coef(0,1);
  l_simplex.set_Z_coef(1,4);
  l_simplex.set_Z_coef(2,3);
  l_simplex.set_B_coef(0,4);
  l_simplex.set_B_coef(1,6);
  l_simplex.set_A_coef(0,0,2);
  l_simplex.set_A_coef(0,1,2);
  l_simplex.set_A_coef(0,2,1);
  l_simplex.set_A_coef(0,3,1);
  l_simplex.set_A_coef(1,0,1);
  l_simplex.set_A_coef(1,1,2);
  l_simplex.set_A_coef(1,2,2);
  l_simplex.set_A_coef(1,4,1);
  l_simplex.define_equation_type(0,simplex::t_equation_type::EQUATION);
  l_simplex.define_equation_type(1,simplex::t_equation_type::EQUATION);

  double l_max = 0;
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
}

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

  simplex::simplex<double> l_simplex(2, // Number of variables : x1 and x2
				     4, // Number of inequations with the form A x <= b
				     0, // Number of equations with the form A x = b
				     0  // Number of inequations with the form A x >= b
				     );
  l_simplex.set_Z_coef(0,1000);
  l_simplex.set_Z_coef(1,1200);
  l_simplex.set_B_coef(0,200);
  l_simplex.set_B_coef(1,60);
  l_simplex.set_B_coef(2,34);
  l_simplex.set_B_coef(3,14);
  l_simplex.set_A_coef(0,0,10);
  l_simplex.set_A_coef(0,1,5);
  l_simplex.set_A_coef(1,0,2);
  l_simplex.set_A_coef(1,1,3);
  l_simplex.set_A_coef(2,0,1);
  l_simplex.set_A_coef(3,1,1);
  l_simplex.define_equation_type(0,simplex::t_equation_type::INEQUATION_LT);
  l_simplex.define_equation_type(1,simplex::t_equation_type::INEQUATION_LT);
  l_simplex.define_equation_type(2,simplex::t_equation_type::INEQUATION_LT);
  l_simplex.define_equation_type(3,simplex::t_equation_type::INEQUATION_LT);

  double l_max = 0;
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
}

//------------------------------------------------------------------------------
int main(int argc,char ** argv)
{
  try
    {
      std::cout << "============ TEST CASE 1 ==============" << std::endl;
      test_case1();
      std::cout << "============ TEST CASE 2 ==============" << std::endl;
      test_case2();
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
//EOF
