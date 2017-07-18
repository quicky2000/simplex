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
#ifndef _SIMPLEX_ARRAY_BASE_H_
#define _SIMPLEX_ARRAY_BASE_H_

namespace simplex
{
  class simplex_array_base
  {
  public:
    inline simplex_array_base(const unsigned int & p_nb_equations,
			      const unsigned int & p_nb_variables
			      );
    inline const unsigned int & get_nb_equations(void)const;
    inline const unsigned int & get_nb_variables(void)const;
  private:
    /**
       Number of lines
    */
    unsigned int m_nb_equations;

    /**
       Number of columns
    */
    unsigned int m_nb_variables;
  };

  //----------------------------------------------------------------------------
  simplex_array_base::simplex_array_base(const unsigned int & p_nb_equations,
					 const unsigned int & p_nb_variables
					 ):
    m_nb_equations(p_nb_equations),
    m_nb_variables(p_nb_variables)
    {
    }

  //----------------------------------------------------------------------------
  const unsigned int & simplex_array_base::get_nb_equations(void)const
    {
      return m_nb_equations;
    }

  //----------------------------------------------------------------------------
  const unsigned int & simplex_array_base::get_nb_variables(void)const
    {
      return m_nb_variables;
    }
}
#endif // _SIMPLEX_ARRAY_BASE_H_
// EOF
