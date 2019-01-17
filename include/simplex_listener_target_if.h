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

#ifndef SIMPLEX_SIMPLEX_LISTENER_TARGET_IF_H
#define SIMPLEX_SIMPLEX_LISTENER_TARGET_IF_H

#include <iostream>

namespace simplex
{
    /**
     * Interface class that should be implement by objects passed as parameters to simplex_listener
     */
    class simplex_listener_target_if
    {
      public:
        virtual std::ostream & display_array(std::ostream & p_stream) const = 0;
      private:
    };
}
#endif //SIMPLEX_SIMPLEX_LISTENER_TARGET_IF_H
