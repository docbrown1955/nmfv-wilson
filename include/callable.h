// This file is part of NMFV_WILSONS.
//
// NMFV_WILSONS is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// NMFV_WILSONS is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with NMFV_WILSONS. If not, see <https://www.gnu.org/licenses/>.

#ifndef CSL_LIB_CALLABLE
#define CSL_LIB_CALLABLE
#include <string>
#include <functional>

namespace c9_nmfv {

template<class ReturnType, class ParamType>
struct Callable {

    Callable(
            std::string_view t_name,
            ReturnType (*t_f)(ParamType const&)
            )
        :name(t_name),
         f(t_f)
    {}
    Callable(
            std::string_view t_name,
            std::function<ReturnType(ParamType const&)> const &t_f
            )
        :name(t_name),
         f(t_f)
    {}

    inline
    ReturnType operator()(ParamType const &params) const {
        return f(params);
    }

    std::string name;
    std::function<ReturnType(ParamType)> f;
};

} // End of namespace c9_nmfv

#endif
