#ifndef FUNCTIONS
#define FUNCTIONS

#include <iostream>
#include <cmath>
#include <type_traits>

template <typename T>
auto normalCDF(T x) -> std::common_type_t <double , T>
{
    using commonType = std::common_type_t<double, T>;
    return static_cast<commonType>(0.5) * (static_cast<commonType>(1) + std::erf(x / std::sqrt(2.0)));
}


template <typename T>
auto max(T x, T y) -> std::common_type_t <double, T>
{
    if (x > y)
    {
        return x - y;
    }
    else
    {
        return 0;
    }
}


#endif