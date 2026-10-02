#ifndef BARRIERSTRUCT
#define BARRIERSTRUCT
#include <type_traits>

template <typename S_t, typename K_t, typename B_t, typename r_t, typename sigma_t, typename T_t>
struct OptionDataBarrier
{
    using commonType = std::common_type_t <double, S_t, K_t, B_t, r_t, sigma_t, T_t>;
    commonType spot {};
    commonType strike {};
    commonType barrier {};
    commonType rate {};
    commonType volatility {};
    commonType maturity {};
};

template <typename S_t, typename K_t, typename B_t, typename r_t, typename sigma_t, typename T_t>
OptionDataBarrier(S_t, K_t, B_t, r_t, sigma_t, T_t) -> OptionDataBarrier<S_t, K_t, B_t, r_t, sigma_t, T_t>;

template <typename T>
struct OptionValueBarrierUp
{
    using Q = std::common_type_t<double, T>;
    Q callUpOut {};
    Q callUpIn {};
    Q putUpOut {};
    Q putUpIn {};
    
};

template <typename T>
struct OptionValueBarrierDown
{
    using Q = std::common_type_t<double, T>;
    Q callDownOut {};
    Q callDownIn {};
    Q putDownOut {};
    Q putDownIn {}; 
};

#endif



