#ifndef HESTON
#define HESTON
#include <type_traits>

template <typename S_t, typename K_t, typename r_t, typename T_t, typename v_t, typename kappa_t, typename theta_t, typename xi_t, typename rho_t>
struct OptionDataHeston
{
    using commonType = std::common_type_t<double, S_t, K_t, r_t, T_t, v_t, kappa_t, theta_t, xi_t, rho_t>;
    commonType spot {};
    commonType strike {};
    commonType rate {};
    commonType maturity {}; 
    commonType variance {};
    commonType kappa {};
    commonType theta{};
    commonType xi{};
    commonType rho{};
};

template <typename S_t, typename K_t, typename r_t, typename T_t, typename v_t, typename kappa_t, typename theta_t, typename xi_t, typename rho_t>
OptionDataHeston(S_t, K_t, r_t, T_t, v_t, kappa_t, theta_t, xi_t, rho_t) -> OptionDataHeston<S_t, K_t, r_t, T_t, v_t, kappa_t, theta_t, xi_t, rho_t>;

#endif