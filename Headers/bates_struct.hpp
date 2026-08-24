#ifndef BATES
#define BATES
#include <type_traits>

template <typename S_t, typename K_t, typename r_t, typename T_t, typename v_t, typename kappa_t, typename theta_t, typename xi_t, typename rho_t, typename lambda_t, typename mu_t, typename j_t>
struct OptionDataBates
{
    using commonType = std::common_type_t<double, S_t, K_t, r_t, T_t, v_t, kappa_t, theta_t, xi_t, rho_t , lambda_t, mu_t, j_t>;
    commonType spot {};
    commonType strike {};
    commonType rate {};
    commonType maturity {}; 
    commonType variance {};
    commonType kappa {};
    commonType theta{};
    commonType xi{};
    commonType rho{};
    commonType lambda {};
    commonType mu_j {};
    commonType sigma_j {};
};

template <typename S_t, typename K_t, typename r_t, typename T_t, typename v_t, typename kappa_t, typename theta_t, typename xi_t, typename rho_t, typename lambda_t, typename mu_t, typename j_t>
OptionDataBates(S_t, K_t, r_t, T_t, v_t, kappa_t, theta_t, xi_t, rho_t, lambda_t, mu_t, j_t) -> OptionDataBates<S_t, K_t, r_t, T_t, v_t, kappa_t, theta_t, xi_t, rho_t, lambda_t, mu_t, j_t>;

#endif












