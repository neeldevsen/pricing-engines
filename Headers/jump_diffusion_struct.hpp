#ifndef JUMPDIFFUSION
#define JUMPDIFFUSION

#include <type_traits>

template <typename S_t, typename K_t, typename r_t, typename sigma_t, typename T_t, typename lambda_t, typename mu_t, typename j_t>
struct OptionDataJD 
{
    using commonType = std::common_type_t <double, S_t, K_t, r_t, sigma_t, T_t, lambda_t, mu_t, j_t>;
    commonType spot {};
    commonType strike {};
    commonType rate {};
    commonType volatility {};
    commonType maturity {};
    commonType jump_intensity {};
    commonType mean_log_jump {};
    commonType jump_size_volatility {};
    
};

template <typename S_t, typename K_t, typename r_t, typename sigma_t, typename T_t, typename lambda_t, typename mu_t, typename j_t>
OptionDataJD(S_t, K_t, r_t, sigma_t, T_t, lambda_t, mu_t, j_t) -> OptionDataJD<S_t, K_t, r_t, sigma_t, T_t, lambda_t, mu_t, j_t>;
#endif

