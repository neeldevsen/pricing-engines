#include <iostream>
#include <cmath>
#include <random>
#include <algorithm>
#include <vector>
#include <numeric>
#include "../../Headers/jump_diffusion_struct.hpp"

template <typename S_t, typename K_t, typename r_t, typename sigma_t, typename T_t, typename lambda_t, typename mu_t, typename j_t>
auto jumpDiffusionExactTerminalSimCall(OptionDataJD<S_t, K_t, r_t, sigma_t, T_t, lambda_t, mu_t, j_t>& data, int M)
{
    using commonType = decltype(data.spot);

    commonType kappa_j {std::exp(data.mean_log_jump + static_cast<commonType>(0.5) * data.jump_size_volatility * data.jump_size_volatility) - static_cast<commonType>(1)};
    std::random_device rd{};
    std::mt19937 mt{rd()};
    std::normal_distribution W {static_cast<commonType>(0), static_cast<commonType>(1)};
    std::poisson_distribution P {data.jump_intensity * data.maturity};
    
    std::vector<commonType> X (M);
    std::vector<commonType> payoffs(M);
    commonType K_0 {std::log(data.spot) + (data.rate - data.jump_intensity * kappa_j - static_cast<commonType>(0.5) * data.volatility * data.volatility) * data.maturity};
    commonType K_1 {data.volatility * std::sqrt(data.maturity)};


    for (int i{0}; i < M; ++i)
    {
        commonType Zj {0};
        commonType deltaN {static_cast<commonType>(P(mt))};  
        if (deltaN == static_cast<commonType>(0))
        {
            Zj = 0;
        }
        else
        {
            std::normal_distribution J {deltaN * data.mean_log_jump, std::sqrt(deltaN) * data.jump_size_volatility };
            Zj = J(mt);
        }
        commonType Z {W(mt)};

        X[i] = K_0 + K_1 * Z + Zj ;
    }

    std::transform(X.begin(), X.end(), payoffs.begin(), [&data](commonType x){return std::max(std::exp(x) - data.strike, static_cast<commonType>(0));});
    return std::exp(-data.rate * data.maturity) * std::accumulate(payoffs.begin(), payoffs.end(), static_cast<commonType>(0)) /static_cast<commonType>(payoffs.size());
}

template <typename S_t, typename K_t, typename r_t, typename sigma_t, typename T_t, typename lambda_t, typename mu_t, typename j_t>
auto jumpDiffusionExactTerminalSimPut(OptionDataJD<S_t, K_t, r_t, sigma_t, T_t, lambda_t, mu_t, j_t>& data, int M)
{
    using commonType = decltype(data.spot);

    commonType kappa_j {std::exp(data.mean_log_jump + static_cast<commonType>(0.5) * data.jump_size_volatility * data.jump_size_volatility) - static_cast<commonType>(1)};
    std::random_device rd{};
    std::mt19937 mt{rd()};
    std::normal_distribution W {static_cast<commonType>(0), static_cast<commonType>(1)};
    std::poisson_distribution P {data.jump_intensity * data.maturity};
    
    std::vector<commonType> X (M);
    std::vector<commonType> payoffs(M);
    commonType K_0 {std::log(data.spot) + (data.rate - data.jump_intensity * kappa_j - static_cast<commonType>(0.5) * data.volatility * data.volatility) * data.maturity};
    commonType K_1 {data.volatility * std::sqrt(data.maturity)};


    for (int i{0}; i < M; ++i)
    {
        commonType Zj {0};
        commonType deltaN {static_cast<commonType>(P(mt))};  
        if (deltaN == static_cast<commonType>(0))
        {
            Zj = 0;
        }
        else
        {
            std::normal_distribution J {deltaN * data.mean_log_jump, std::sqrt(deltaN) * data.jump_size_volatility };
            Zj = J(mt);
        }
        commonType Z {W(mt)};

        X[i] = K_0 + K_1 * Z + Zj ;
    }

    std::transform(X.begin(), X.end(), payoffs.begin(), [&data](commonType x){return std::max(data.strike - std::exp(x), static_cast<commonType>(0));});
    return std::exp(-data.rate * data.maturity) * std::accumulate(payoffs.begin(), payoffs.end(), static_cast<commonType>(0)) /static_cast<commonType>(payoffs.size());
}

int main()
{
    OptionDataJD values {100, 100, 0.05 ,0.20, 1, 0.75, -0.10, 0.25};

    std::cout << jumpDiffusionExactTerminalSimCall(values, 10000000) << '\n';
    std::cout << jumpDiffusionExactTerminalSimPut(values, 10000000) << '\n';

}
