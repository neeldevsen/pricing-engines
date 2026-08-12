#include <iostream>
#include <cmath>
#include <random>
#include <algorithm>
#include <vector>
#include <numeric>
#include "../../Headers/jump_diffusion_struct.hpp"
#include "../../Headers/functions.hpp"

template <typename S_t, typename K_t, typename r_t, typename sigma_t, typename T_t, typename lambda_t, typename mu_t, typename j_t>
auto jumpDiffusionPoissonConditionedMCCall(OptionDataJD<S_t, K_t, r_t, sigma_t, T_t, lambda_t, mu_t, j_t>& data, int K_max, int M)
{
    using commonType = decltype(data.spot);

    commonType kappa_j {std::exp(data.mean_log_jump + static_cast<commonType>(0.5) * data.jump_size_volatility * data.jump_size_volatility) - static_cast<commonType>(1)};
    std::random_device rd{};
    std::mt19937 mt{rd()};
    std::normal_distribution W {static_cast<commonType>(0), static_cast<commonType>(1)};;
    

    std::vector<commonType> payoffs(M);
    commonType K_0 {std::log(data.spot) + (data.rate - data.jump_intensity * kappa_j - static_cast<commonType>(0.5) * data.volatility * data.volatility) * data.maturity};
    commonType K_1 {data.volatility * std::sqrt(data.maturity)};

    std::vector<commonType> probabilities (K_max + 1);
    
    commonType exp_factor {std::exp(-data.jump_intensity * data.maturity)}; 

    probabilities[0] = exp_factor;

    for (int k{1}; k <= K_max; ++k)
    {
        probabilities[k] = probabilities[k-1] * data.jump_intensity * data.maturity / k;
    }


    for (int i{0}; i < M; ++i)
    {
        std::vector<commonType> C_weighted (K_max + 1);
   
        for (int k{0}; k <= K_max; ++k)
        {
            commonType Z {W(mt)};
            commonType Zj {0};
            if (k != 0)
            {
                std::normal_distribution J {k * data.mean_log_jump, std::sqrt(k) * data.jump_size_volatility };
                Zj = J(mt);
            }
            else
            {
                Zj = 0;
            }
      
            C_weighted[k] = std::max(std::exp(K_0 + K_1 * Z + Zj) - data.strike , static_cast<commonType>(0)) * probabilities[k];

        }
        payoffs[i] = std::accumulate(C_weighted.begin(), C_weighted.end(), static_cast<commonType>(0));
    }

    return std::exp(-data.rate * data.maturity) * std::accumulate(payoffs.begin(), payoffs.end(), static_cast<commonType>(0)) /static_cast<commonType>(payoffs.size());
}

template <typename S_t, typename K_t, typename r_t, typename sigma_t, typename T_t, typename lambda_t, typename mu_t, typename j_t>
auto jumpDiffusionPoissonConditionedMCPut(OptionDataJD<S_t, K_t, r_t, sigma_t, T_t, lambda_t, mu_t, j_t>& data, int K_max, int M)
{
    using commonType = decltype(data.spot);

    commonType kappa_j {std::exp(data.mean_log_jump + static_cast<commonType>(0.5) * data.jump_size_volatility * data.jump_size_volatility) - static_cast<commonType>(1)};
    std::random_device rd{};
    std::mt19937 mt{rd()};
    std::normal_distribution W {static_cast<commonType>(0), static_cast<commonType>(1)};;
    

    std::vector<commonType> payoffs(M);
    commonType K_0 {std::log(data.spot) + (data.rate - data.jump_intensity * kappa_j - static_cast<commonType>(0.5) * data.volatility * data.volatility) * data.maturity};
    commonType K_1 {data.volatility * std::sqrt(data.maturity)};

    std::vector<commonType> probabilities (K_max + 1);
    
    commonType exp_factor {std::exp(-data.jump_intensity * data.maturity)}; 

    probabilities[0] = exp_factor;

    for (int k{1}; k <= K_max; ++k)
    {
        probabilities[k] = probabilities[k-1] * data.jump_intensity * data.maturity / k;
    }


    for (int i{0}; i < M; ++i)
    {
        std::vector<commonType> P_weighted (K_max + 1);
   
        for (int k{0}; k <= K_max; ++k)
        {
            commonType Z {W(mt)};
            commonType Zj {0};
            if (k != 0)
            {
                std::normal_distribution J {k * data.mean_log_jump, std::sqrt(k) * data.jump_size_volatility };
                Zj = J(mt);
            }
            else
            {
                Zj = 0;
            }
      
            P_weighted[k] = std::max(data.strike - std::exp(K_0 + K_1 * Z + Zj), static_cast<commonType>(0)) * probabilities[k];

        }
        payoffs[i] = std::accumulate(P_weighted.begin(), P_weighted.end(), static_cast<commonType>(0));
    }

    return std::exp(-data.rate * data.maturity) * std::accumulate(payoffs.begin(), payoffs.end(), static_cast<commonType>(0)) /static_cast<commonType>(payoffs.size());
}

int main()
{
    OptionDataJD values {100, 100, 0.05 ,0.20, 1, 0.75, -0.10, 0.25};

    std::cout << jumpDiffusionPoissonConditionedMCCall(values, 8, 1000000) << '\n';
    std::cout << jumpDiffusionPoissonConditionedMCPut(values, 8, 1000000) << '\n';

}