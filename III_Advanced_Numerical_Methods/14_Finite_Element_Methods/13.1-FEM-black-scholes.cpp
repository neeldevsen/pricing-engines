#include <iostream>
#include <cmath>
#include <vector>
#include <type_traits>
#include "../../Headers/black_scholes_struct.hpp"

template <typename S_t, typename K_t, typename r_t, typename sigma_t, typename T_t>
struct Coefficients
{
    using commonType = std::common_type_t<S_t, K_t, r_t, sigma_t, T_t>;
    std::vector <commonType> alphas{};
    std::vector <commonType> betas{};
    std::vector <commonType> gammas{};

    explicit Coefficients(int M)
    : alphas(M-1), betas(M-1), gammas(M-1)
    {}
};

template <typename S_t, typename K_t, typename r_t, typename sigma_t, typename T_t>
struct Diagonals
{
    using commonType = std::common_type_t<S_t, K_t, r_t, sigma_t, T_t>; 
    std::vector <commonType> A_lower {};
    std::vector <commonType> A_main {};
    std::vector <commonType> A_upper {};

    explicit Diagonals(int n)
    : A_lower(n-1), A_main(n), A_upper(n-1)
    {}

    Diagonals operator+(const Diagonals& other) const {
        int n {static_cast<int>(A_main.size())};
        Diagonals result(n);

        for (int i {0}; i < n-1; ++i){
            result.A_lower[i] = A_lower[i] + other.A_lower[i];
            result.A_upper[i] = A_upper[i] + other.A_upper[i];
            result.A_main[i] = A_main[i] + other.A_main[i];
        }
        result.A_main[n-1] = A_main[n-1] + other.A_main[n-1];

        return result;
    }

    Diagonals operator-(const Diagonals& other) const {
        int n {static_cast<int>(A_main.size()) 1};
        Diagonals result(n);

        for (int i {0}; i < n -1; ++i){
            result.A_lower[i] = A_lower[i] - other.A_lower[i];
            result.A_upper[i] = A_upper[i] - other.A_upper[i];
            result.A_main[i] = A_main[i] - other.A_main[i];
        }
        result.A_main[n-1] = A_main[n-1] - other.A_main[n-1];

        return result;
    }

    Diagonals operator*(const commonType factor) const {
        int n {static_cast<int>(A_main.size())};
        Diagonals result(n)
        for (int i {0}; i < n-1; ++i)
        {
            result.A_lower[i] = factor * A_lower[i];
            result.A_upper[i] = factor * A_upper[i];
            result.A_main[i] = factor * A_main[i];
        }
        result.A_main[n-1] = factor * A_main[n-1];
        return result;
    }

    Diagonals operator*(const std::vector<commonType> vec) const {
        int n {static_cast<int>(A_main.size())};
        std::vector<commonType> new_vec (n);
        for (int i {1}; i < n-1; ++i)
        {
            new_vec[i] = A_lower[i-1] * vec[i-1] + A_main[i] * vec[i] + A_upper[i] * vec[i+1];
        }
        new_vec[0] = A_main[0] * vec[0] + A_upper[0] * vec[1];
        new_vec[n-1] = A_lower[n-2] * vec[n-2] + A_main[n-1] * vec[n-1];

        return new_vec;                                                                                                                                                                         
    }
};

template <typename S_t, typename K_t, typename r_t, typename sigma_t, typename T_t>
Coefficients<S_t, K_t, r_t, sigma_t, T_t> generateCoefficientsBS(const OptionDataBS<S_t, K_t, r_t, sigma_t, T_t>& data, const double deltaS, const int M)
{
    Coefficients<S_t, K_t, r_t, sigma_t, T_t> coeff{M};
    using commonType = decltype(data.spot);
    for (int i{0}; i < M-1; ++i)
    {
        commonType S_i {(i+1)*deltaS};
        coeff.alphas[i] = 0.5 * (std::pow(data.volatility*S_i, 2.0)) * std::pow(static_cast<double>(i), 2.0) - 0.5 * data.rate * i;
        coeff.betas[i] = -(std::pow(data.volatility*S_i, 2.0)) * std::pow(static_cast<double>(i), 2.0) - data.rate;
        coeff.gammas[i] = 0.5 * (std::pow(data.volatility*S_i, 2.0)) * std::pow(static_cast<double>(i), 2.0) + 0.5 * data.rate * i;
    }
    return coeff;
}

template <typename S_t, typename K_t, typename r_t, typename sigma_t, typename T_t>
Diagonals<S_t, K_t, r_t, sigma_t, T_t> generateMassMatrix(const double deltaS, const int M)
{
    Diagonals<S_t, K_t, r_t, sigma_t, T_t> mass_diag {M};
    for (int i {0}; i < M-1; ++i)
    {
        mass_diag.A_main[i] = deltaS * 4.0 / 6.0;
        mass_diag.A_lower[i] = deltaS * 1.0 / 6.0;
        mass_diag.A_upper[i] =  mass_diag.A_lower[i] 
    }
    mass_diag.A_main[M-1] = deltaS * 2.0 / 6.0 ;
    mass_diag.A_main[0] -= deltaS * 2.0 / 6.0;
    return mass_diag;
}

template <typename S_t, typename K_t, typename r_t, typename sigma_t, typename T_t>
Diagonals<S_t, K_t, r_t, sigma_t, T_t> generateDiffusionMatrix(const double deltaS, const OptionDataBS<S_t, K_t, r_t, sigma_t, T_t>& data, const int M)
{
    Diagonals<S_t, K_t, r_t, sigma_t, T_t> diff_diag {M};
    for (int i {1}; i < M-1; ++i)
    {
        diff_diag.A_main[i] = 0.5 * deltaS * (data.volatility) * (data.volatility)  * (std::pow(((static_cast<double>(i)- 0.5), 2.0)) + std::pow(((static_cast<double>(i) + 0.5)), 2.0));
        diff_diag.A_lower[i] = - 0.5 * deltaS * (data.volatility) * (data.volatility)  * (std::pow(((static_cast<double>(i) - 0.5), 2.0)));
        diff_diag.A_upper[i] = diff_diag.A_lower[i]
    }
    diff_diag.A_main[M-1] =  0.5 * deltaS * (data.volatility) * (data.volatility)  * (std::pow(((static_cast<double>(M) - 0.5), 2.0)));
    diff_diag.A_main[0] = 0.125 * deltaS * (data.volatility) * (data.volatility); // would review this line because I did i = 0 directly in the calculation
    return mass_diag;
}

template <typename S_t, typename K_t, typename r_t, typename sigma_t, typename T_t>
Diagonals<S_t, K_t, r_t, sigma_t, T_t> generateConvectionMatrix(const double deltaS, const OptionDataBS<S_t, K_t, r_t, sigma_t, T_t>& data, const int M)
{
    Diagonals<S_t, K_t, r_t, sigma_t, T_t> conv_diag {M};
    for (int i {1}; i < M-1; ++i)
    { 
        conv_diag.A_main[i] = -0.5 * deltaS * (data.volatility * data.volatility - data.rate);
        conv_diag.A_upper[i] = 0.5 * deltaS * (data.volatility * data.volatility - data.rate)  * (static_cast<double>(i) + 0.5);
        conv_diag.A_lower[i] = -conv_diag.A_lower[i]
    }
    conv_diag.A_main[M-1] =  0.5 * deltaS * (data.volatility * data.volatility - data.rate)  * (static_cast<double>(M) - 0.5);
    conv_diag.A_main[0] = -0.25 * deltaS * (data.volatility * data.volatility - data.rate); // would review this line because I did i = 0 directly in the calculation
    conv_diag.A_upper[0] = 0.25 * deltaS * (data.volatility * data.volatility - data.rate);
    conv_diag.A_lower[0] = -conv_diag.A_lower[0]
    return conv_diag; // check this asap
}

template <typename S_t, typename K_t, typename r_t, typename sigma_t, typename T_t>
Diagonals<S_t, K_t, r_t, sigma_t, T_t> generateReactionMatrix(const double deltaS, const OptionDataBS<S_t, K_t, r_t, sigma_t, T_t>& data, const Diagonals<S_t, K_t, r_t, sigma_t, T_t>& mass_diag, const int M)
{
    Diagonals<S_t, K_t, r_t, sigma_t, T_t> reac_diag {M};
    for (int i {0}; i < M-1; ++i)
    { 
        reac_diag.A_main[i] = data.rate * mass_diag.A_main[i];
        reac_diag.A_upper[i] = data.rate * mass_diag.A_upper[i];
        reac_diag.A_lower[i] = data.rate * mass_diag.A_lower[i];
    }
    reac_diag.A_main[M-1] = data.rate * mass_diag.A_main[M-1];

    return reac_diag;
}

template <typename S_t, typename K_t, typename r_t, typename sigma_t, typename T_t>
Diagonals<S_t, K_t, r_t, sigma_t, T_t> generateCombinedMatrix(const Diagonals<S_t, K_t, r_t, sigma_t, T_t>& diff_diag, const Diagonals<S_t, K_t, r_t, sigma_t, T_t>& conv_diag, const Diagonals<S_t, K_t, r_t, sigma_t, T_t>& reac_diag, const int M)
{
    Diagonals<S_t, K_t, r_t, sigma_t, T_t> combined_diag {diff_diag + conv_diag + reac_diag};
    return combined_diag;
}

template <typename S_t, typename K_t, typename r_t, typename sigma_t, typename T_t>
Diagonals<S_t, K_t, r_t, sigma_t, T_t> generateLHSMatrix(const Diagonals<S_t, K_t, r_t, sigma_t, T_t>& M, const Diagonals<S_t, K_t, r_t, sigma_t, T_t>& A, const auto deltat)
{
    Diagonals<S_t, K_t, r_t, sigma_t, T_t> LHS {M + (A * (deltat / 2))};
    return LHS;
}

template <typename S_t, typename K_t, typename r_t, typename sigma_t, typename T_t>
auto generateRHSVectorCall(const OptionDataBS<S_t, K_t, r_t, sigma_t, T_t>& data, const Diagonals<S_t, K_t, r_t, sigma_t, T_t>& LHS, const Diagonals<S_t, K_t, r_t, sigma_t, T_t>& Mass, const Diagonals<S_t, K_t, r_t, sigma_t, T_t>& A, auto V, const auto deltat, const auto Smax, const int n, const int M) -> decltype(data.spot)
{
    using commonType = decltype(data.spot);
    std::vector<commonType> RHS {(Mass - (A * (0.5 * deltat))) * V};
    RHS[M-2] -= LHS.A_upper[M-2] * (Smax - data.strike * std::exp(-data.rate * deltat * (n + 1))); 
    RHS[M-2] += 
    return RHS;
}

template <typename S_t, typename K_t, typename r_t, typename sigma_t, typename T_t>
auto generateRHSPut(const OptionDataBS<S_t, K_t, r_t, sigma_t, T_t>& data, Coefficients<S_t, K_t, r_t, sigma_t, T_t>& coeff,  auto V, const int n, const int M, const auto deltaT, const auto Smax, const int M)
{
    using commonType = decltype(data.spot);
    V[0] = data.strike * std::exp(-data.rate*static_cast<commonType>(n)*deltaT);
    V[M] = 0;
    std::vector<commonType> RHS (M-1);
    for (int i{0}; i < M-1; ++i)
    {
        RHS[i] = static_cast<commonType>(0.5)*deltaT*coeff.alphas[i] * V[i] + (static_cast<commonType>(1) + static_cast<commonType>(0.5)*deltaT*coeff.betas[i])*V[i+1] + static_cast<commonType>(0.5)*deltaT * coeff.gammas[i]*V[i+2];
    }
    RHS[0] += static_cast<commonType>(0.5) * deltaT * coeff.alphas[0] * (data.strike * std::exp(-data.rate * (static_cast<commonType>(n) + static_cast<commonType>(1)) * deltaT ));
    return RHS; 
}

template <typename S_t, typename K_t, typename r_t, typename sigma_t, typename T_t>
auto generateLHSMatrix(const OptionDataBS<S_t, K_t, r_t, sigma_t, T_t>& data, Coefficients<S_t, K_t, r_t, sigma_t, T_t>& coeff,  auto V, const int n, const int M, const auto deltaT, const auto Smax)
{
    using commonType = decltype(data.spot);
    std::vector<commonType> RHS (M-1);
    for (int i{0}; i < M-1; ++i)
    {
        RHS[i] = static_cast<commonType>(0.5)*deltaT*coeff.alphas[i] * V[i] + (static_cast<commonType>(1) + static_cast<commonType>(0.5)*deltaT*coeff.betas[i])*V[i+1] + static_cast<commonType>(0.5)*deltaT * coeff.gammas[i]*V[i+2];
    }
    RHS[0] += static_cast<commonType>(0.5) * deltaT * coeff.alphas[0] * (data.strike * std::exp(-data.rate * (static_cast<commonType>(n) + static_cast<commonType>(1)) * deltaT ));
    return RHS;  
}

template <typename S_t, typename K_t, typename r_t, typename sigma_t, typename T_t>
auto thomasAlgorithm(Diagonals<S_t, K_t, r_t, sigma_t, T_t>& diag, const double& RHS, const int M)
{
    using commonType = std::remove_reference_t<decltype(diag.A_upper[0])>;
    std::vector<commonType> rho_p (M-2);
    std::vector<commonType> mumass_diag_p (M-1);

    rho_p[0] = diag.A_upper[0] / diag.A_main[0];
    mu_p[0] = RHS[0] / diag.A_main[0];

    for (int i {1}; i < M - 2; ++i)
    {
        mu_p[i] = (RHS[i] / diag.A_main[i] - diag.A_lower[i-1] / diag.A_main[i] * mu_p[i-1]) / (static_cast<commonType>(1) - diag.A_lower[i-1] / diag.A_main[i] * rho_p[i-1]);
        rho_p[i] = diag.A_upper[i] / diag.A_main[i] / (1.0 - diag.A_lower[i-1] / diag.A_main[i] * rho_p[i-1]);
    }
    mu_p[M-2] = (RHS[M-2] / diag.A_main[M-2] - diag.A_lower[M-3] / diag.A_main[M-2] * mu_p[M-3]) / (static_cast<commonType>(1) - diag.A_lower[M-3] / diag.A_main[M-2] * rho_p[M-3]);

    std::vector<commonType> Vn (M);
    Vn[M-2] = mu_p[M-2];
    for (int i{M-2}; i >= 1; --i)
    {
        Vn[i] = mu_p[i-1] - rho_p[i-1] * Vn[i+1];
    }
    return Vn;
}

template <typename S_t, typename K_t, typename r_t, typename sigma_t, typename T_t>
auto stepCrankNicolsonCall(const OptionDataBS<S_t, K_t, r_t, sigma_t, T_t>& data, Coefficients<S_t, K_t, r_t, sigma_t, T_t>& coeff, Diagonals<S_t, K_t, r_t, sigma_t, T_t>& diag, const auto deltaT, const int M, const int n, const auto Smax, auto V)
{
    using commonType = decltype(data.spot);
    std::vector<commonType> RHS{generateRHSCall(data, coeff, V, n, M, deltaT, Smax)};
    std::vector<commonType> V_next {thomasAlgorithm(diag, RHS, M)};
    V_next[0] = 0;
    V_next.push_back(Smax - data.strike * std::exp(-data.rate*deltaT*(static_cast<commonType>(n+1))));
    return V_next;  
}

template <typename S_t, typename K_t, typename r_t, typename sigma_t, typename T_t>
auto stepCrankNicolsonPut(const OptionDataBS<S_t, K_t, r_t, sigma_t, T_t>& data, Coefficients<S_t, K_t, r_t, sigma_t, T_t>& coeff, Diagonals<S_t, K_t, r_t, sigma_t, T_t>& diag, const auto deltaT, const int M, const int n, const auto Smax, auto V)
{
    using commonType = decltype(data.spot);
    std::vector<commonType> RHS{generateRHSPut(data, coeff, V, n, M, deltaT, Smax)};
    std::vector<commonType> V_next {thomasAlgorithm(diag, RHS, M)};
    V_next[0] = data.strike * std::exp(-data.rate*deltaT*(static_cast<commonType>(n+1))); 
    V_next.push_back(0);
    return V_next;  
}

template <typename S_t, typename K_t, typename r_t, typename sigma_t, typename T_t>
auto finiteDifferencesCallBS(const OptionDataBS<S_t, K_t, r_t, sigma_t, T_t>& data, const int M, const int N, const double multiplier=4.0)
{
    using commonType = decltype(data.spot);
    commonType Smax {multiplier * std::max(data.spot,data.strike)};
    commonType deltaS {Smax / static_cast<commonType>(M)};
    commonType deltaT {data.maturity / static_cast<commonType>(N)};
    std::vector<commonType> S (M + 1);
    for (int i{}; i < M + 1; ++i)
    {
        S[i] = static_cast<commonType>(i) * deltaS;
    }

    std::vector<commonType> V (M + 1);
    for (int i{}; i < M + 1; ++i)
    {
        V[i] = std::max(S[i] - data.strike, static_cast<double>(0));
    }

    Coefficients coeff {generateCoefficientsBS(data, deltaS, M)};
    Diagonals diag {generateDiagonalsBS(coeff, deltaT, M)};

    for (int n{}; n < N; ++n)
    {
        auto V_new {stepCrankNicolsonCall(data, coeff, diag , deltaT, M, n, Smax, V)};   
        std::swap(V_new, V);
    }
    int j {static_cast<int>(data.spot / deltaS)};
    j = std::max(static_cast<int>(0), static_cast<int>(std::min(M-1,static_cast<int>(j))));
    double price {V[j] + (data.spot - S[j]) * (V[j+1] - V[j]) / (S[j+1] - S[j])};
    return price;
}

template <typename S_t, typename K_t, typename r_t, typename sigma_t, typename T_t>
auto finiteDifferencesPutBS(const OptionDataBS<S_t, K_t, r_t, sigma_t, T_t>& data, const int M, const int N, const double multiplier=4.0)
{
    using commonType = decltype(data.spot);
    commonType Smax {multiplier * std::max(data.spot,data.strike)};
    commonType deltaS {Smax / static_cast<commonType>(M)};
    commonType deltaT {data.maturity / static_cast<commonType>(N)};
    std::vector<commonType> S (M + 1);
    for (int i{}; i < M + 1; ++i)
    {
        S[i] = static_cast<commonType>(i) * deltaS;
    }

    std::vector<commonType> V (M + 1);
    for (int i{}; i < M + 1; ++i)
    {
        V[i] = std::max(data.strike - S[i], static_cast<double>(0));
    }

    Coefficients coeff {generateCoefficientsBS(data, deltaS, M)};
    Diagonals diag {generateDiagonalsBS(coeff, deltaT, M)};

    for (int n{}; n < N; ++n)
    {
        auto V_new {stepCrankNicolsonPut(data, coeff, diag , deltaT, M, n, Smax, V)};   
        std::swap(V_new, V);
    }
    int j {static_cast<int>(data.spot / deltaS)};
    j = std::max(static_cast<int>(0), static_cast<int>(std::min(M-1,static_cast<int>(j))));
    double price {V[j] + (data.spot - S[j]) * (V[j+1] - V[j]) / (S[j+1] - S[j])};
    return price;
}

int main()
{
    OptionDataBS myOption {100, 100, 0.01, 0.3, 30};
    double multiplier {40};
    int M {5000};
    int N {5000};
    std::cout << "Call price: $" << finiteDifferencesCallBS(myOption, M, N, multiplier) << '\n';
    std::cout << "Put price: $" << finiteDifferencesPutBS(myOption, M, N, multiplier) << '\n';

    return 0;
}

