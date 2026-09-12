# 7 Jump Diffusion

Jump Diffusion basically accounts for the price of the asset to randomly jump from one price to another, and not with the typical GBM behavior. The equation of an asset is as shown:
$$
\frac{dS_t}{S_t} = (r-\lambda \kappa_j)dt + \sigma dW_t + (J-1)dN_t
$$
Where $N_T \sim Poisson(\lambda T)$ ,  $\kappa_j = \exp(\mu + \frac{1}{2}\sigma^2_j$) , where $\lambda$ is the mean jump intensity, $J = \exp(\mu + \sigma_j  Z_j)$  where $Z_j \sim \mathcal{N}(0,1)$ ,  $\mu$ is the mean log jump of the asset and $\sigma_j$ is the jump size volatility of the asset.


