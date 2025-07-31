# Final Project for MA5770 - Bayesian Statistics

**Name:** Prateek Singh, Xinyun Liu, Manoj Kumar Surabhi, Xiaolin Liu  
**Project Title:** Bayesian Mixture Normal Model of Plasma Glucose Data  

## Section 1: Introduction  
The study analyzes plasma glucose levels from 532 females living near Phoenix, Arizona, who were tested for diabetes. The goal is to use a Bayesian approach with a mixture of normal distributions, each with different means and variances, to capture the data's skewness and outliers, where a single normal distribution may not adequately represent the data.

---

## Section 2: Statistical Analysis  
The distribution of glucose levels is observed by the histogram and kernel density function of the data. The measure of central tendencies (mean, median, and mode) and dispersion (variance and standard deviation) are calculated. A Bayesian mixture model approximates the data distribution using a mixture of two normal models. The Gibbs and JAGS samplers are used to obtain samples to approximate the posterior distribution of the data.  

**Figure 1** - Histogram and Kernel Density Estimation and $normal(\mu = 121.0,\ \sigma = 31.0)$  
![](media/image4.png)  

**Summary Statistics:**  
- Mean: 121.0  
- Median: 115.0  
- Mode: 100.0  
- Variance: 961.0  
- Standard Deviation: 31.0  

---

## Section 3: Statistical Method  
The Bayesian mixture normal model is defined as:  
$$
p\left( y_{i}|\pi,\Theta_{1},\Theta_{2},\sigma_{1}^{2},\sigma_{2}^{2} \right) = \pi \ast dnorm\left( y_{i};\theta_{1},\sigma_{1}^{2} \right) + (1 - \pi) \ast dnorm\left( y_{i};\theta_{2},\sigma_{2}^{2} \right)
$$

**Prior Distributions:**  
- $\pi \sim Beta(\alpha, \beta)$  
- $\theta_j \sim Normal(\mu_0, \tau_0^2)$ for $j = 1, 2$  
- $\sigma_j^2 \sim Inverse-Gamma\left(\frac{\upsilon_0}{2}, \frac{\upsilon_0 \sigma_0^2}{2}\right)$ for $j = 1, 2$  

---

## Section 4: Gibbs Sampler and MCMC Diagnostic  
**Initialization:**  
- $\alpha = 1, \beta = 1, \mu_0 = 120, \tau_0^2 = 200, \nu_0 = 10, \sigma_0^2 = 1000$  
- Latent variables $x_i$ initialized based on observed data.  

**Diagnostics:**  
- **Trace Plots:** Stable mixing observed.  
- **Autocorrelation:** Low values (0.1628 and 0.1419 for $\theta_{(1)}^{(s)}$).  
- **Effective Sample Size (ESS):** 134 and 94 for split chains.  

**Figure 2** - Histogram and Kernel Density Estimation for $\theta_{(1)}^{(s)}$  
![](media/image9.png)  

**Table 1** - Summary of Mean and 95% Confidence Interval for $\theta_{(1)}^{(s)}$  
| $\theta_{(1)}^{(s)}$ | Mean  | 95% CI (Lower) | 95% CI (Upper) |  
|----------------------|-------|----------------|----------------|  
| First Half           | 104.0 | 100.4          | 107.7          |  
| Last Half            | 103.9 | 100.4          | 108.0          |  

---

## Section 5: JAGS and MCMC Diagnostic  
**Model Specification in JAGS:**  
```jags
model {
  for (i in 1:n) {
    X[i] ~ dcat(pi_vals[])
    y[i] ~ dnorm(ifelse(X[i] == 1, theta1, theta2), ifelse(X[i] == 1, 1 / s21, 1 / s22))
  }
  pi ~ dbeta(alpha, beta)
  theta1 ~ dnorm(mu0, 1 / s21)
  theta2 ~ dnorm(mu0, 1 / s22)
  tau1 ~ dgamma(v0/2, v0 * s20 / 2)
  tau2 ~ dgamma(v0/2, v0 * s20 / 2)
  s21 <- 1 / tau1
  s22 <- 1 / tau2
}
