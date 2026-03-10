#set page(
  paper: "a4",
  flipped: true,
  margin: (x: 1.5cm, y: 1.5cm),
  columns: 3,
)

#set text(font: "Libertinus Serif", size: 9pt)
#set heading(numbering: "1.", outlined: false)

#show heading.where(level: 1): it => block(
  text(fill: rgb("#003366"), weight: "bold", size: 12pt, it),
  below: 0.8em,
  above: 1.2em,
)
#show heading.where(level: 2): it => block(
  text(fill: rgb("#0055A4"), weight: "bold", size: 10pt, it),
  below: 0.6em,
  above: 1em,
)

#align(center)[
  #text(size: 16pt, weight: "bold", fill: rgb("#003366"))[MIDA 1: Complete Cheatsheet]\
  #text(size: 10pt)[Model Identification and Data Analysis -- Prof. Simone Formentin]\
  #v(1em)
]

= 1. Stochastic Processes
A stochastic process is a sequence of random variables indexed by time.

*Weak Description:*
- *Expected Value:* $m_v = E[v(t)]$ describes the average behavior.
- *Covariance Function:* Tells us how data spread around the expected value and correlate with each other over a time lag $tau$.
  $ gamma_v(tau, t) = E[(v(t) - E[v(t)])(v(t-tau) - E[v(t-tau)])] $

*Stationarity:*
A stochastic process is stationary if its statistical properties do not depend on time:
$ E[v(t)] = m (forall t), quad gamma_v(tau, t) = gamma_v(tau) $

*Properties of Covariance Function:*
For a stationary process:
+ $gamma_y(0) >= 0$ (This is the variance)
+ $|gamma_y(tau)| <= gamma_y(0)$ for $tau != 0$
+ $gamma_y(-tau) = gamma_y(tau)$ (Even function)

*White Noise (WN):*
$e(t) ~ text("WN")(mu, lambda^2)$. An uncorrelated process.
- Mean: $mu = E[e(t)]$
- Variance: $gamma_e(0) = lambda^2$
- Covariance: $gamma_e(tau) = 0$ for all $tau != 0$.

*Debias/Detrend Process:*
To compute covariance easily when $E[y(t)] != 0$, define the debiased process:
$ tilde(y)(t) = y(t) - E[y(t)] $
Then $gamma_y(tau) = gamma_(tilde(y))(tau)$.

= 2. Time Series Models (ARMA)
*Moving Average -- MA(q):*
$ y(t) = c_0 e(t) + c_1 e(t-1) + dots.c + c_q e(t-q) $
All MA processes are strictly stationary.

*Autoregressive -- AR(p):*
$ y(t) = a_1 y(t-1) + a_2 y(t-2) + dots.c + a_p y(t-p) + e(t) $
Stationary if all poles are strictly inside the unit circle ($|z| < 1$).

*Theorem of the Gain:*
At steady state, for $y(t) = G(z)u(t)$:
$ E[y(t)] = G(1) dot E[u(t)] $

= 3. Frequency Domain Analysis
The spectrum maps the signal's properties into the frequency domain.
$ Gamma_y(omega) = sum_(tau = -oo)^(+oo) gamma_y(tau) e^(-j omega tau) $

*Main Properties of the Spectrum:*
1. $Gamma_y(omega)$ is real.
2. $Gamma_y(omega) >= 0$.
3. $Gamma_y(omega) = Gamma_y(-omega)$ (Even function).
4. Periodic with period $T = 2pi$.
5. For a white noise process: $Gamma_e(omega) = lambda^2 quad forall omega$.

*Fundamental Theorem of Spectral Analysis:*
If $y(t) = W(z)u(t)$ and $u(t)$ is stationary:
$ Gamma_y(omega) = |W(e^(j omega))|^2 Gamma_u(omega) $

*Inverse Formula (Variance computation):*
$ gamma_y(0) = 1 / (2pi) integral_(-pi)^pi Gamma_y(omega) dif omega $

*Linearity:*
- Scalar multiple: if $y = a x$, then $Gamma_y = a^2 Gamma_x$.
- Uncorrelated sum: if $z = x + y$, $x,y$ uncorr., then $Gamma_z = Gamma_x + Gamma_y$.

*Theorem of the Gain:*
For a stable $W(z)$: $E[y(t)] = W(1) dot E[u(t)]$. Useful for de-biasing.

= 4. Optimal Prediction (Kolmogorov-Wiener)
Goal: Predict future values based on past observations.

*Canonical Representation:*
To design the optimal predictor, the model $y(t) = (C(z))/(A(z)) e(t)$ must be in canonical form:
1. Numerator $C(z)$ and Denominator $A(z)$ are *monic* (coefficient of highest power of $z$ is 1).
2. *Coprime* (no common factors; simplify if any).
3. Have the *same degree*.
4. *Poles and Zeroes* must be strictly inside the unit disk ($|z| < 1$).

*One-Step Predictor ($k=1$):*
$ hat(y)(t|t-1) = (C(z) - A(z)) / C(z) y(t) $
- Prediction Error: $epsilon(t|t-1) = y(t) - hat(y)(t|t-1) = e(t)$
- The variance of the 1-step prediction error is minimal and equals $lambda^2$.

*1-Step ARMA Shortcuts:*
- From noise: $hat(y)(t|t-1) = (C(z) - A(z))/A(z) e(t)$
- From data: $hat(y)(t|t-1) = (C(z) - A(z))/C(z) y(t)$
- Prediction error (whitening filter): $epsilon(t|t-1) = A(z)/C(z) y(t)$

*Multi-step AR(1):* $hat(y)(t+r|t) = a^r y(t)$ — converges to $E[y]$ as $r -> oo$.

*k-Step Predictor ($k >= 2$):*
Use polynomial long division to expand the transfer function for $k$ steps:
$ C(z) / A(z) = E(z) + z^(-k) F(z) / A(z) $
Then, the optimal $k$-step predictor from data is:
$ hat(y)(t|t-k) = F(z) / C(z) y(t) $
- Prediction Error Variance: $E[epsilon(t|t-k)^2] = lambda^2 (1 + e_1^2 + e_2^2 + dots.c + e_(k-1)^2)$.
- As $k -> oo$, the prediction tends to the expected value $E[y(t)]$.

= 5. System Identification
Find parameters $theta$ and variance $lambda^2$ that best fit the data.

*Asymptotic Analysis ($oo$ samples):*
Minimize the theoretical variance of the 1-step prediction error:
$ overline(J)(theta) = E[(y(t) - hat(y)(t|t-1, theta))^2] $
$ theta^* = op("argmin")_theta overline(J)(theta), quad lambda_*^2 = overline(J)(theta^*) $

*Quasi-Newton Update:*
Approximate Hessian (always PSD):
$ H approx 2/N sum (partial epsilon)/(partial theta) ((partial epsilon)/(partial theta))^T $
$ theta^((i+1)) = theta^((i)) - H^(-1) [sum epsilon dot (partial epsilon) / (partial theta)] $

*Auxiliary Signals* (for computing $(partial epsilon)/(partial theta)$):
$ alpha(t) = -1/C(z) y(t), quad beta(t) = -1/C(z) u(t), quad gamma(t) = -1/C(z) epsilon(t) $

*Four Cases of PEM Convergence:*
- $cal(S) in cal(M)$, $Delta$ singleton: $hat(theta)_N -> theta^0$ (unique correct solution).
- $cal(S) in cal(M)$, $Delta$ not singleton: converges to one in $Delta$ (over-parameterization).
- $cal(S) in.not cal(M)$, $Delta$ singleton: $hat(theta)_N -> theta^*$ (best proxy).
- $cal(S) in.not cal(M)$, $Delta$ not singleton: no guarantees on which minimum.

*Finite Samples ($N$ samples):*
$ J_N(theta) = 1 / N sum_(t=1)^N (y(t) - hat(y)(t|t-1, theta))^2 $

If the prediction error $epsilon(t|t-1)$ evaluated at $theta^*$ is white noise, the system model matches the true generating process ($cal(S) in cal(M)$). If not, verify if the coefficients of the residual correlation are reasonably small.

= 6. The Box-Jenkins Method
A checklist to go from raw data to a production model.

== Step 1: Identification
- *Check Stationarity:* Plot data (`df.plot()`). Use the Augmented Dickey-Fuller (ADF) test.
  - $H_0$: Time series is non-stationary.
  - If $p$-value $< 0.05$, reject $H_0$ (Data is stationary).
- *Transformations:* If non-stationary, apply Differencing ($Delta y(t) = y(t) - y(t-1)$) or logarithms to stabilize variance. (ARIMA becomes ARMA).
- *Pick Orders (p, q):*
  - *ACF (Autocorrelation Function):* Cuts off after lag $q$ for an MA($q$) process.
  - *PACF (Partial Autocorrelation Function):* Cuts off after lag $p$ for an AR($p$) process.

== Step 2: Estimation
- Train the model (`model.fit()`).
- Compare candidate models using objective criteria:
  - *AIC (Akaike Information Criterion):* Better for predictive models.
  - *BIC (Bayesian Information Criterion):* Heavily penalizes complexity; chooses simpler, explanatory models. Lower is better for both.

== Step 3: Model Diagnostics
Check if the model captured all systematic information (residuals should be White Gaussian Noise).
- *MAE (Mean Absolute Error):* Useful accuracy statistic.
- *Ljung-Box Test:* Checks residual autocorrelation.
  - $H_0$: Residuals are uncorrelated.
  - *Prob(Q) $> 0.05$:* Good! Fail to reject $H_0$. Residuals are white noise.
- *Jarque-Bera Test:* Checks normality.
  - *Prob(JB) $> 0.05$:* Residuals look normally distributed.
- *Plots:* Standardized residuals, Q-Q plot (should be a straight line), and Correlogram (ACF of residuals should stay within confidence bounds).

*Train/Validation Split:*
Partition data into identification (training) and validation sets. Avoids overfitting but wastes data.

== Step 4: Production
Once the diagnostics confirm the model is adequate, use it to forecast (`results.get_forecast()`). Forecasts for larger horizons will smoothly decay toward the historical mean.
