// Compact cheatsheet — all key formulas, procedures, and decision rules
// Designed for quick reference during exercises and exams

#import "../prelude.typ": *

#set page(flipped: true, paper: "a4", margin: (x: 1.5cm, y: 1.5cm))
#set text(font: "Libertinus Serif", size: 9pt)
#set par(leading: 0.5em, spacing: 0.65em)
#set block(spacing: 0.65em)
#set heading(numbering: none)
#set enum(indent: 0pt, spacing: 0.45em)
#set list(indent: 0pt, spacing: 0.45em)

#show heading.where(level: 1): it => block(
  text(fill: rgb("#003366"), weight: "bold", size: 14pt, it),
  below: 0.6em,
  above: 0.8em,
)
#show heading.where(level: 2): it => block(
  text(fill: rgb("#003366"), weight: "bold", size: 11pt, it),
  below: 0.5em,
  above: 0.8em,
)
#show heading.where(level: 3): it => block(
  text(fill: rgb("#0055A4"), weight: "bold", size: 9.5pt, it),
  below: 0.4em,
  above: 0.6em,
)
#show heading.where(level: 4): set text(size: 9pt)

#let hline() = line(length: 100%, stroke: 0.4pt + luma(160))
#let bb(body) = text(weight: "bold", body)

// Global table styling
#set table(
  stroke: (x: none, y: 0.4pt + luma(200)),
  inset: (x: 6pt, y: 4pt),
  fill: (_, row) => {
    if row == 0 { rgb("#d8e2ee") } else if calc.odd(row) { rgb("#f0f3f8") }
  },
)

#align(center)[
  #text(size: 16pt, weight: "bold", fill: rgb("#003366"))[MIDA 1: Complete Cheatsheet]\
  #text(size: 10pt)[Model Identification and Data Analysis]\
  #v(0.5em)
]

#columns(2, gutter: 30pt)[

  == Models

  #hline()

  === At a glance

  #table(
    columns: (auto, auto, 1fr, auto, auto),
    table.header([*Model*], [*Definition*], [*$W(z)$*], [*Poles/Zeros*], [*Memory*]),
    [WN], [$e(t)$], [$1$], [—], [None],
    [MA($n$)], [$sum_(j=0)^n c_j e(t-j)$], [$C(z)$], [$n$ zeros, all-zeros], [Finite ($n$)],
    [AR($m$)], [$sum_(i=1)^m a_i y(t-i) + e(t)$], [$1\/A(z)$], [$m$ poles, all-poles], [Infinite],
    [ARMA], [$dots.c$], [$C(z)\/A(z)$], [$m$ poles + $n$ zeros], [Infinite],
    [ARMAX], [$dots.c + B(z)u(t)$], [$C(z)\/A(z)$, $B(z)\/A(z)$], [$m$ poles + $n$,$p$ zeros], [Infinite],
  )

  === Properties of Covariance Function
  #table(
    columns: (auto, 1fr),
    [*Non-negative variance*], [$gamma(0) >= 0$],
    [*Maximum at lag 0*], [$|gamma(tau)| <= gamma(0)$ for $tau eq.not 0$],
    [*Even function*], [$gamma(-tau) = gamma(tau)$],
  )

  === WN

  #table(
    columns: (auto, 1fr),
    [*Mean*], [$m_e = mu$],
    [*Variance*], [$gamma(0) = lambda^2$],
    [*Covariance*], [$gamma(tau) = 0 space forall tau eq.not 0$],
    [*Spectrum*], [$Gamma(omega) = lambda^2$],
  )

  === MA($n$): $y(t) = c_0 e(t) + c_1 e(t-1) + dots.c + c_n e(t-n)$

  #table(
    columns: (auto, 1fr),
    [*Mean*], [$m_y = mu sum_(i=0)^n c_i$],
    [*Monic*], [Set $c_0 = 1$ → $n+1$ free params],
    [*Covariance*],
    [$gamma(tau) = lambda^2 sum_(i=0)^(n-|tau|) c_i c_(i+|tau|)$, $space gamma(tau) = 0$ for $|tau| > n$],

    [*MA(1) spectrum*], [$Gamma(omega) = (1+c^2) lambda^2 + 2c lambda^2 cos(omega)$],
  )

  === AR($m$): $y(t) = a_1 y(t-1) + dots.c + a_m y(t-m) + e(t)$

  #table(
    columns: (auto, 1fr),
    [*Stability*], [All roots of $A(z)$ strictly inside unit disk ($|z_i| < 1$)],
    [*Key property*], [$EE[y(t-tau) e(t)] = 0 space forall tau > 0$],
    [*AR→MA∞*], [$1/(1-a z^(-1)) = sum_(k=0)^oo a^k z^(-k)$ (iff $|a|<1$)],
    [*AR(1) cov.*], [$gamma(0) = lambda^2 / (1-a^2), quad gamma(tau) = a dot gamma(tau-1)$ for $|tau| >= 1$],
    [*AR(1) spectrum*], [$Gamma(omega) = lambda^2 / (1 + a^2 - 2a cos(omega))$],
    [*Yule-Walker*], [$gamma(tau) = sum_(i=1)^m a_i gamma(tau-i)$ for $|tau| >= 1$],
    [*Matrix form*], [$va(a) = Gamma^(-1) va(gamma), quad lambda^2 = gamma(0) - sum_(i=1)^m a_i gamma(i)$],
  )

  === ARMA($m,n$): $y(t) = sum_(i=1)^m a_i y(t-i) + sum_(j=0)^n c_j e(t-j)$

  #table(
    columns: (auto, 1fr),
    [*$W(z)$*], [$C(z)/A(z)$],
    [*ARMA(1,1) $gamma(0)$*], [$((1+c_1^2)(1-a^2) + 2a c_1)/(1-a^2) lambda^2$],
    [*ARMA(1,1) $gamma(1)$*], [$a gamma(0) + c_1 lambda^2$],
    [*Tail ($|tau| >= 2$)*], [$gamma(tau) = a dot gamma(tau-1)$ (decays like AR)],
  )

  === ARMAX($m,n,p,k_0$) = $sum_(i=1)^m a_i y(t-i) + sum_(j=0)^n c_j e(t-j) + sum_(l=0)^p b_l u(t-k_0-l)$

  #table(
    columns: (auto, 1fr),
    [*Definition*], [$y(t) = B(z)/A(z) z^(-k_0) u(t) + C(z)/A(z) e(t)$],
    [*Mean (gain)*], [$m_y = W(1) dot mu = C(1)/A(1) dot mu$],
  )

  == Frequency Domain

  #hline()

  #bb[Spectral density:] (real, $>=0$, even, $2pi$-periodic)
  $ Gamma(omega) = sum_(tau=-oo)^(+oo) gamma(tau) e^(-j omega tau) $

  #bb[Inverse:]
  $ gamma(tau) = 1/(2pi) integral_(-pi)^pi Gamma(omega) e^(j omega tau) dif omega $
  #bb[Variance:]
  $ gamma(0) = 1/(2pi) integral_(-pi)^pi Gamma(omega) dif omega $

  #bb[Spectral factorization:] if $y(t) = W(z) v(t)$, $W(z)$ stable:
  $ Gamma_y (omega) = |W(e^(j omega))|^2 Gamma_v (omega) $

  #bb[Complex spectrum:]
  $ Phi(z) = sum gamma(tau) z^(-tau), quad Gamma(omega) = Phi(e^(j omega)) $

  #bb[Cross-spectrum:]
  $ Phi_(u y)(z) = W(z) Phi_(u u)(z), quad Phi_(y y) = W(z) W(z^(-1)) Phi_(u u) $

  #bb[Graphical method:] $overline(a), overline(b)$ = distances from zero/pole to $e^(j omega)$:
  $ Gamma_y (omega) = |overline(a) / overline(b)|^2 lambda^2 $

  === Spectrum Properties & Qualitative Analysis
  #table(
    columns: (auto, 1fr),
    [*Main properties*], [Real, $Gamma(omega) >= 0$, even ($Gamma(omega) = Gamma(-omega)$), $2pi$-periodic],
    [*Poles near unit disk*], [Spectrum is high (approaches infinity on the disk)],
    [*Zeros near unit disk*], [Spectrum is low (is zero on the disk)],
  )

  #colbreak()
  == Canonical Representation

  #hline()

  $y(t) = hat(W)(z) e(t) = C(z)/A(z) e(t)$ is canonical iff:

  #enum(
    [#bb[Monic:] $C(z), A(z)$ have leading coeff $= 1$],
    [#bb[Null relative degree:] $deg C = deg A$ (no time shift)],
    [#bb[Coprime:] $C(z), A(z)$ share no common factors],
    [#bb[Min.\ phase:] roots of $C(z)$ inside $|z| <= 1$; roots of $A(z)$ inside $|z| < 1$],
  )

  #bb[All-pass filter:] $|q|>1$, preserves spectrum ($T(z)T(z^(-1)) = 1$):
  $ T(z) = 1/q dot (z-q)/(z - 1/q) $


  == Prediction (Kolmogorov-Wiener)

  === Debias Process (Non-zero expected value)
  For processes with a non-null expected value, define the debias processes:
  $ tilde(y)(t) = y(t) - EE[y(t)], quad tilde(e)(t) = e(t) - EE[e(t)] $
  Compute the predictor using the debias process, then add the expected value back.

  #hline()

  Given $y(t) = C(z)/A(z) e(t)$ in #bb[canonical form]:

  === Step 1 — Diophantine equation (long division)

  Divide $C(z)$ by $A(z)$ for $k$ steps:
  $ C(z) = E(z) A(z) + z^(-k) F(z) $

  - $E(z)$: quotient, $deg = k-1$ #h(1fr) - $F(z)$: remainder, $deg = max(m, n)-1$

  === Step 2 — Predictor from noise
  $ hat(y)(t+k|t) = F(z)/A(z) e(t) $

  === Step 3 — Predictor from data
  Apply whitening filter $e(t) = A(z)/C(z) y(t)$:
  $ #rect(inset: 4pt, stroke: 0.5pt)[$hat(y)(t+k|t) = F(z)/C(z) y(t)$] $

  === Prediction error
  $ epsilon(t+k|t) = E(z) e(t+k) = sum_(i=0)^(k-1) e_i e(t+k-i) $
  $ "Var"[epsilon] = (e_0^2 + e_1^2 + dots.c + e_(k-1)^2) lambda^2 $

  #bb[Bounds:] $lambda^2 <= "Var"[epsilon(t+k|t)] <= "Var"[y(t)]$

  === Quick results

  #table(
    columns: (1fr, 1fr, 1fr),
    table.header([*Model*], [*$hat(y)(t+1|t)$*], [*Var*]),
    [AR(1)], [$a y(t)$], [$lambda^2$],
    [MA(1)], [$c/(1+c z^(-1)) y(t)$], [$lambda^2$],
  )

  === Non-zero mean
  $ hat(y)(t+k|t) = F(z)/C(z) y(t) + (1 - F(1)/C(1)) m_y $

  === ARMAX $k$-step
  $ hat(y)(t+k|t) = F(z)/C(z) y(t) + (E(z)B(z))/C(z) z^(-k_0) u(t+k) $

  ARMAX 1-step ($E=1$, $F = C - A$):
  $ hat(y)(t+1|t) = (C(z)-A(z))/C(z) y(t) + B(z)/C(z) z^(-k_0) u(t+1) $


  == Model Identification



  #hline()

  === Procedure

  + Design experiment (input, $N$ samples)
  + Select model class $cal(M)(theta)$
  + Define cost: $cal(J)_N (theta) = 1/N sum_(t=1)^N epsilon(t, theta)^2$
  + Minimize: $hat(theta)_N = argmin_theta cal(J)_N (theta)$
  + Validate model

  === ARX vs ARMAX

  #table(
    columns: (auto, 1fr, 1fr),
    table.header([], [*ARX (Least Squares)*], [*ARMAX (Iterative)*]),
    [*Pred.\ error*], [$epsilon = y(t) - phi(t)^TT theta$], [$epsilon = A(z)/C(z) y(t) - B(z)/C(z) u(t-d)$],
    [*Solution*], [Closed-form (OLS)], [Iterative (Newton / gradient)],
    [*Formula*],
    [$hat(theta) = [sum phi phi^TT]^(-1) sum phi y$],
    [$theta^((i+1)) = theta^((i)) - H^(-1) nabla cal(J)$],

    [*Bias*], [Unbiased (if true ARX)], [Asympt. unbiased],
    [*Variance*], [$lambda^2 [sum phi phi^TT]^(-1)$], [—],
    [*Pitfall*], [Biased if true ARMAX], [Local minima only],
    [*Init*], [—], [Use ARX LS ($C(z)=1$)],
  )

  Regressor: $phi(t) = mat(y(t-1), dots.c, y(t-m), u(t-d), dots.c, u(t-d-p+1))^TT$

  Noise estimate: $hat(lambda)^2 = cal(J)_N (hat(theta)_N)$

  === Persistent Excitation

  $u(t)$ is PE of order $n$ iff $R_N = 1/N sum phi(t) phi(t)^TT succ 0$

  #table(
    columns: (1fr, auto),
    table.header([*Input type*], [*PE order*]),
    [White noise], [Any order],
    [Sum of $k$ sinusoids], [$2k$],
    [Single sinusoid at $omega_0$], [$2$ only],
  )

  == Asymptotic Analysis (Infinite Samples)
  To identify the optimal parameters $theta^*$ by minimizing the theoretical prediction error variance $bar(J) = EE[epsilon(t|t-1)^2]$ for a family model $cal(M)$:
  $ epsilon(t|t-1) = y(t) - hat(y)(t|t-1) = A_m(z)/C_m(z) y(t) $
  Find $theta^* = argmin_theta bar(J)$ and $lambda^2_* = bar(J)(theta^*)$.


  == Model Validation

  #hline()

  === Anderson's whiteness test

  $
    hat(rho)_epsilon (tau) = hat(gamma)_epsilon (tau) / hat(gamma)_epsilon (0), quad hat(gamma)_epsilon (tau) = 1/N sum_(t=1)^(N-tau) epsilon(t) epsilon(t+tau)
  $

  Under $H_0$ (WN): $hat(rho)_epsilon (tau) tilde cal(N)(0, 1\/N)$

  #bb[95% confidence band:]
  $ |hat(rho)_epsilon (tau)| <= 1.96 / sqrt(N) quad forall tau = 1, dots, tau_max $

  === Cross-correlation test (ARMAX)

  $ hat(rho)_(epsilon u)(tau) tilde cal(N)(0, 1\/N) $
  Checks correctness of $B(z)/A(z)$ part.

  === Model order selection

  #table(
    columns: (auto, 1fr),
    table.header([*Criterion*], [*Formula*]),
    [FPE], [$"FPE"(n) = (N+n)/(N-n) cal(J)_N (hat(theta))$],
    [AIC], [$"AIC"(n) = ln(cal(J)_N) + 2n/N$],
    [MDL], [$"MDL"(n) = ln(cal(J)_N) + (n ln N)/(2N)$],
  )

  MDL penalizes more → simpler models. #h(1fr) Choose model with *lowest* criterion value.


  == Time Series Analysis — COR/PARCOR

  #hline()

  === Decision table

  #table(
    columns: (auto, 1fr, 1fr),
    table.header([*Model*], [*COR $hat(rho)(tau)$*], [*PARCOR $hat(alpha)(n)$*]),
    [MA($n$)], [Cuts off at $tau = n$], [Tails off],
    [AR($m$)], [Tails off], [Cuts off at $n = m$],
    [ARMA], [Tails off], [Tails off],
  )

  #bb[PARCOR:] $alpha(n) = a_n^((n))$ (last coeff of AR($n$) fit) — use Durbin-Levinson recursion.
  #colbreak()
  === Durbin-Levinson recursion

  $ a_1^((1)) = gamma(1)/gamma(0), quad lambda^2_((1)) = gamma(0)(1 - rho(1)^2) $

  $ a_(n+1)^((n+1)) = 1/lambda^2_((n)) (gamma(n+1) - sum_(i=1)^n a_i^((n)) gamma(n+1-i)) $
  $ a_i^((n+1)) = a_i^((n)) - a_(n+1)^((n+1)) a_(n+1-i)^((n)) $
  $ lambda^2_((n+1)) = (1 - (a_(n+1)^((n+1)))^2) lambda^2_((n)) $

  === Complete workflow

  + Preprocess: remove trend ($hat(k)t + hat(m)$ via OLS), deseasonalize
  + COR → if cuts off at $n_c$ → try MA($n_c$)
  + PARCOR → if cuts off at $n_a$ → try AR($n_a$)
  + Neither cuts off → try ARMA
  + Estimate params (LS for AR, iterative for MA/ARMA)
  + Validate: whiteness test + FPE/AIC/MDL
  + Select best model


  == Non-parametric Identification

  #hline()

  #table(
    columns: (auto, 1fr, auto, auto),
    table.header([*Estimator*], [*Formula*], [*Correct*], [*Consistent*]),
    [Sample mean], [$hat(mu)_N = 1/N sum y(t)$], [Yes], [Yes],
    [Sample cov.], [$hat(gamma)_N (tau) = 1/(N-|tau|) sum y(t) y(t+tau)$], [Yes], [Yes],
    [Periodogram], [$hat(Gamma)^"per" = 1/N |sum y(t) e^(-j omega t)|^2$], [Asymp.], [*No*],
  )

  #bb[Bartlett:] split into $K$ segments, average periodograms → variance $div K$.


  == Data Preprocessing

  #hline()

  #table(
    columns: (auto, 1fr),
    table.header([*Step*], [*Method*]),
    [Detrend], [Fit $hat(k)t + hat(m)$ via OLS, subtract: $tilde(y)(t) = y(t) - (hat(m) + hat(k)t)$],
    [Deseasonalize], [$hat(s)(t) = 1/M sum_(h=0)^(M-1) y(t + h T)$, subtract, model residual],
    [ARIMA($m,d,n$)], [Differencing $Delta = 1 - z^(-1)$, then $A(z) Delta^d y(t) = C(z) e(t)$],
    [Sampling], [$T_s <= pi / omega_max$, rule of thumb: $T_s approx T_"rise" / 10$],
  )
  #colbreak()

  == Box-Jenkins Method

  #hline()

  From raw data → production model in 3 phases:

  === Phase 1: Identification

  + #bb[Stationarity check:]
    - Plot time series, inspect visually
    - Augmented Dickey-Fuller test: $p < 0.05$ → reject non-stationarity
    - If non-stationary: apply differencing ($Delta y(t) = y(t) - y(t-1)$) or log transform
  + #bb[Order selection ($p,q$):]
    - Plot ACF and PACF
    - ACF cuts off at $n$ → MA($n$); PACF cuts off at $m$ → AR($m$)
    - Neither → ARMA, choose orders via FPE/AIC/MDL

  === Phase 2: Estimation

  - Fit model coefficients to data (LS for AR, MLE/iterative for MA/ARMA)
  - Compare models via AIC/BIC: lower = better
  - AIC → better predictive models; BIC → simpler explanatory models

  === Phase 3: Model Diagnostics

  - #bb[Residual whiteness:] ACF of residuals within $plus.minus 1.96/sqrt(N)$
  - #bb[MAE:] Mean Absolute Error measures mean difference between predictions and true values.
  #table(
    columns: (auto, 1fr),
    [*Ljung-Box (Q)*], [Prob(Q) > 0.05 $arrow.r$ Fail to reject $H_0$ (residuals are white noise ✓)],
    [*Jarque-Bera (JB)*], [Prob(JB) > 0.05 $arrow.r$ Fail to reject $H_0$ (residuals normally distributed ✓)],
  )
  - If residuals not white ($p < 0.05$) $arrow.r$ go back to Phase 1.
  - If model is OK $arrow.r$ ready for production (prediction/forecasting).
]

