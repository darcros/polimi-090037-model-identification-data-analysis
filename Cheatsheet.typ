// MIDA 1 — Comprehensive Cheatsheet
// Color-coded: Red=Theorems, Orange=Definitions, Green=Properties, Fuchsia=Examples, Teal=Remarks

#import "prelude.typ": *

#set page(flipped: true, paper: "a4", margin: (x: 1.2cm, y: 1.2cm))
#set text(font: "Libertinus Serif", size: 8pt)
#set par(leading: 0.4em, spacing: 0.5em)
#set block(spacing: 0.5em)
#set heading(numbering: none, outlined: false)
#set enum(indent: 0pt, spacing: 0.35em)
#set list(indent: 0pt, spacing: 0.35em)

#show heading.where(level: 1): it => block(
  text(fill: rgb("#003366"), weight: "bold", size: 12pt, it),
  below: 0.4em,
  above: 0.6em,
)
#show heading.where(level: 2): it => block(
  text(fill: rgb("#003366"), weight: "bold", size: 9.5pt, it),
  below: 0.3em,
  above: 0.5em,
)
#show heading.where(level: 3): it => block(
  text(fill: rgb("#0055A4"), weight: "bold", size: 8.5pt, it),
  below: 0.2em,
  above: 0.4em,
)

// ── Compact color boxes matching chapter colors ──
#let cbox(title, body, color) = block(
  stroke: (left: 2.5pt + color),
  inset: (left: 6pt, y: 3pt, right: 4pt),
  width: 100%,
  above: 0.4em,
  below: 0.4em,
  [#text(fill: color, weight: "bold", size: 7.5pt)[#title] #body],
)
#let thm(title, body) = cbox(title, body, red.darken(20%))
#let defn(title, body) = cbox(title, body, orange)
#let prop(title, body) = cbox(title, body, green.darken(40%))
#let exmp(title, body) = cbox(title, body, fuchsia.darken(10%))
#let rmk(title, body) = cbox(title, body, rgb("#118D8D"))
#let lem(title, body) = cbox(title, body, teal.darken(10%))

#let bb(body) = text(weight: "bold", body)
#let hline() = line(length: 100%, stroke: 0.3pt + luma(180))

// Global table styling
#set table(
  stroke: (x: none, y: 0.3pt + luma(200)),
  inset: (x: 5pt, y: 3pt),
  fill: (_, row) => if row == 0 { rgb("#d8e2ee") } else if calc.odd(row) { rgb("#f0f3f8") },
)

#align(center)[
  #text(size: 14pt, weight: "bold", fill: rgb("#003366"))[MIDA 1 — Complete Cheatsheet]\
  #text(size: 8pt)[Model Identification and Data Analysis — All Theorems, Definitions & Formulas]
  #v(0.3em)
]

#columns(3, gutter: 14pt)[

  // ╔══════════════════════════════════════════════╗
  // ║  1. STOCHASTIC PROCESSES                     ║
  // ╚══════════════════════════════════════════════╝

  = Stochastic Processes
  #hline()

  #defn([SP & Weak Description], [
    $m(t) = EE[v(t)]$, #h(1em) $gamma(t_1, t_2) = EE[(v(t_1)-m_1)(v(t_2)-m_2)]$
  ])

  #defn([Stationarity (SSP)], [
    $m(t) = m space forall t$, #h(0.5em) $gamma(t_1, t_2) = gamma(tau)$ depends only on $tau = t_1 - t_2$.
  ])

  #prop([Covariance of SSPs], [
    - $gamma(0) >= 0$ (variance)
    - $|gamma(tau)| <= gamma(0)$ for all $tau != 0$
    - $gamma(-tau) = gamma(tau)$ (even function)
  ])

  #defn([Normalized covariance], [
    $rho(tau) = gamma(tau) / gamma(0)$, #h(0.5em) $|rho(tau)| <= 1$
  ])

  #thm([Gaussian + Weakly stat. $=>$ Strongly stat.], [])

  #defn([Ergodic process], [
    $lim_(N->oo) 1/N sum_(t=1)^N (dot) -> EE[dot]$ — time avgs $=$ ensemble avgs.
  ])

  #thm([Wold decomposition], [
    Any SSP: $v(t) = w(t) + d(t)$, $w(t) = sum_(k=0)^oo c_k e(t-k)$ (PND $= MA(oo)$), $d(t)$ purely deterministic, $gamma_(w d)(tau)=0 space forall tau$.
  ])

  #defn([Shift operator], [
    $z^(-1) x(t) = x(t-1)$, $z x(t) = x(t+1)$. Linear, composable. #h(1em) *Note:* in exercises, do not mix positive with negative operators.
  ])

  #thm([Asymptotic stability], [
    $W(z)$ asymp. stable iff all poles strictly inside unit circle $|z_i| < 1$.
  ])

  #thm([Stationarity of output], [
    If $e(t)$ is SSP and $W(z)$ asymp. stable, then $y(t) = W(z)e(t)$ is SSP.
  ])

  #defn([Discrete-time transfer function], [
    Transfer function of SP is $W(z)$ such that $y(t) = W(z) e(t)$ where $e(t) ~ WN$.
  ])

  #thm([ARMA($m,n$) stationarity], [
    ARMA($m,n$) is SSP iff $W(z) = C(z)/A(z)$ is asymptotically stable (all poles $|z_i|<1$).
  ])

  // ╔══════════════════════════════════════════════╗
  // ║  2. MODELS                                    ║
  // ╚══════════════════════════════════════════════╝

  = Models
  #hline()

  == Model Summary
  #table(
    columns: (auto, auto, 1fr, auto),
    table.header([*Model*], [*$W(z)$*], [*Poles / Zeros*], [*Memory*]),
    [WN], [$1$], [—], [None],
    [MA($n$)], [$C(z)$], [$n$ zeros (all-zeros)], [Finite],
    [AR($m$)], [$1\/A(z)$], [$m$ poles (all-poles)], [$oo$],
    [ARMA], [$C(z)\/A(z)$], [$m$ poles + $n$ zeros], [$oo$],
    [ARMAX], [$C(z)\/A(z)$, $B(z)\/A(z)$], [$m$+$n$+$p$], [$oo$],
  )

  == White Noise

  #defn([WN], [
    $e(t) ~ WN(mu, lambda^2)$: $m_e = mu$, $gamma(0) = lambda^2$, $gamma(tau) = 0 space forall tau!=0$, $Gamma(omega) = lambda^2$.
  ])

  == MA($n$): $y(t) = sum_(j=0)^n c_j e(t-j)$

  #prop([MA properties], [
    - $m_y = mu sum_(i=0)^n c_i$ #h(1em) (monic: $c_0=1$, $n+1$ free params)
    - $gamma(tau) = lambda^2 sum_(i=0)^(n-|tau|) c_i c_(i+|tau|)$, #h(0.5em) $gamma(tau)=0$ for $|tau|>n$
    - Always strictly stationary
    - MA($oo$) requires $sum c_i^2 < oo$ (square-summable) for well-defined variance
  ])

  #exmp([MA(1) spectrum], [
    $Gamma(omega) = (1+c^2) lambda^2 + 2c lambda^2 cos(omega)$
  ])

  == AR($m$): $y(t) = sum_(i=1)^m a_i y(t-i) + e(t)$

  #prop([AR properties], [
    - Stable iff all roots of $A(z)$ strictly inside unit disk
    - $EE[y(t-tau)e(t)] = 0 space forall tau > 0$
    - $m_y = 0$ (for zero-mean WN input)
  ])

  #defn([Yule-Walker equations], [
    $gamma(tau) = sum_(i=1)^m a_i gamma(tau-i)$ for $|tau| >= 1$.\
    Matrix: $va(a) = Gamma^(-1) va(gamma)$, $lambda^2 = gamma(0) - sum_(i=1)^m a_i gamma(i)$
  ])

  #exmp([AR(1) results], [
    $gamma(0) = lambda^2/(1-a^2)$, $gamma(tau) = a dot gamma(tau-1)$\
    $Gamma(omega) = lambda^2 / (1+a^2-2a cos omega)$
  ])

  #thm([AR(1) as MA($oo$)], [
    $1/(1-a z^(-1)) = sum_(k=0)^oo a^k z^(-k)$ converges iff $|a|<1$.
  ])

  == ARMA($m,n$)

  #defn([ARMA], [
    $y(t) = sum_(i=1)^m a_i y(t-i) + sum_(j=0)^n c_j e(t-j)$, $W(z) = C(z)/A(z)$
  ])

  #exmp([ARMA(1,1)], [
    $gamma(0) = ((1+c_1^2)(1-a^2)+2a c_1)/(1-a^2) lambda^2$\
    $gamma(1) = a gamma(0) + c_1 lambda^2$, $gamma(tau) = a gamma(tau-1)$ for $|tau|>=2$
  ])

  #thm([ARMA as MA($oo$)], [
    Any stable ARMA: $v(t) = sum_(k=0)^oo w_k e(t-k)$ via long division $C(z) div A(z)$.
  ])

  == ARMAX($m,n,p,k_0$)

  #defn([ARMAX], [
    $y(t) = B(z)/A(z) z^(-k_0) u(t) + C(z)/A(z) e(t)$\
    $theta = (a_1,...,a_m, b_0,...,b_(p-1), c_1,...,c_n)^TT$
  ])

  #defn([Unbiased model], [
    $tilde(y)(t) = y(t) - m_y$, $tilde(e)(t) = e(t) - mu$. Then $gamma_y(tau) = gamma_(tilde(y))(tau)$.
  ])

  // ╔══════════════════════════════════════════════╗
  // ║  3. FREQUENCY DOMAIN                         ║
  // ╚══════════════════════════════════════════════╝

  = Frequency Domain
  #hline()

  #defn([Spectral density], [
    $Gamma(omega) = sum_(tau=-oo)^(+oo) gamma(tau) e^(-j omega tau)$
  ])

  #prop([Spectrum properties], [
    Real, $Gamma(omega) >= 0$, even, $2pi$-periodic.
    Exists iff $sum |gamma(tau)| < oo$.
  ])

  #defn([Inverse DFT], [
    $gamma(tau) = 1/(2pi) integral_(-pi)^pi Gamma(omega) e^(j omega tau) dif omega$
  ])

  #thm([Variance as area], [
    $gamma(0) = "Var"[y] = 1/(2pi) integral_(-pi)^pi Gamma(omega) dif omega$
  ])

  #lem([Shannon bound], [
    $T_min = 2$ samples, $omega_max = pi$.
  ])

  #prop([Linearity of spectrum], [
    Scalar: $y = a x => Gamma_y = a^2 Gamma_x$.\
    Uncorr. sum: $z = x+y => Gamma_z = Gamma_x + Gamma_y$.
  ])

  #prop([Euler formula identities], [
    $e^(j omega) = cos(omega) + j sin(omega)$\
    $e^(-j omega) = cos(omega) - j sin(omega)$\
    $e^(-j omega) + e^(j omega) = 2cos(omega)$\
    $e^(j omega) - e^(-j omega) = 2j sin(omega)$
  ])

  #thm([Spectral factorization], [
    If $y(t) = W(z) v(t)$, $W(z)$ stable:\
    $Gamma_y(omega) = |W(e^(j omega))|^2 Gamma_v(omega)$
  ])

  #thm([Theorem of the gain], [
    $EE[y(t)] = W(1) dot EE[u(t)]$ — holds even outside canonical form.
  ])

  #rmk([De-biasing procedure], [
    + Find $EE[y] = W(1) mu$
    + Define $tilde(y) = y - EE[y]$, $tilde(e) = e - mu$
    + $gamma_(tilde(y))(tau) = gamma_y(tau)$
  ])

  #defn([Complex spectrum], [
    $Phi(z) = sum gamma(tau) z^(-tau)$, $Gamma(omega) = Phi(e^(j omega))$
  ])

  #defn([Cross-spectrum], [
    $Phi_(u y)(z) = W(z) Phi_(u u)(z)$\
    $Phi_(y y)(z) = W(z) W(z^(-1)) Phi_(u u)(z)$
  ])

  #rmk([Graphical method], [
    $overline(a)$, $overline(b)$ = distances from zero/pole to $e^(j omega)$:\
    $Gamma_y(omega) = |overline(a)/overline(b)|^2 lambda^2$.\
    Poles near unit circle $=>$ high spectrum; zeros near $=>$ low.
  ])

  // ╔══════════════════════════════════════════════╗
  // ║  4. CANONICAL REPRESENTATION                  ║
  // ╚══════════════════════════════════════════════╝

  = Canonical Representation
  #hline()

  #thm([Canonical form], [
    $y(t) = hat(W)(z)e(t) = C(z)/A(z) e(t)$ is canonical iff:
    + *Monic*: $C(z), A(z)$ have leading coeff $= 1$
    + *Same degree*: $deg C = deg A$ (null relative degree)
    + *Coprime*: no common factors
    + *Min. phase*: roots of $C(z)$ in $|z| <= 1$; roots of $A(z)$ in $|z| < 1$
  ])

  #prop([Sources of redundancy (4 transformations)], [
    + Scale by $alpha/alpha$: $W -> W/alpha$, noise $-> WN(0, alpha^2 lambda^2)$ — resolved by *Monic* condition
    + Multiply by $z^n / z^n$: shifts noise in time — resolved by *Null relative degree*
    + Multiply by $(z-p)/(z-p)$ with $|p|<1$: adds removable factors — resolved by *Coprime* condition
    + Multiply by all-pass $T(z)$ with $|q|>1$: moves poles in/out unit circle — resolved by *Min. phase*

    Removing all redundancy yields *unique canonical representation*.
  ])

  #thm([All-pass filter], [
    $T(z) = 1/q dot (z-q)/(z-1/q)$, $|q|>1$. Then $T(z)T(z^(-1))=1$.\
    Preserves spectrum: $tilde(Phi)(z) = Phi(z)$.
  ])

  #exmp([Canonical form example], [
    ARMA(1,1) with $|c| <= 1$: already canonical.\
    If $|c|>1$: use $hat(W)(z) = (1+1/c dot z^(-1))/(1-a z^(-1))$, $hat(e) ~ WN(0, c^2 lambda^2)$.
  ])

  // ╔══════════════════════════════════════════════╗
  // ║  5. PREDICTION                                ║
  // ╚══════════════════════════════════════════════╝

  = Prediction
  #hline()

  Given $y(t) = C(z)/A(z) e(t)$ in *canonical form*, $e(t) ~ WN(0, lambda^2)$.

  == Step 1 — Diophantine equation

  #defn([Long division], [
    Divide $C(z)$ by $A(z)$ for $k$ steps:
    $C(z) = E(z) A(z) + z^(-k) F(z)$\
    $E(z)$: degree $k-1$. $F(z)$: degree $max(m, n)-1$.
  ])

  == Step 2 — Predictor

  #defn([Predictor from noise], [$hat(y)(t+k|t) = F(z)/A(z) e(t)$])

  #defn([Whitening filter], [$e(t) = A(z)/C(z) y(t)$])

  #defn([Predictor from data], [
    $#rect(inset: 3pt, stroke: 0.4pt)[$ hat(y)(t+k|t) = F(z)/C(z) y(t) $]$
  ])

  == Prediction Error

  #prop([Prediction error], [
    $epsilon(t+k|t) = E(z) e(t+k) = sum_(i=0)^(k-1) e_i e(t+k-i)$\
    $"Var"[epsilon] = (e_0^2 + e_1^2 + dots + e_(k-1)^2) lambda^2$\
    *Bounds:* $lambda^2 <= "Var"[epsilon] <= "Var"[y(t)]$
  ])

  == 1-step Shortcuts ($k=1$: $E=1$, $F=C-A$)

  #prop([1-step ARMA predictor], [
    From noise: $hat(y)(t|t-1) = (C(z)-A(z))/A(z) e(t)$\
    From data: $hat(y)(t|t-1) = (C(z)-A(z))/C(z) y(t)$\
    Pred. error (whitening): $epsilon(t|t-1) = A(z)/C(z) y(t) = e(t)$
  ])
  #prop([ARMA($m,n$) mean & representation], [
    Mean: $m_y = frac(sum_(i=0)^n c_i, 1 - sum_(i=1)^m a_i) m_e$
    Any AR($m$) with $|a_i|<1$ $=>$ MA($oo$) with $c_k = a^k$ (asymptotic representation)
  ])
  == Quick Reference

  #table(
    columns: (auto, 1fr, auto),
    table.header([*Model*], [*$hat(y)(t+1|t)$*], [*Var*]),
    [AR(1)], [$a y(t)$], [$lambda^2$],
    [AR($n$)], [$a_1 y(t-1) + dots + a_n y(t-n)$], [$lambda^2$],
    [MA(1)], [$c/(1+c z^(-1)) y(t)$], [$lambda^2$],
    [AR(1) $r$-step], [$a^r y(t)$ → $EE[y]$ as $r->oo$], [grows],
  )

  == Non-zero Mean

  #prop([Prediction with bias], [
    $hat(y)(t+k|t) = F(z)/C(z) y(t) + (1 - F(1)/C(1)) m_y$
  ])

  == ARMAX Prediction

  #prop([ARMAX $k$-step], [
    $hat(y)(t+k|t) = F(z)/C(z) y(t) + (E(z) B(z))/C(z) z^(-k_0) u(t+k)$\
    1-step: $hat(y)(t+1|t) = (C(z)-A(z))/C(z) y(t) + B(z)/C(z) z^(-k_0) u(t+1)$
  ])


  // ╔══════════════════════════════════════════════╗
  // ║  6. MODEL IDENTIFICATION                      ║
  // ╚══════════════════════════════════════════════╝

  = Model Identification
  #hline()

  == Procedure
  + Design experiment (input $u$, $N$ samples)
  + Select model class $cal(M)(theta)$
  + Define cost: $cal(J)_N(theta) = 1/N sum_(t=1)^N epsilon(t, theta)^2$
  + Minimize: $hat(theta)_N = argmin_theta cal(J)_N(theta)$
  + Validate model

  #defn([PEM framework], [
    $epsilon(t, theta) = y(t) - hat(y)(t|t-1,theta)$\
    ARX: $epsilon = A(z)y(t) - B(z)u(t-d)$ — linear in $theta$\
    ARMAX: $epsilon = A(z)/C(z) y(t) - B(z)/C(z) u(t-d)$ — nonlinear in $theta$
  ])

  == ARX Identification (Least Squares)

  #defn([Regressor], [
    $phi(t) = [y(t-1), dots, y(t-m), u(t-d), dots, u(t-d-p+1)]^TT$\
    Predictor: $hat(y)(t|t-1) = phi(t)^TT theta$
  ])

  #thm([OLS formula], [
    $hat(theta)_N = [sum_(t=1)^N phi(t)phi(t)^TT]^(-1) sum_(t=1)^N phi(t) y(t)$\
    Hessian: $2/N sum phi(t)phi(t)^TT gt.eq.slant 0$ (always PSD).
  ])

  #thm([Unbiasedness of LS], [
    If true system is ARX: $EE[hat(theta)_N] = theta_0$
  ])

  #thm([Variance of LS], [
    $"Var"[hat(theta)_N] = lambda^2 [sum_(t=1)^N phi(t) phi(t)^TT]^(-1)$ — decreases as $1/N$.
  ])

  #rmk([Noise variance estimate], [
    $hat(lambda)^2 = cal(J)_N(hat(theta)_N) = 1/N sum epsilon(t, hat(theta)_N)^2$
  ])

  == Persistent Excitation

  #defn([PE of order $n$], [
    $R_N = 1/N sum phi(t) phi(t)^TT succ 0$ for large $N$.
  ])

  #table(
    columns: (1fr, auto),
    table.header([*Input type*], [*PE order*]),
    [White noise], [Any],
    [Sum of $k$ sinusoids], [$2k$],
    [Single sinusoid], [$2$ only],
    [PRBS], [High order],
  )

  #defn([Identifiability], [
    + *Structural*: different $theta$ $=>$ different models
    + *Experimental*: input PE of sufficient order
  ])

  == ARMAX Identification (Iterative PEM)

  #rmk([Key difference], [
    $hat(y)(t|t-1)$ is nonlinear in $theta$ (through $C(z)$) $=>$ $cal(J)_N$ non-quadratic $=>$ multiple local minima $=>$ iterative methods required.
  ])

  #prop([Quasi-Newton update], [
    $theta^((i+1)) = theta^((i)) - [sum pdv(epsilon, theta)(pdv(epsilon, theta))^TT]^(-1) [sum epsilon dot pdv(epsilon, theta)]$
  ])

  #prop([Gradient], [
    $pdv(cal(J)_N, theta) = 2/N sum epsilon(t) pdv(epsilon(t), theta)$
  ])

  #prop([Approx. Hessian (always PSD)], [
    $pdv(cal(J)_N, theta, 2) approx 2/N sum pdv(epsilon, theta)(pdv(epsilon, theta))^TT$
  ])

  #defn([Auxiliary signals], [
    $alpha(t) = -1/C(z) y(t)$, $beta(t) = -1/C(z) u(t)$, $gamma(t) = -1/C(z) epsilon(t)$\
    $pdv(epsilon, theta) = vec(alpha(t-1), dots, alpha(t-m), beta(t-d), dots, beta(t-d-p+1), gamma(t-1), dots, gamma(t-n))$
  ])

  #rmk([Initialization], [
    Start with ARX LS estimate (set $C(z)=1$), or use multiple random inits.
  ])

  == ARX vs ARMAX Comparison

  #table(
    columns: (auto, 1fr, 1fr),
    table.header([], [*ARX (LS)*], [*ARMAX (Iterative)*]),
    [*Solution*], [Closed-form (OLS)], [Iterative (Newton)],
    [*Bias*], [Unbiased (if true ARX)], [Asymp. unbiased],
    [*Variance*], [$lambda^2 [sum phi phi^TT]^(-1)$], [—],
    [*Pitfall*], [Biased if true ARMAX], [Local minima],
    [*Init*], [—], [Use ARX LS],
  )

  == Asymptotic Analysis ($N -> oo$)

  #thm([Asymptotic PEM], [
    $hat(theta)_N ->_(N->oo) theta^*$ where $theta^*$ min. $overline(cal(J))(theta) = EE[epsilon(t, theta)^2]$.\
    If $cal(S) in cal(M)$: $theta^* = theta_0$.
  ])

  #thm([PEM convergence ($cal(S) in cal(M)$)], [
    $EE[epsilon(t, theta)^2] = lambda^2 + EE[(hat(y)(theta^0) - hat(y)(theta))^2] >= lambda^2$\
    Equality iff $theta = theta^0$ $=>$ $hat(theta)_N -> theta^0$.
  ])

  #rmk([Innovation form], [
    At $theta^*$: $epsilon(t, theta_0) = e(t) ~ WN(0, lambda^2)$ — basis for validation.
  ])

  == Four Cases of PEM Convergence

  #table(
    columns: (auto, 1fr, 1fr),
    table.header([], [$Delta$ singleton], [$Delta$ not singleton]),
    [$cal(S) in cal(M)$], [$hat(theta)_N -> theta^0$ (unique)], [One of $Delta$ (over-param.)],
    [$cal(S) in.not cal(M)$], [$hat(theta)_N -> theta^*$ (best proxy)], [No guarantee],
  )

  == Variance & Confidence

  #thm([Variance of PEM estimates], [
    Asymptotically: $sqrt(N)(hat(theta)_N - theta_0) -> cal(N)(0, P)$\
    CI: $hat(theta)_(N,i) plus.minus z_(alpha\/2) sqrt("Var"[hat(theta)_(N,i)])$
  ])

  #rmk([Bias in ARX with colored noise], [
    If true system is ARMAX but we fit ARX: LS is *biased* because $C(z)/A(z) e(t)$ correlates with regressor.
  ])

  // ╔══════════════════════════════════════════════╗
  // ║  7. MODEL VALIDATION                          ║
  // ╚══════════════════════════════════════════════╝

  = Model Validation
  #hline()

  #defn([Anderson's whiteness test], [
    $hat(rho)_epsilon(tau) = hat(gamma)_epsilon(tau)/hat(gamma)_epsilon(0)$, under $H_0$ (WN): $hat(rho) ~ cal(N)(0, 1/N)$\
    *95% band:* $|hat(rho)_epsilon(tau)| <= 1.96/sqrt(N)$ for $tau = 1, dots, T_max$
  ])

  #defn([Cross-correlation test (ARMAX)], [
    $hat(rho)_(epsilon u)(tau) ~ cal(N)(0, 1/N)$ — checks $B(z)/A(z)$ part.\
    If significant: wrong $B(z)/A(z)$, wrong delay, missing dynamics.
  ])

  #rmk([Which test checks what], [
    Anderson (auto-corr.) $=>$ checks noise model $C(z)/A(z)$.\
    Cross-corr. $=>$ checks I/O model $B(z)/A(z)$.
  ])

  == Model Order Selection

  #table(
    columns: (auto, 1fr),
    table.header([*Criterion*], [*Formula*]),
    [FPE], [$"FPE"(n) = (N+n)/(N-n) cal(J)_N(hat(theta))$],
    [AIC], [$"AIC"(n) = ln(cal(J)_N) + 2n/N$],
    [MDL], [$"MDL"(n) = ln(cal(J)_N) + (n ln N)/(2N)$],
  )

  #table(
    columns: (auto, 1fr),
    table.header([*Criterion*], [*When to use*]),
    [FPE/AIC], [When $S in cal(M)$ or prefer flexibility],
    [MDL], [True system in model class; more conservative],
  )

  #rmk([Selection rules], [
    MDL penalizes more $=>$ simpler models. FPE and AIC *asymptotically equivalent* for large $N$. Choose model with *lowest* criterion.
  ])

  #prop([Why tests fail in practice], [
    + *Finite-sample correlation:* Spurious correlations vanish only as $N -> oo$
    + *Estimation error:* $hat(theta)_N != theta^0$ exactly; residuals retain true system structure
    + *Underspecification:* True system $S in.not cal(M)$; residuals contain lost information
    + *Non-stationarity:* Real data violates constant statistics assumption
    + *Multiple comparisons:* Testing $T_max$ lags; expect $approx 0.05 T_max$ false positives
  ])

  == Cross-validation

  #rmk([Train/validation split], [
    Partition data into identification (training) and validation sets. Standard $k$-fold CV *not suitable* for time series (breaks temporal corr.). Use time-series-aware: forward-chaining or rolling windows.
  ])


  // ╔══════════════════════════════════════════════╗
  // ║  8. NON-PARAMETRIC IDENTIFICATION             ║
  // ╚══════════════════════════════════════════════╝

  = Non-parametric Identification
  #hline()

  #defn([Correctness], [$EE[hat(T)_N] = T^*$. Asymptotic: $EE[hat(T)_N] ->_(N->oo) T^*$.])

  #defn([Consistency], [$"Var"[hat(T)_N] ->_(N->oo) 0$])

  #table(
    columns: (auto, 1fr, auto, auto),
    table.header([*Estimator*], [*Formula*], [*Correct*], [*Consistent*]),
    [Sample mean], [$hat(mu)_N = 1/N sum y(t)$], [Yes], [Yes],
    [Sample cov.], [$hat(gamma)_N(tau) = 1/(N-|tau|) sum y(t)y(t+tau)$], [Yes], [Yes],
    [Alt. cov.], [$hat(gamma)'_N(tau) = 1/N sum (y-hat(mu))(y_tau - hat(mu))$], [Biased], [Yes],
    [Periodogram], [$hat(Gamma)^"per" = 1/N |sum y(t)e^(-j omega t)|^2$], [Asymp.], [*No*],
  )
  #rmk([Alt. covariance guarantee], [
    $hat(gamma)'_N(tau)$ biased but guarantees positive semidefinite Toeplitz matrix (useful for factorization).
  ])
  #rmk([Periodogram inconsistency], [
    Variance $arrow.r.not 0$ as $N -> oo$. Remains $prop Gamma(omega)^2$.
  ])

  #defn([Bartlett's method], [
    Split into $K$ segments of length $M = N/K$, average periodograms $=>$ variance $div K$, but lower freq. resolution.
  ])

  #defn([Windowing method], [
    $hat(Gamma)^"win"(omega) = sum w(tau) hat(gamma)(tau) e^(-j omega tau)$\
    Window: $w(0)=1$, symmetric, finite support.
  ])

  // ╔══════════════════════════════════════════════╗
  // ║  9. DATA PREPROCESSING                        ║
  // ╚══════════════════════════════════════════════╝

  = Data Preprocessing
  #hline()

  #table(
    columns: (auto, 1fr),
    table.header([*Step*], [*Method*]),
    [Detrend], [Fit $hat(k)t + hat(m)$ via OLS, subtract: $tilde(y)(t) = y(t)-(hat(m)+hat(k)t)$],
    [Deseason.], [$hat(s)(t) = 1/M sum_(h=0)^(M-1) y(t+h T)$; subtract, model residual],
    [ARIMA], [$Delta = 1-z^(-1)$; $A(z)Delta^d y(t) = C(z)e(t)$],
    [Sampling], [$T_s <= pi/omega_max$; rule of thumb: $T_s approx T_"rise"/10$],
    [Outliers], [$> plus.minus 3 sigma$: remove/interpolate, robust loss, or trim],
  )

  #defn([ARIMA($m,d,n$)], [
    $A(z) Delta^d y(t) = C(z) e(t)$. Most common: $d=1$.
  ])

  #defn([CARIMA/ARIMAX], [
    $A(z) Delta y(t) = B(z) u(t-d) + C(z)/Delta e(t)$
  ])

  // ╔══════════════════════════════════════════════╗
  // ║  10. TIME SERIES ANALYSIS                     ║
  // ╚══════════════════════════════════════════════╝

  = Time Series Analysis
  #hline()

  == COR/PARCOR Decision Table

  #table(
    columns: (auto, 1fr, 1fr),
    table.header([*Model*], [*COR $hat(rho)(tau)$*], [*PARCOR $hat(alpha)(n)$*]),
    [MA($n$)], [Cuts off at $tau = n$], [Tails off],
    [AR($m$)], [Tails off], [Cuts off at $n = m$],
    [ARMA], [Tails off], [Tails off],
  )

  #defn([PARCOR function], [
    $alpha(n) = a_n^((n))$ — last coeff of AR($n$) fit from Durbin-Levinson.\
    $|alpha(n)| <= 1$. For $AR(n_0)$: $alpha(n)=0$ for $n > n_0$.
  ])

  == Durbin-Levinson Recursion

  #defn([Durbin-Levinson], [
    *Init:* $a_1^((1)) = gamma(1)/gamma(0)$, $lambda^2_((1)) = gamma(0)(1-rho(1)^2)$\
    *Update:*\
    $a_(n+1)^((n+1)) = 1/lambda^2_((n)) (gamma(n+1) - sum_(i=1)^n a_i^((n)) gamma(n+1-i))$\
    $a_i^((n+1)) = a_i^((n)) - a_(n+1)^((n+1)) a_(n+1-i)^((n))$\
    $lambda^2_((n+1)) = (1-(a_(n+1)^((n+1)))^2) lambda^2_((n))$
  ])

  == Complete Workflow

  #prop([Step-by-step], [
    + Preprocess: detrend ($hat(k)t+hat(m)$ via OLS), deseasonalize
    + COR $=>$ cuts at $n_c$ $=>$ try MA($n_c$)
    + PARCOR $=>$ cuts at $n_a$ $=>$ try AR($n_a$)
    + Neither cuts off $=>$ try ARMA
    + Estimate params (LS for AR, iterative for MA/ARMA)
    + Validate: whiteness test + FPE/AIC/MDL
    + Select best model
  ])

  // ╔══════════════════════════════════════════════╗
  // ║  11. BOX-JENKINS METHOD                       ║
  // ╚══════════════════════════════════════════════╝

  = Box-Jenkins Method
  #hline()

  #prop([Phase 1: Identification], [
    + Visual stationarity check + ADF test ($p < 0.05$ $=>$ stationary)
    + Non-stationary $=>$ differencing ($Delta y$) or log transform
    + ACF cuts at $n$ $=>$ MA($n$); PACF cuts at $m$ $=>$ AR($m$); neither $=>$ ARMA
  ])

  #rmk([Phase 2: Estimation], [
    LS for AR, iterative for MA/ARMA. Compare: AIC (predictive) vs MDL/BIC (parsimonious). Lower $=$ better.
  ])

  #table(
    columns: (auto, 1fr),
    table.header([*Test*], [*Interpretation*]),
    [Ljung-Box (Q)], [Prob(Q) $> 0.05$ $=>$ residuals are WN ✓],
    [Jarque-Bera (JB)], [Prob(JB) $> 0.05$ $=>$ residuals normal ✓],
    [Residual ACF], [Within $plus.minus 1.96/sqrt(N)$ ✓],
  )

  #rmk([Decision], [
    Residuals not white ($p < 0.05$) $=>$ back to Phase 1. Model OK $=>$ forecast.
  ])

  // ╔══════════════════════════════════════════════╗
  // ║  EXERCISE FORMULAS QUICK REFERENCE            ║
  // ╚══════════════════════════════════════════════╝

  = Exercise Formulas
  #hline()

  == Computing $gamma(tau)$ for AR/ARMA

  #exmp([General recipe], [
    + Write $y(t) = a y(t-1) + dots + c_0 e(t) + dots$
    + Multiply both sides by $y(t-tau)$, take $EE[dot]$
    + Use $EE[y(t-tau)e(t)] = 0$ for $tau > 0$
    + For $tau = 0$: get $gamma(0)$ involving $lambda^2$
    + For $tau >= 1$: recursive Yule-Walker equations
    + Solve the system
  ])

  == Spectrum Computation

  #exmp([Two methods], [
    *Analytical:* $Gamma(omega) = |W(e^(j omega))|^2 lambda^2$. Expand $|dot|^2$ using Euler.\
    *Graphical:* Place poles/zeros in $CC$. At each $omega$, measure distances $overline(a)_i$ (from zeros) and $overline(b)_j$ (from poles) to $e^(j omega)$: $Gamma(omega) = (product |overline(a)_i|^2)/(product |overline(b)_j|^2) lambda^2$
  ])

  == Prediction Step-by-Step

  #exmp([Complete recipe], [
    + Verify canonical form (monic, coprime, same degree, min-phase)
    + If not canonical: apply all-pass filter
    + *1-step:* $E=1$, $F=C-A$, $hat(y)(t|t-1) = (C-A)/C y(t)$
    + *$k$-step:* Long division $C(z) div A(z)$ for $k$ steps $=>$ $E(z)$, $F(z)$
    + From data: $hat(y)(t+k|t) = F(z)/C(z) y(t)$
    + Non-zero mean: add bias term $(1-F(1)/C(1))m_y$
    + Error variance: $(sum_(i=0)^(k-1) e_i^2) lambda^2$
  ])

  == Identification Step-by-Step

  #exmp([ARX/ARMAX recipe], [
    *ARX:* Build $phi(t)$ from data, apply OLS formula directly.\
    *ARMAX:*
    + Init $theta^((0))$ from ARX LS (set $C(z)=1$)
    + Compute $epsilon(t, theta^((i)))$ and auxiliary signals $alpha, beta, gamma$
    + Update $theta^((i+1))$ via Quasi-Newton
    + Repeat until $||theta^((i+1)) - theta^((i))|| <$ tol
    + $hat(lambda)^2 = cal(J)_N(hat(theta)_N)$
  ])

  == Validation Step-by-Step

  #exmp([Validation recipe], [
    + Compute residuals $epsilon(t) = y(t) - hat(y)(t|t-1, hat(theta))$
    + Anderson test: $|hat(rho)_epsilon(tau)| <= 1.96/sqrt(N)$?
    + For ARMAX: cross-correlation $|hat(rho)_(epsilon u)(tau)| <= 1.96/sqrt(N)$?
    + Compare models: FPE/AIC/MDL — *lowest wins*
    + If "white" fails: increase model order or change class
  ])

]

