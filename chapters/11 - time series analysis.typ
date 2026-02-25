#import "/prelude.typ": *

= Time series analysis

== Introduction

A *time series* is a sequence of data ordered by a discrete index (usually time). There is no exogenous variable — the goal is to build a model that describes the dynamics of the observed variable from data alone.

#remark[
  The data sequence is assumed to be a portion of a realization of a stationary stochastic process. The model describes the observed variable as the output of a dynamical system fed by a remote (not accessible) white noise:
  #{
    import fletcher: diagram, edge, node
    figure(
      diagram(
        node-shape: "rect",
        node-stroke: 1pt,
        edge((0, 0), "r", "-|>")[$eta(t) tilde WN$],
        node((1, 0))[$cal(S)$],
        edge("r", "-|>")[$v(t)$],
      ),
    )
  }

  Models and identification algorithms:
  - *AR* $arrow.r$ Least Squares (or Yule-Walker)
  - *MA, ARMA* $arrow.r$ Maximum Likelihood (iterative)
]

== Model type and structure selection

#problem[
  Determine whether the most appropriate model is MA, AR, or ARMA, and estimate the appropriate order(s).
]

#remark[
  An exhaustive strategy would estimate all possible $AR(n_a)$, $MA(n_c)$, $ARMA(n_a, n_c)$ up to given orders and compare them using FPE/AIC/MDL.

  A smarter approach exploits properties of:
  - the *correlation function* (COR) — for MA model identification
  - the *partial correlation function* (PARCOR) — for AR model identification
]

=== Using the correlation function for MA identification

#remark[
  For an $MA(n)$ process: $gamma(tau) = 0$ for $|tau| > n$.

  Therefore, if the sample correlation function $hat(rho)(tau) = hat(gamma)(tau) / hat(gamma)(0)$ "cuts off" at lag $n$ (is within the confidence band $plus.minus 1.96 / sqrt(N)$ for $tau > n$), then an $MA(n)$ model is appropriate.

  This analysis also directly provides the model order $n$.
]

== Yule-Walker equations

#definition(title: "Yule-Walker equations for AR(n)")[
  For an $AR(n)$ model, the covariance function satisfies:
  $
    vec(a_1, a_2, dots.v, a_n) = mat(gamma(0), gamma(1), dots, gamma(n-1); gamma(1), gamma(0), dots, gamma(n-2); dots.v, , dots.down, dots.v; gamma(n-1), gamma(n-2), dots, gamma(0))^(-1) vec(gamma(1), gamma(2), dots.v, gamma(n))
  $

  $ lambda^2 = gamma(0) - a_1 gamma(1) - a_2 gamma(2) - dots - a_n gamma(n) $
]

#remark[
  The matrix in the Yule-Walker equations has a *Toeplitz* structure (elements on each diagonal are equal).

  The YW estimates are asymptotically equivalent to LS estimates: as $N -> infinity$, $hat(theta)_N^("YW") -> hat(theta)_N^("LS")$.
]

== Durbin-Levinson algorithm

The Durbin-Levinson algorithm allows estimation of $AR(n)$ parameters *iteratively*, starting from $AR(n-1)$, requiring only scalar inversions. This makes it efficient for model order selection.

#definition(title: "Durbin-Levinson recursion")[
  Given the parameters of an $AR(n)$ model ($a_i^((n))$, $lambda^2_((n))$), the $AR(n+1)$ parameters are:

  $ a_(n+1)^((n+1)) = 1/lambda^2_((n)) (gamma(n+1) - sum_(i=1)^n a_i^((n)) gamma(n+1-i)) $

  $ a_i^((n+1)) = a_i^((n)) - a_(n+1)^((n+1)) a_(n+1-i)^((n)), quad i = 1, 2, dots, n $

  $ lambda^2_((n+1)) = (1 - (a_(n+1)^((n+1)))^2) lambda^2_((n)) $
]

#remark[
  Starting point ($n=1$):
  $ a_1^((1)) = gamma(1)/gamma(0), quad lambda^2_((1)) = gamma(0) - gamma(1)^2/gamma(0) = gamma(0)(1 - rho(1)^2) $
]

#example(title: "From AR(1) to AR(2)")[
  Given $a_1^((1))$ and $lambda^2_((1))$:

  $ a_2^((2)) = 1/lambda^2_((1))(gamma(2) - a_1^((1)) gamma(1)) $

  $ a_1^((2)) = a_1^((1)) - a_2^((2)) a_1^((1)) $

  $ lambda^2_((2)) = lambda^2_((1))(1 - (a_2^((2)))^2) $

  If the process is actually $AR(1)$, then $gamma(2) = a_1^((1)) gamma(1)$, so $a_2^((2)) = 0$ and $lambda^2_((2)) = lambda^2_((1))$ — the AR(2) model cannot improve over AR(1).
]

== PARCOR function

#definition(title: "Partial Correlation (PARCOR) function")[
  The PARCOR function is defined as:
  $ alpha(n) = a_n^((n)), quad n = 1, 2, 3, dots $

  It is the *last coefficient* of each AR model of increasing order, obtained from the Durbin-Levinson recursion.
]

#properties(title: "Key properties")[
  - $|alpha(n)| <= 1 quad forall n$ (necessary for stability)
  - For an $AR(n_0)$ process: $alpha(n) = 0$ for all $n > n_0$
  - The PARCOR function "cuts off" at the true AR order, analogous to how the correlation function cuts off at the true MA order
]

#remark(title: "Using PARCOR for AR model selection")[
  If the sample PARCOR function $hat(alpha)(n)$ is within the confidence band $plus.minus 1.96/sqrt(N)$ for all $n > n_0$, then an $AR(n_0)$ model is appropriate.
]

== Practical workflow for time series analysis

#note-box(title: "Step-by-step procedure")[
  + *Preprocessing*: remove trends, seasonality, check stationarity
  + *Compute sample COR*: if it cuts off at lag $n_c$ $arrow.r$ try $MA(n_c)$
  + *Compute PARCOR* (via Durbin-Levinson): if it cuts off at order $n_a$ $arrow.r$ try $AR(n_a)$
  + *If neither cuts off clearly*: try ARMA models with orders suggested by COR/PARCOR
  + *Estimate parameters*: LS for AR, ML (iterative) for MA/ARMA
  + *Validate*: whiteness test on residuals, FPE/AIC/MDL comparison
  + *Select the best model*: balance fit quality and parsimony
]

#remark(title: "COR and PARCOR summary")[
  #table(
    columns: (auto, auto, auto),
    [*Model*], [*COR $hat(rho)(tau)$*], [*PARCOR $hat(alpha)(n)$*],
    [MA(n)], [Cuts off at $tau = n$], [Tails off],
    [AR(m)], [Tails off], [Cuts off at $n = m$],
    [ARMA(m,n)], [Tails off], [Tails off],
  )
]
