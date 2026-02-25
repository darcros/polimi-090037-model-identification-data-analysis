#import "/prelude.typ": *

= Prediction

#problem(title: "Prediction problem")[
  Given
  - a SSP $y(t)$ in canonical form
    $ y(t) = C(z) / A(z) e(t), quad e(t) ~ WN(mu, lambda^2) $
    with $C(z)$, $A(z)$, $mu$, $lambda^2$ known.

  - a specific realization $y(0), y(1), dots y(t)$

  How can we predict future values of the process at a generic time instant $t+k$ ?


  #figure(
    cetz.canvas({
      import suiji: gen-rng
      import cetz.draw: *
      import cetz-plot: *
      import "../util.typ": random-series

      let t = 10.0
      let tk = 15.0

      let rng = gen-rng(1000)
      let (rng, realization) = random-series(rng, t, from: 0, to: t)

      plot.plot(
        size: (10, 2),
        axis-style: "school-book",
        x-min: -1,
        x-max: tk + 1,
        x-tick-step: none,
        x-ticks: ((t, $t$), (tk, $t+k$)),
        x-label: [$t$],
        y-tick-step: none,
        y-label: [$y(t)$],
        {
          plot.add(label: "realization", realization, line: "spline")
          plot.annotate({
            line((t, 0.0), realization.last(), stroke: (dash: "dashed"))
            content((tk, 0.5), [?])
          })
        },
      )
    }),
  )
]

#definition(title: "Predictor")[
  $hat(y) (t+k|t)$ is the predictor of SSP $y(t)$: an estimate of the future value $y(t+k)$ based on all available observations up to time $t$.
]

#definition(title: "Trivial predictor")[
  $hat(y) (t+k|t) = EE[y(t)] = m_y$
]

#definition(title: "Linear predictor")[
  An "optimal" linear predictor contains both information about the model information $chevron.l alpha_0, alpha_1, dots chevron.r$ and the data $chevron.l y(t), y(t-1), dots chevron.r$.
  $
    hat(y)(t+k|t) & = alpha_0 y(t) + alpha_1 y(t - 1) + dots \
                  & = sum_(i=0)^(+infinity) (alpha_i) y(t-i)
  $

  with $alpha_i in RR$ such that $sum_(i=0)^(+infinity) alpha_i < +infinity$

  Note that $alpha_i$ encodes information from the *model* while $y(t-i)$ provides information from the *data*.
]

#definition(title: "Prediction error")[
  The prediction error of $hat(y)$ is
  $ epsilon(t+k|t) = y(t+k) - hat(y) (t+k|t) $
] <def:prediction-error>

#definition(title: "Mean squared prediction error")[
  The MSP of $hat(y)$ is
  $
    cal(J) (hat(y)) & = EE[epsilon(t+k|t)^2] \
                    & = EE[(y(t+k) - hat(y) (t+k|t))^2]
  $
]

== Optimal predictor

#problem(title: "Optimal predictor problem")[
  Given
  - a SSP $y(t)$ in canonical form
    $ y(t) = C(z) / A(z) e(t), quad e(t) ~ WN(0, lambda^2) $
    with $C(z)$, $A(z)$, $mu$, $lambda^2$ known.

  - a specific realization $y(0), y(1), dots y(t)$

  we would like to find the predictor $hat(y) (t+k|t)$ that minimizes $cal(J)(hat(y))$

  Since we work with linear (ARMA/ARMAX) models, we will restrict ourselves to linear predictors.
  So we want to find the parameters $hat(alpha) = mat(hat(alpha)_1, hat(alpha)_2, dots)^T$ for a linear predictor, such that the mean squared prediction error is minimized
  $ hat(alpha) = argmin_alpha cal(J) (alpha) $
]

#note-box[
  For now we are considering only processes with zero mean.
]

#solution[
  Remember that
  $
      y(t) & = W(z)e(t) \
           & = w_0 e(t) + w_1 e(t-1) + w_2 e(t-2) + dots \
           & = sum_(i=0)^(+infinity) w_i e(t-i) \
    y(t-1) & = sum_(i=0)^(+infinity) w_i e(t-i-1) \
    y(t-2) & = sum_(i=0)^(+infinity) w_i e(t-i-2) \
           & dots.v
  $

  So we can rewrite our linear predictor in the following way
  $
    hat(y)(t+k|t) & = alpha_0 y(t) + alpha_1 y(t - 1) + dots \
                  & = alpha_0 sum_(i=0)^(+infinity) w_i e(t-i) + alpha_1 sum_(i=0)^(+infinity) w_i e(t-i-1) + dots \
                  & = sum_(i=0)^(+infinity) beta_i e(t-i)
  $

  // TODO: what is the definition of beta_i?
  // how exactly do we re-arrange the w_i and alpha_i to get beta_i?

  Then, we reformulate the optimization problem wrt $alpha$ in a new optimization problem wrt $beta$
  $ cal(V) = EE[(y(t+k) - hat(y)(t+k|t, beta))^2] $
  we are now looking for $hat(beta) = argmin_beta cal(V) (beta)$

  Let's consider $y(t+k)$ alone, we can split it into two parts
  $
    y(t+k) & = sum_(i=0)^(+infinity) w_i e(t-i-k) \
           & =
             underbrace(
               sum_(j=0)^(k-1) w_j e(t+k-j),
               "future"\ ("because" k-j > 0)
             ) +
             underbrace(
               sum_(i=0)^(+infinity) w_(i+k) e(t-i),
               "present" ("for" i = 0)\ "and past" ("for" i > 0)
             )
  $

  Now lets consider the full $cal(V)(beta)$:
  #set math.equation(numbering: none)
  #block(width: 100%, inset: 0pt)[
    $
      cal(V)(beta) & = EE[(y(t+k) - hat(y)(t+k|t, beta))^2] \
                   & = EE[( sum_(j=0)^(k-1) w_j e(t+k-j) +
                         sum_(i=0)^(+infinity) w_(k+i) e(t-i) -
                         sum_(i=0)^(+infinity) beta_i e(t-i) )^2] \
                   & = EE[ underbrace((sum_(j=0)^(k-1) w_j e(t+k-j))^2, "not a function of" beta) +
                       (sum_(i=0)^(+infinity) (w_(k+i) - beta_i) e(t-i))^2 + underbrace(
                         cancel(2 (sum_(j=0)^(k-1) w_j e(t+k-j)) (sum_(i=0)^(+infinity) (w_(k+i) - beta_i) e(t-i))),
                         "all uncorrelated because WN at different time instants"
                       ) ]
    $
  ]
  This means that $hat(beta)_i = w_(k+i), i = 0, 1, dots$

  So our optimal predictor is
  $ hat(y)(t+k|t) = sum_(i=0)^(+infinity) w_(k+i) e(t-i) $
]

#remark[
  This solution presents some practical challenges, because:
  - it requires samples of $e(t)$ (we can only sample $y(t)$)
  - it requires the use of an infinite number of coefficients
]

#problem[
  We want to find a predictor $hat(y)$ as before, but with a finite number of coefficients.
]

#solution[
  Consider the operatorial representation of $y(t)$:
  $ y(t) = W(z)e(t) = C(z) / A(z) e(t) $

  We can rewrite the transfer function as
  $
    W(z)
    = C(z) / A(z)
    = (1 + c_1 z^(-1) + c_2 z^(-2) + dots + c_n z^(-n)) / (1 + a_1 z^(-1) + a_2 z^(-2) + dots + a_m z^(-m))
    = w_0 + w_1 z^(-1) + w_2 z^(-2) + dots
  $

  but this sequence is infinite.

  Lets instead consider a finite number $k$ of division steps
  $
    longdivision(
      divisor: A(z),
      quotient: E(z),
      steps: C(z) \ dots.v \ z^(-k) F(z)
    )
  $

  After $k$ division steps we get a quotient $E(z)$ which is a polynomial of order $k$ and a rest $z^(-k) F(z)$
  $ W(z) = C(z) / A(z) = E(z) + (z^(-k) F(z)) / A(z) $

  #definition(title: "Diophantine equation")[
    $ C(z) = E(z)A(z) + z^(-k) F(z) $
    This identity, obtained by rewriting the long division, is called the *Diophantine equation*. $E(z)$ has degree $k-1$ and $F(z)$ has degree $max(m, n) - 1$.
  ]

  If we substitute $W(z) = E(z) + (z^(-k) F(z)) / A(z)$ into $y(t+k)$ we can observe that
  $
    y(t+k) & = W(z) e(t+k) \
           & = (E(z) + (z^(-k) F(z)) / A(z)) e(t+k) \
           & = E(z) e(t+k) + F(z) / A(z) z^(-k) e(t+k) \
           & =
             underbrace(E(z) e(t+k), "(a)") +
             underbrace(F(z) / A(z) e(t), "(b)")
  $

  #[
    #set enum(numbering: "(a)")
    + $E(z) e(t+k)$ depends on *future* samples of $e(t)$ and is unpredictable at time $t$
    + $F(z)/A(z) e(t)$ depends on *past and present* samples of $e(t)$ and is computable at time $t$ --- this is our predictor
  ]

  $ y(t+k) = #[prediction error] + #[prediction at time $t$] $

  #definition(title: "Predictor from noise")[
    Given a SSP $y(t) = W(z)e(t)$
    where
    - $W(z)$ in normal form
    - $e(t) ~ WN(0, lambda^2)$
    the _predictor from noise_ is

    $ hat(y)(t+k|t) = F(z) / A(z) e(t) $

    where $F(z)$, $A(z)$ are from $W(z) = E(z) + (z^(-k) F(z)) / A(z)$.

    $hat(y)$ is an $ARMA$ process.
  ]
]

#remark[
  This predictor is better than the first one as it has a finite number of coefficients, however it is still dependent on samples of the noise $e(t)$ which we don't have.
]

#problem[
  We want to find a predictor $hat(y)$ as before, but
  - with a finite number of coefficients
  - _dependent on samples of $y(t)$, not $e(t)$_
]

#definition(title: "Whitening filter")[
  For a SSP $y(t)$ with transfer function $W(z)$
  $ W(z) = C(z) / A(z) $
  the _whitening filter_ (or _inverse filter_) is
  $ W^(-1)(z) = A(z) / C(z) $

  The whitening filter "de-colors" $y(t)$ and turns it into white noise
  #{
    import fletcher: diagram, edge, node
    figure(
      diagram(
        node-shape: "rect",
        node-stroke: 1pt,
        edge((0, 0), "r", "-|>")[$e(t)$],
        node((1, 0))[$W(z)$],
        edge((1, 0), "rrr", "-|>")[$y(t) = C(z) / A(z) e(t)$],
        node((4, 0))[$X(z)$],
        edge("rrr", "-|>")[$C(z) / A(z) A(z) / C(z) e(t) = e(t)$],
      ),
    )
  }
]

#remark[
  $
    y(t) = C(z) / A(z) e(t) <==> e(t) = A(z) / C(z) y(t)
  $
]

#solution[
  We apply the whitening filter to the predictor $hat(y)$

  $
    hat(y)(t+k|t)
    = F(z) / A(z) e(t)
    = F(z) / cancel(A(z)) cancel(A(z)) / C(z) y(t)
    = F(z) / C(z) y(t)
  $

  #definition(title: "Predictor from data")[
    Given a SSP $y(t) = W(z)e(t)$
    where
    - $W(z)$ in normal form
    - $e(t) ~ WN(0, lambda^2)$
    the _predictor from data_ is

    $ hat(y)(t+k|t) = F(z) / C(z) y(t) $

    where $F(z)$, $C(z)$ are from $W(z) = E(z) + (z^(-k) F(z)) / A(z)$.

    $hat(y)$ is an ARMA process (since $F(z)/C(z)$ is a rational transfer function applied to the stationary process $y(t)$).
  ]
]

#remark(title: "Remark: Predictors in practice")[
  How can we use the predictor from data in practice?
  Since our Stochastic process is stationary, we have that:
  $ hat(y)(t+k|t) = F(z) / C(z) y(t) arrow.r.double.long hat(y)(t|t-k) = F(z) / C(z) y(t-k) $
  This means we can get a predicition for the next possible value of our time-series. and we can compute it with data that we have, $y(t-k)$.
]

#remark(title: "Remark: worsening predictions")[
  The quality of the prediction gets worse with increasing $k$ since going forward will mean considering less samples $chevron.l y(0), dots, y(t) chevron.r arrow.r.long^(k arrow.t) chevron.l y(k), y(t) chevron.r$ since we don't have any data after $y(t)$
]
== Analysis of the prediction error

#remark(title: "Naïve (stupid) predictor")[
  The simplest predictor is the *naïve predictor*: $hat(y)(t+k|t) = y(t)$ (just repeat the last known value).

  This is generally outperformed by the optimal predictor, which exploits the model structure.
]

Let's consider the prediction error of the optimal predictor
$
  epsilon(t+k|t) & = y(t+k) - hat(y) (t+k|t) \
                 & = W(z)e(t+k) - hat(y) (t+k|t) \
                 & = (E(z) + F(z) / A(z) z^(-k)) e(t+k) - F(z) / A(z) e(t) \
                 & = E(z)e(t+k) + F(z) / A(z) z^(-k) e(t+k) - F(z) / A(z) e(t) \
                 & = E(z)e(t+k) + cancel(F(z) / A(z) e(t)) - cancel(F(z) / A(z) e(t)) \
                 & = E(z)e(t+k) \
                 & = cal(e)_o e(t) + cal(e)_1 e(t-1) + dots + cal(e)_(k-1) e(t-k+1)
$

So the variance of $epsilon(t|t+k)$ is
$
  EE[epsilon(t+k|k)^2] = (cal(e)_0^2 + cal(e)_1^2 + dots + cal(e)_(k-1)^2) lambda^2
$

#remark[
  If $k -> +infinity$ then $E(z)$ becomes the $MA(infinity)$ representation of $y(t)$.
  This means that $ EE[epsilon(t+k|t)^2] = EE[E(z)e(t+k)] -->_(k -> +infinity) EE[y(t)] $
]

#figure(
  cetz.canvas({
    import suiji: gen-rng
    import cetz.draw: *
    import cetz-plot: *
    import "../util.typ": random-series

    // variance function implementation
    let lambda2 = 1.0
    // exponentially decaying coefficients (typical of a stable causal system)
    let e(i) = if i == 0 { 0.35 } else { 1 / calc.pow(2, i * 0.45) }
    let var(k) = {
      let coefficients = array.range(0, k).map(e)
      return coefficients.map(e => calc.pow(e, 2)).sum(default: 0.0) * lambda2
    }
    let k-max = 10
    let points = array.range(1, k-max + 1).map(k => (k, var(k)))

    // estimate of the asymptotic value
    let asymptotic-value = var(50)

    // variance function math expression
    let var-math-expr(k) = {
      let coefficients = array.range(0, k + 1).map(k => $cal(e)_#k$).join($+$)
      $(#coefficients) lambda^2$
    }
    let var-ticks = array
      .range(0, 3)
      .map(k => {
        let (_x, y) = points.at(k)
        return (y, var-math-expr(k))
      })

    plot.plot(
      size: (10, 5),
      axis-style: "school-book",
      x-label: [$k$],
      x-min: 0.0,
      x-max: k-max,
      x-tick-step: 1.0,
      y-label: [$EE[epsilon(t+k|k)^2]$],
      y-min: 0.0,
      y-max: asymptotic-value * 1.1,
      y-tick-step: none,
      y-ticks: (..var-ticks, (asymptotic-value, $"var"[y(t)]$)),
      {
        plot.add(points, mark: "x", style: (stroke: (dash: "dashed", paint: blue)))
        plot.add-hline(asymptotic-value, style: (stroke: (dash: "dashed")))
      },
    )
  }),
)

Now lets consider the trivial predictor $hat(y)(t+k|t) = EE[y(t)] = m_y$

It's prediction error will be
$ epsilon(t+k|t) = y(t+k) - hat(y)(t+k) = y(t+k) - m_y $
and the variance
$ EE[epsilon(t+k|t)^2] = EE[(y(t+k) - m_y)^2] = "var"[y(t)] $

It follows that the optimal predictor always does at least as well as the trivial one:
- *Lower bound*: for $k=1$, the prediction error is $e(t+1)$, whose variance is $lambda^2 = "var"[e(t)]$.
- *Upper bound*: as $k -> infinity$, the predictor loses all information and the error approaches $"var"[y(t)]$ (the trivial predictor).
$ "var"[e(t)] <= "var"[epsilon(t+k|t)] < "var"[y(t)] $

== Initialization of the predictor

#example[
  Consider a predictor $ hat(y)(t+k|t) = 1 / (1 + c z^(-1)) y(t) $

  We can write
  $
    & (1 + c z^(-1)) hat(y)(t+k|t) = y(t) \
    & hat(y)(t+k|t) + c z^(-1) hat(y)(t+k|t) = y(t) \
    & hat(y)(t+k|t) + c hat(y)(t+k-1|t-1) = y(t) \
    & hat(y)(t+k|t) = - c hat(y)(t+k-1|t-1) + y(t)
  $

  Suppose we start collecting samples of $y(t)$ at $t=0$, then we would run into a problem:
  in order to calculate $hat(y)(k+0|0)$
  we would need $hat(y)(k-1|-1)$ which depends on $y(-1)$ which we do not have.
]

#problem[
  How to initialize the predictor? That is, what value do we assign to $hat(y)(k-1|-1)$ ?
]

#solution[
  We can assign any value, since $F(z) / C(z)$ si asymptotically stable, so any error will vanish rapidly.

  A common solution is to assign $hat(y)(k-1|-1) = m_y$.
]

== Prediction of processes with non-zero mean

#problem[
  We want to find the optimal predictor for a precess
  $ y(t) = C(z) / A(z) e(t), quad e(t) ~ WN(mu, lambda^2) $
  with $mu != 0$.
]

#remark[
  $y$ has mean value
  $ EE[y(t)] = m_y = W(1) mu = C(1) / A(1) mu $
]

#solution[
  We consider the de-biased processes
  $
    tilde(y)(t) & = y(t) - m_y = W(z) tilde(e)(t) \
    tilde(e)(t) & = e(t) - mu
  $
  and the de-biased predictor from data
  $
    hat(tilde(y))(t+k|t) = F(z) / C(z) tilde(y)(t)
  $

  If we flip the definition of $tilde(y)$ we get
  $ y(t) = tilde(y)(t) + m_y $
  which means that
  $ hat(y)(t+k|t) = hat(tilde(y))(t+k|t) + m_y $

  If we expand the expression, we get
  $
    hat(y)(t+k|t) & = hat(tilde(y))(t+k|t) + m_y \
                  & = F(z) / C(z) tilde(y)(t) + m_y \
                  & = F(z) / C(z) (y(t) - m_y) + m_y \
                  & = F(z) / C(z) y(t) - F(z) / C(z) m_y + m_y
  $

  Since $m_y$ is a constant, applying a transfer function $G(z)$ to a constant yields $G(1) dot m_y$ (the DC gain):
  $ F(z) / C(z) m_y = F(1) / C(1) m_y $

  And we can finally write
  $ hat(y)(t+k|t) = F(z) / C(z) y(t) + "bias" $
  where $"bias" = (1 - F(1) / C(1)) m_y$
]

== Prediction examples

#example(title: "1-step prediction for AR(1)")[
  Let $y(t) = a y(t-1) + e(t)$, $e(t) tilde WN(0, lambda^2)$, $|a| < 1$.

  Here $W(z) = 1/(1 - a z^(-1))$, so $C(z) = 1$, $A(z) = 1 - a z^(-1)$.

  For $k=1$: $C(z) = E(z) A(z) + z^(-1) F(z)$, with $E(z) = 1$, $F(z) = a$.

  *Predictor from noise:* $hat(y)(t+1|t) = a/(1-a z^(-1)) e(t) = a y(t) + a e(t) - a e(t) = a y(t)...$

  More directly: $hat(y)(t+1|t) = F(z)/C(z) y(t) = a y(t)$

  This makes intuitive sense: the best prediction for an AR(1) is simply $a$ times the current value.

  *Prediction error variance:* $"Var"[epsilon(t+1|t)] = lambda^2$ (just the noise variance).
]

#example(title: "1-step prediction for MA(1)")[
  Let $y(t) = e(t) + c e(t-1)$, $e(t) tilde WN(0, lambda^2)$.

  $W(z) = 1 + c z^(-1)$, $A(z) = 1$, $C(z) = 1 + c z^(-1)$.

  For $k=1$: $C(z) = E(z) dot 1 + z^(-1) F(z) arrow.r.double E(z) = 1, F(z) = c$.

  *Predictor from data:* $hat(y)(t+1|t) = c/(1+c z^(-1)) y(t)$

  *Prediction error:* $epsilon(t+1|t) = E(z) e(t+1) = e(t+1)$, so $"Var"[epsilon] = lambda^2$.
]

#caution-box(title: "Non-canonical form warning")[
  The prediction formulas assume the process is in *canonical form*. If the representation is not canonical (e.g., $C(z)$ has roots outside the unit circle), the predictor may be *unstable*.

  Always convert to canonical form before computing the predictor.
]

== Prediction for ARMAX models

#problem(title: "ARMAX prediction")[
  Given an ARMAX model:
  $ y(t) = B(z)/A(z) z^(-k_0) u(t) + C(z)/A(z) e(t) $

  Find the optimal $k$-step predictor.
]

#solution[
  Using the Diophantine equation $C(z) = E(z) A(z) + z^(-k) F(z)$ as before:
  $
    hat(y)(t+k|t) = F(z)/C(z) y(t) + (E(z) B(z))/C(z) z^(-k_0) u(t+k)
  $

  The predictor now has two components:
  + A *feedback* term from past outputs: $F(z)/C(z) y(t)$
  + A *feedforward* term from known inputs: $(E(z) B(z))/C(z) z^(-k_0) u(t+k)$

  Note: future values of $u(t)$ are assumed known (since $u$ is the controlled input).
]

#remark(title: "1-step ARMAX prediction")[
  For $k=1$: $E(z) = 1$, $F(z) = C(z) - A(z)$, and the predictor simplifies to:
  $ hat(y)(t+1|t) = (C(z) - A(z))/C(z) y(t) + B(z)/C(z) z^(-k_0) u(t+1) $
]
