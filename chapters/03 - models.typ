#import "/prelude.typ": *

= Models

== White noise

#definition(title: "White noise")[
  $e(t) ~ WN(mu, lambda^2)$ iff
  - $m_e = mu$
  - $gamma(0) = lambda^2$
  - $gamma(tau) = 0 quad forall tau != 0$

  White noises are SSPs.
]

#remark[
  The constant realization is admissible as a realization for a white noise but it is very low probability.
]

#properties(title: "spectral domain")[
  Let $e(t) ~ WN(mu, lambda^2)$

  Then $Gamma_e(omega) = lambda^2, forall omega in RR$
]

#figure(
  cetz.canvas({
    import suiji: gen-rng
    import cetz.draw: *
    import cetz-plot: *
    import "../util.typ": random-series
    let rng = gen-rng(1)
    let (rng, realization) = random-series(rng, 20, from: -1, to: 20)

    let gamma-func(tau, lambda2: 1.0) = {
      tau = calc.abs(tau) // gamma is even function

      return if tau == 0 {
        lambda2
      } else {
        0
      }
    }

    // inclusive range
    let range(from, to) = array.range(from, to + 1)

    let expected_value = range(-5, +5).map(x => (x, 0.6))
    let values = range(-5, +5).map(x => (x, gamma-func(x)))

    plot.plot(
      size: (10, 4),
      axis-style: "school-book",
      x-tick-step: 1,
      x-label: [$tau$],
      y-tick-step: none,
      y-label: none,
      {
        plot.add(label: [$ m_(e)(t) = mu quad forall t $], expected_value)
        plot.add(
          label: [$
            gamma_(e)(tau) = cases(
              lambda^2 quad & tau = 0,
              0 quad forall & tau != 0
            )
          $],
          mark: "x",
          style: (stroke: (dash: "dashed")),
          values,
        )
      },
    )
  }),
  caption: [Plot of $gamma(tau)$],
)

== Moving averages

#definition(title: [Moving Average of order $n$])[
  Let $e(t) ~ WN(mu, lambda^2)$

  Then $y(t) = c_0 e(t) + c_1 e(t-1) + dots + c_n e(t-n) ~ MA(n)$

  Moving averages are SSPs.
]

#properties[
  $m_y = mu dot sum_(i=0)^n c_i$
  $
    gamma(tau) = lambda^2 dot cases(
      c_0 c_0 + c_1 c_1 + dots + c_n c_(n) & quad tau = 0,
      c_0 c_1 + c_1 c_2 + dots + c_n c_(n-1) & quad tau = 1,
      c_0 c_2 + c_1 c_3 + dots + c_n c_(n-2) & quad tau = 2,
      dots.v,
      c_0 c_n & quad tau = n,
      0 & quad tau > n
    )
  $

  In general $ gamma(tau) = sum_(i=0)^(n-tau) c_i c_(i + tau) $
] <prop:moving-averages>

#figure(
  cetz.canvas({
    import cetz.draw: *
    import cetz-plot: *

    // MA(2) with c0=1, c1=0.8, c2=0.3, lambda^2=1
    let c = (1.0, 0.8, 0.3)
    let lambda2 = 1.0
    let n = 2

    // Compute gamma(tau) = lambda^2 * sum_{i=0}^{n-tau} c_i * c_{i+tau}
    let gamma-func(tau) = {
      tau = calc.abs(tau)
      if tau > n { return 0.0 }
      let s = 0.0
      for i in range(0, n - tau + 1) {
        s += c.at(i) * c.at(i + tau)
      }
      return lambda2 * s
    }

    let range(from, to) = array.range(from, to + 1)
    let values = range(-5, +5).map(x => (x, gamma-func(x)))

    plot.plot(
      size: (10, 3.5),
      axis-style: "school-book",
      x-tick-step: 1,
      x-label: [$tau$],
      y-tick-step: none,
      y-ticks: ((gamma-func(0), $gamma(0)$), (gamma-func(1), $gamma(1)$), (gamma-func(2), $gamma(2)$)),
      y-label: [$gamma(tau)$],
      {
        plot.add(
          values,
          mark: "x",
          style: (stroke: (dash: "dashed")),
        )
      },
    )
  }),
  caption: [Covariance function of an $MA(2)$ process ($c_0 = 1, c_1 = 0.8, c_2 = 0.3$). Note that $gamma(tau) = 0$ for $|tau| > n = 2$: the memory of the process is *finite*.],
)

#note-box[
  - $gamma (tau)$ has $n - tau + 1$ non-zero terms.
  - Moving averages are called *colored processes*
  - For $forall tau > n$ the model acts as white noise
]

#caution-box(title: "Parameter redundancy and monic convention")[
  The $MA(n)$ process has $n+2$ parameters: $c_0, c_1, dots, c_n, lambda^2$.

  This representation is *redundant*: the process with parameters
  $ tilde(c)_i = alpha c_i, quad tilde(lambda)^2 = lambda^2 / alpha^2 $
  has identical mean and covariance function.

  To avoid this, $c_0$ is set to 1 ($C(z)$ is a *monic* polynomial).
  The process then has $n+1$ free parameters: $c_1, dots, c_n, lambda^2$.
]

=== Generalized MA (GMA) process

#definition(title: "Generalized MA process")[
  Consider
  $ v(t) = c_0 eta(t) + c_1 eta(t-1) + dots + c_n eta(t-n) $
  where $eta(dot)$ is stationary but *not necessarily white*.

  The cross-terms in the variance calculation can no longer be neglected:
  $ "Var"[v(t)] = sum_(i=0)^n sum_(j=0)^n c_i c_j gamma_(eta eta)(i-j) $

  The GMA process is still stationary --- passing a stationary process through an MA filter preserves stationarity.
]

#properties(title: "zeros and poles of MA model's transfer function")[
  All $MA(n)$ processes have
  - $n$ non-trivial zeros
  - $n$ poles, all lying at the origin

  For this reason $AR(n)$ processes are called "all-zeros" processes.
]

=== Moving averages of order $infinity$
With the objective of modelling every kind of stochastic process, including those that have memory of the previous inputs and previous outputs, we define processes which have memory of any input up to that point because any previous output will also be a linear combination of previous inputs.

#definition(title: [Moving Average of order $infinity$])[
  Let $e(t) ~ WN(mu, lambda^2)$

  Then, under the assumption that $sum_(i=0)^infinity c_i < infinity$
  $ y(t) = c_0 e(t) + c_1 e(t-1) + dots = sum_(i=0)^infinity c_i e(t) ~ MA(infinity) $

  Moving averages are SSPs.
]

#properties[
  A generalized version of @prop:moving-averages also holds for $MA(infinity)$
  - $m_y = mu dot sum_(i=0)^infinity c_i$

  - $gamma(tau) = sum_(i=0)^infinity c_i c_(i + tau)$
]

== Autoregressive models

#definition(title: "Autoregressive model")[
  Let $e(t) ~ WN(mu, lambda^2)$

  Then $y(t) = e(t) + a_1 y(t-1) + a_2 y(t-2) + dots + a_n y(t-m) ~ AR(m)$
]

#remark[
  The AR process is the *steady-state solution* of the difference equation
  $ v(t) = a_1 v(t-1) + dots + a_m v(t-m) + e(t) $
  If $W(z) = 1/A(z)$ is asymptotically stable (all poles of $A(z)$ strictly inside the unit circle), then $v(t)$ converges to a unique SSP regardless of the initial condition.
]

#remark[
  With autoregressive models we can have $gamma(tau) != 0$ for $tau -> infinity$ with a finite number of coefficients.
  - This is not possible with $MA(n)$ models.
  - This is possible with $MA(infinity)$ models, but working with an infinite number of coefficients is cumbersome.

  Any $AR(m)$ model with stable poles can be expressed as an $MA(infinity)$ model.
]

#definition(title: "Steady state solution for an AR model")[
  Let $y(t) ~ AR(m)$ model, we can expand its definition

  $
    y(t) & = a y (t-1) + e(t) \
         & = a^2 y(t-2) + a dot e(t-1)+ e(t) \
         & = a^2 y(t-2) + a dot e(t-1)+ e(t) \
         & = dots \
         & = a^0e(t) + a^1 e(t-1) + a^2 e(t-2) + dots +a^m e(t-m) + dots + a^{t-t_0}y(t_0)
  $
  The term $a^{t-t_0}y(t_0)$ is called *initial condition*, the rest is the *steady state solution*
]

#properties(title: "zeros and poles of AR model's transfer function")[
  All $AR(m)$ processes have
  - $m$ zeros, all lying at the origin
  - $m$ non-trivial poles

  For this reason $AR(m)$ processes are called "all-poles" processes.

  The poles of $A(z)$ are found by solving $1 + a_1 z^(-1) + dots + a_m z^(-m) = 0$ (or equivalently $z^m + a_1 z^(m-1) + dots + a_m = 0$).
]

#caution-box[
  Considering the initial condition in the expression of the process, namely the combined $a^{t-t_0}y(t_0)$ can make the SP a non stationary one in the initial transient.
]

Hence, to find the mean and covariance of an $AR(n)$ model we must first check that it is a SSP. //In the case we're uncertain about the model type, we check the stability or ensure it's a $MA(n)$ process.

We can do that by calculating its transfer function, finding its poles and applying @thm:stationarity.

Then we can calculate the mean value $m$ and covariance function $gamma(tau)$ by taking advantage of its recursive nature.

#example[
  Let $e(t) ~ WN(0, lambda^2)$ and consider
  $ y(t) = a y(t-1) + e(t) ~ AR(1) $

  Calculating the mean value $m_y$
  $
    m_y & = EE[y(t)] \
        & = EE[a y(t-1) + e(t)] \
        & = a EE[y(t-1)] + EE[e(t)] \
        & = a EE[y(t)] + EE[e(t)] \
        & = a m_y + 0 \
  $
  $
    =>
    m_y & = a m_y \
    m_y & = 0
  $

  Calculating variance $gamma(0)$
  $
    gamma(0) & = EE[(y(t) - m_y)(y(t-0) - m_y)] \
             & = EE[(y(t))^2] \
             & = EE[(a y(t-1) + e(t))^2] \
             & = EE[a^2 y(t-1)^2 + e(t)^2 + 2 a y(t-1) e(t)] \
             & = a^2 EE[y(t-1)^2] + EE[e(t)^2] + 2 a underbrace(cancel(EE[y(t-1) e(t)]), #[see @thm:null-expected-value]) \
             & = a^2 EE[y(t-1)^2] + EE[e(t)^2] \
             & = a^2 EE[y(t)^2] + EE[e(t)^2] \
             & = gamma_y(0) + lambda^2
  $
  $
    => gamma_y(0) & = a^2 gamma_y(0) + lambda^2 \
       gamma_y(0) & = lambda^2 / (1 - a^2)
  $

  Calculating covariance $gamma(tau)$
  $
    gamma(1) &= EE[(y(t) - m_y)(y(t-1) - m_y)] \
    &= EE[y(t)y(t-1)] \
    &= EE[(a y(t-1) + e(t)) y(t-1)] \
    &= EE[a y(t-1)^2 + e(t)y(t-1)] \
    &= a EE[y(t-1)^2] + underbrace(cancel(EE[e(t)y(t-1)]), #[Noise is in the future for \ the considered outcome]) \
    &= a gamma(0)
  $
  $
    gamma(2) & = EE[(y(t) - m_y)(y(t-2) - m_y)] \
             & = EE[y(t)y(t-2)] \
             & = EE[(a y(t-1) + e(t)) y(t-2)] \
             & = a EE[t(t-1)y(t-2)] + cancel(EE[e(t)y(t-2)]) \
             & = a gamma(1)
  $

  We get a set of recursive equations:
  $
    cases(
      gamma(0) = 1 / (1-a^2) lambda^2,
      gamma(tau) = a gamma(tau - 1) quad |tau| >= 1,
    )
  $

  #figure(
    cetz.canvas({
      import suiji: gen-rng
      import cetz.draw: *
      import cetz-plot: *
      import "../util.typ": random-series
      let rng = gen-rng(1)
      let (rng, realization) = random-series(rng, 20, from: -1, to: 20)

      let gamma-func(tau, a: 0.5, lambda2: 1.0) = {
        tau = calc.abs(tau) // gamma is even function

        return if tau == 0 {
          1 / (1 - calc.pow(a, 2)) * lambda2
        } else {
          a * gamma-func(tau - 1, a: a, lambda2: lambda2)
        }
      }

      // inclusive range
      let range(from, to) = array.range(from, to + 1)

      let values-pos-a = range(-5, +5).map(x => (x, gamma-func(x, a: +0.5)))
      let values-neg-a = range(-5, +5).map(x => (x, gamma-func(x, a: -0.5)))

      plot.plot(
        size: (10, 4),
        axis-style: "school-book",
        x-tick-step: 1,
        x-label: [$tau$],
        y-tick-step: none,
        y-label: none,
        {
          plot.add(
            label: [$gamma(tau), |a|<1, a>0$],
            values-pos-a,
            mark: "x",
            line: "spline",
            style: (stroke: (dash: "dashed")),
          )
          plot.add(
            label: [$gamma(tau), |a|<1, a<0$],
            values-neg-a,
            mark: "x",
            line: "spline",
            style: (stroke: (dash: "dashed")),
          )
        },
      )
    }),
    caption: [Plot of $gamma(tau)$],
  )
]

#remark(title: "AR(1) realization smoothness")[
  The sign of $a$ determines the character of AR(1) realizations:
  - $a > 0$: $gamma(tau) = a^(|tau|) gamma(0)$ is always positive and decays monotonically $arrow.r$ *smooth* realizations (consecutive values tend to be similar)
  - $a < 0$: $gamma(tau)$ alternates in sign $arrow.r$ *nervous* realizations (consecutive values tend to oscillate)
]

=== Yule-Walker equations

#definition(title: "Yule-Walker equations")[
  For an $AR(m)$ process, the covariance function satisfies the *Yule-Walker equations*:
  $ gamma(tau) = a_1 gamma(tau - 1) + a_2 gamma(tau - 2) + dots + a_m gamma(tau - m) quad "for" |tau| >= 1 $

  These are obtained by multiplying the AR equation by $v(t - tau)$ and taking expectations.
]

#example(title: "Yule-Walker for AR(1)")[
  $
    cases(
      gamma(0) = a gamma(1) + lambda^2 quad & ("from" tau=0)","
                                              gamma(tau) = a gamma(tau - 1) quad & "for" |tau| >= 1
    )
  $
  Solving: $gamma(0) = lambda^2 / (1-a^2)$, $space gamma(tau) = a^(|tau|) gamma(0)$.
]

#example(title: "Yule-Walker for AR(2)")[
  $
    cases(
      gamma(0) = a_1 gamma(1) + a_2 gamma(2) + lambda^2,
      gamma(1) = a_1 gamma(0) + a_2 gamma(1),
      gamma(tau) = a_1 gamma(tau - 1) + a_2 gamma(tau - 2) quad "for" |tau| >= 2
    )
  $
  From the second equation: $gamma(1) = a_1 / (1-a_2) gamma(0)$.
]

=== AR(1) to MA(infinity) equivalence

#theorem(title: "AR(1) as MA(infinity) via long division")[
  Consider $v(t) = a v(t-1) + e(t)$, i.e., $v(t) = 1/(1-a z^(-1)) e(t)$.

  Performing polynomial long division of $1$ by $(1 - a z^(-1))$:
  $ 1/(1 - a z^(-1)) = 1 + a z^(-1) + a^2 z^(-2) + dots $

  Therefore: $v(t) = e(t) + a e(t-1) + a^2 e(t-2) + dots = sum_(k=0)^infinity a^k e(t-k) tilde MA(infinity)$

  This converges iff $|a| < 1$ (the pole is inside the unit circle).
]

#theorem[
  Let $e(t) ~ WN(0, lambda^2)$ and $y(t) = e(t) + a y(t-1) ~ AR(1)$

  Then $EE[y(t-tau) e(t)] = 0, forall tau > 0$
] <thm:null-expected-value>

#remark[
  @thm:null-expected-value is also valid for $AR(n)$ (trust me bro).
]

#proof[
  Remember that $y(t) = e(t) + a e(t-1) + a^2 e(t-2) + dots$

  Then $y(t-tau) = e(t-tau) + a e(t-tau-1) + a^2 e(t-tau-2) + dots$

  $
    EE[e(t) y(t-tau)] & = EE[e(t) (e(t-tau) + a e(t-tau-1) + a^2 e(t-tau-2) + dots)] \
                      & = EE[e(t)e(t-tau)] + a EE[e(t)e(t-tau-1)] + + a EE[e(t)e(t-tau-2)] + dots \
                      & = 0 + 0 + 0 + dots
  $
]

== ARMA models
The modeling power of $MA(infinity)$ is unmatched by $AR(1)$ models but we can combine an $AR(m)$ and a $MA(n)$ models to get an $ARMA(m, n)$ model to model our time series with the least amount of coefficients possible.

=== ARMA to MA(infinity) equivalence

#theorem(title: "ARMA as MA(infinity)")[
  Any $ARMA(m, n)$ process with $W(z) = C(z)/A(z)$ asymptotically stable can be written as:
  $ v(t) = W(z) e(t) = sum_(k=0)^infinity w_k e(t-k) tilde MA(infinity) $
  where $w_k$ are the coefficients of the impulse response (obtained by long division of $C(z)$ by $A(z)$).

  The series converges iff all poles of $A(z)$ are strictly inside the unit circle.
]

=== ARMA as series connection

#remark[
  An $ARMA(m, n)$ can be seen as a *series connection* of an AR part and an MA part:
  $ v(t) = C(z)/A(z) e(t) = underbrace(C(z), "MA") dot underbrace(1/A(z), "AR") e(t) $
  The order of the series connection does not matter (LTI systems commute).
]

#definition(title: "ARMA model")[
  Let $e(t) ~ WN(mu, lambda^2)$

  Then $y(t) ~ ARMA(m, n)$ iff
  $
    y(t) = & underbrace(a_1 y(t-1) + a_2 y(t-2) + dots + a_m y(t-m), AR(m)) + \
           & underbrace(c_0 e(t) + c_1 e(t-1) + dots + c_n e(t-n), MA(n)) \
  $
]

#example(title: "ARMA(1,1) covariance computation")[
  Let $v(t) = a v(t-1) + c_0 e(t) + c_1 e(t-1)$, $e(t) tilde WN(0, lambda^2)$, with $c_0 = 1$.

  *Mean value:* $EE[v(t)] = a EE[v(t)] + 0 arrow.r.double m_v = 0$

  *Variance* ($tau = 0$): Multiply defining equation by $v(t)$ and take expectations:
  $ gamma(0) = a gamma(1) + lambda^2(1 + c_1^2 + 2 a c_1) $

  *Covariance* ($tau = 1$): $gamma(1) = a gamma(0) + c_1 lambda^2$

  *Covariance* ($tau >= 2$): $gamma(tau) = a gamma(tau - 1)$ (same recursive structure as AR)

  Solving the system:
  $ gamma(0) = ((1 + c_1^2)(1-a^2) + 2 a c_1) / (1 - a^2) lambda^2, quad gamma(1) = a gamma(0) + c_1 lambda^2 $
]

#figure(
  cetz.canvas({
    import cetz.draw: *
    import cetz-plot: *

    // ARMA(1,1) with a=0.7, c1=0.5, lambda^2=1
    let a = 0.7
    let c1 = 0.5
    let lambda2 = 1.0

    // gamma(0) = ((1 + c1^2)(1 - a^2) + 2*a*c1) / (1 - a^2) * lambda^2
    let g0 = ((1 + c1 * c1) * (1 - a * a) + 2 * a * c1) / (1 - a * a) * lambda2
    // gamma(1) = a * gamma(0) + c1 * lambda^2
    let g1 = a * g0 + c1 * lambda2
    // gamma(tau) = a * gamma(tau-1) for tau >= 2

    let gamma-func(tau) = {
      tau = calc.abs(tau)
      if tau == 0 { return g0 }
      if tau == 1 { return g1 }
      // recursive: gamma(tau) = a^(tau-1) * gamma(1)
      return calc.pow(a, tau - 1) * g1
    }

    let range(from, to) = array.range(from, to + 1)
    let values = range(-8, +8).map(x => (x, gamma-func(x)))

    plot.plot(
      size: (10, 3.5),
      axis-style: "school-book",
      x-tick-step: 1,
      x-label: [$tau$],
      y-tick-step: none,
      y-ticks: ((g0, $gamma(0)$), (g1, $gamma(1)$)),
      y-label: [$gamma(tau)$],
      {
        plot.add(
          values,
          mark: "x",
          line: "spline",
          style: (stroke: (dash: "dashed")),
        )
      },
    )
  }),
  caption: [Covariance function of an $ARMA(1, 1)$ process ($a = 0.7, c_1 = 0.5$). For $|tau| >= 2$, $gamma(tau) = a dot gamma(tau-1)$: the tail decays exponentially like an AR(1), but the initial values $gamma(0), gamma(1)$ are shaped by the MA part.],
)

=== Other variants of ARMA

If we also consider the input of the system (the eXogenous part) we get an $ARMAX$ model.

#definition(title: "ARMAX model")[
  Let $e(t) ~ WN(mu, lambda^2)$

  Then $y(t) ~ ARMAX(m, n, k, p)$ iff
  $
    y(t) = & underbrace(a_1 y(t-1) + a_2 y(t-2) + dots + a_m y(t-m), AR(m)) + \
           & underbrace(c_0 e(t) + c_1 e(t-1) + dots + c_n e(t-n), MA(n)) + \
           & underbrace(b_0 u(t-k) + b_1 u(t-k-1) + dots + b_n u(t-k-p), "X"(k, p))
  $

  Where
  - $k$: pure input/output delay
  - $p$: order of the exogenous part
]

#remark(title: "ARMAX block diagram")[
  The ARMAX model can be represented as:
  #{
    import fletcher: diagram, edge, node
    figure(
      diagram(
        node-shape: "rect",
        node-stroke: 1pt,
        spacing: (3em, 2em),
        edge((-1, 0), "r", "-|>")[$u(t)$],
        node((0, 0))[$B(z)/A(z) z^(-k)$],
        edge("r", "-|>"),
        node((1, 0), shape: "circle", inset: 2pt, sym.plus),
        edge((0, -1), "r,d", "-|>"),
        node((0, -1))[$C(z)/A(z)$],
        edge((-1, -1), "r", "-|>")[$e(t)$],
        edge((1, 0), "r", "-|>")[$y(t)$],
      ),
    )
  }
  The output is the sum of the filtered input (deterministic) and filtered noise (stochastic).
]

#remark[
  The input/output delay $k$ of an $ARMAX(m, n, k, p)$ process is visible in the step response graph of the eXogenous part.
  #figure(
    cetz.canvas({
      import cetz.draw: *
      import cetz-plot: *
      plot.plot(
        size: (6, 4),
        axis-style: "school-book",
        x-min: -1,
        x-max: +5,
        x-tick-step: none,
        x-label: $t$,
        y-min: -2,
        y-max: +2,
        y-tick-step: none,
        {
          plot.add(label: "step", ((-1, 0), (1, 0), (1, 1), (5, 1)))

          plot.add(
            label: "step response",
            domain: (-1, 5),
            sample-at: (2,), // add a sample to keep the line sharp
            x => if x < 2 {
              0
            } else {
              1 - calc.exp(-1.75 * (x - 2))
            },
          )

          plot.annotate({
            cetz.decorations.flat-brace(name: "k-brace", (1, 0), (2, 0), flip: true)
            content("k-brace.content")[$k$]
          })
        },
      )
    }),
  )
]

If we apply a non-linear function to an $ARMAX$ model we get a $"N-ARMAX"$ model.

#definition(title: "N-ARMAX model")[
  Let $e(t) ~ WN(mu, lambda^2)$

  Let $f$ be any non-linear function (polynomials, splines, wavelengths, neural networks etc.)

  Then $y(t) ~ "N-ARMAX"(m, n, k, p)$ iff
  $
    y(t) = f( quad & a_1 y(t-1), a_2 y(t-2), dots, a_m y(t-m), \
                   & c_0 e(t), c_1 e(t-1), dots, c_n e(t-n), \
                   & b_0 u(t-k), b_1 u(t-k-1), dots, b_n u(t-k-p))
                     quad )
  $
]

=== Unbiased models
When calculating the covariance function of a process, it's often useful to define the unbiased version of the process.

#definition(title: "Unbiased model")[
  Let $y(t) ~ ARMA(m, n)$, $e(t) ~ WN(m_e, lambda^2)$

  $tilde(e)(t) = e(t) - m_e$

  $tilde(y)(t) = y(t) - m_y$ is the unbiased version of the process
]

We can prove that the covariance function of the unbiased version of a process is the same as the covariance function of the original process.

#theorem(title: "Covariance function of unbiased model")[
  $EE[e(t)e(t- tau)] = 0$ if $EE[e(t)] = 0$

  $EE[ underbrace((e(t) - m_e), tilde(e)(t)) underbrace(e(t- tau) - EE[e(t-tau)], tilde(e)(t- tau))] = 0 forall tau != 0$

  Therefore

  $gamma_y (tau) &= EE[y(t) - m_y)(y(t-tau) - E[y(t-tau))] \
  &= EE[tilde(y)(t) tilde(y)(t-tau)] \ &= gamma_tilde(y)(tau)$
]

#theorem(title: "of the gain")[

  Using the unbiased version of the process we can get an equivalent representation
  #{
    import fletcher: diagram, edge, node
    figure(
      diagram(
        node-shape: "rect",
        node-stroke: 1pt,
        edge((0, 0), "r", "-|>")[$tilde(e)(t)$],
        edge((1, 1), "u", "-|>")[$m_e$],
        edge((1, 0), "r", "-|>")[$e(t)$],

        node((2, 0))[$G(z)$],
        edge((2, 0), "rr", "-|>")[$y(t)$],
        node((1, 0), shape: "circle", inset: 2pt, sym.plus.minus),
      ),
    )
  }
  #{
    import fletcher: diagram, edge, node
    figure(
      diagram(
        node-shape: "rect",
        node-stroke: 1pt,
        edge((0, 0), "rr", "-|>")[$tilde(e)(t)$],

        node((2, 0))[$G(z)$],
        edge((2, 0), "rrr", "-|>")[$tilde(y)(t)$],
        edge((0, 1), "rr", "-|>")[$u= m_e$],
        edge((2, 1), "rrr", "-")[$U = m_y = EE[y(t)]$],
        edge((5, 1), "u", "-|>"),

        node((2, 1))[$G(z)$],

        node((5, 0), shape: "circle", inset: 2pt, sym.plus.minus),
        edge((5, 0), "r", "-|>")[$y(t)$],
      ),
    )
  }

  We can use this to compute the expected value of our stochastic process

  At *steady state* we have that

  $EE[y(t)] = lim_(t->1) = G(z) dot m_e = G(1) dot U$

]
