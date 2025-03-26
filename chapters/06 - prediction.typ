#import "/prelude.typ": *

= Prediction

#problem(title: "Prediction problem")[
  Given
  - a SSP $y(t)$ in canonical form
    $ y(t) = C(z) / A(z) e(t), quad e(t) ~ WN(mu, lambda^2) $
    with $C(z)$, $A(z)$, $mu$, $lambda^2$ known.

  - a specific realization $y(0), y(1), dots y(t)$

  How can we predict future values of the process at a generic time instant $t+k$ ?

  // FIXME: confusing variable naming
  // - t is both the variable in y(t) and also the time index in t+k
  // - in the plot t is both the x axis label and a mark in the x axis

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
  $hat(y) (t+k|t)$ is the predictor of SSP $y(t)$.

  // TODO: explain
]

#definition(title: "Trivial predictor")[
  $hat(y) (t+k|t) = EE[y(t)] = m_y$
]

#definition(title: "Linear predictor")[
  An "optimal" linear predictor contains both information about the model information $angle.l alpha_0, alpha_1, dots angle.r$ and the data $ angle.l y(t), y(t-1), dots angle.r$.
  $
    hat(y)(t+k|t) &= alpha_0 y(t) + alpha_1 y(t - 1) + dots \
    &= sum_(i=0)^(+infinity) (alpha_i) y(t-i)
  $

  with $alpha_i in RR$ such that $sum_(i=0)^(+infinity) alpha_i < +infinity$

  // TODO: note that alpha_i is info from the model and y(t-i) is info from the data
]

#definition(title: "Prediction error")[
  The prediction error of $hat(y)$ is
  $ epsilon(t+k|t) = y(t+k) - hat(y) (t+k|t) $
] <def:prediction-error>

#definition(title: "Mean squared prediction error")[
  The MSP of $hat(y)$ is
  $
    cal(J) (hat(y)) &= EE[epsilon(t+k|t)^2] \
    &= EE[(y(t+k) - hat(y) (t+k|t))^2]
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
    y(t) &= W(z)e(t) \ &= w_0 e(t) + w_1 e(t-1) + w_2 e(t-2) + dots \
    &= sum_(i=0)^(+infinity) w_i e(t-i) \
    y(t-1) &= sum_(i=0)^(+infinity) w_i e(t-i-1) \
    y(t-2) &= sum_(i=0)^(+infinity) w_i e(t-i-2) \
    &dots.v
  $

  So we can rewrite our linear predictor in the following way
  $
    hat(y)(t+k|t) &= alpha_0 y(t) + alpha_1 y(t - 1) + dots \
    &= alpha_0 sum_(i=0)^(+infinity) w_i e(t-i) + alpha_1 sum_(i=0)^(+infinity) w_i e(t-i-1) + dots \
    &= sum_(i=0)^(+infinity) beta_i e(t-i)
  $

  // TODO: what is the definition of beta_i?
  // how exactly do we re-arrange the w_i and alpha_i to get beta_i?

  Then, we reformulate the optimization problem wrt $alpha$ in a new optimization problem wrt $beta$
  $ cal(V) = EE[(y(t+k) - hat(y)(t+k|t, beta))^2] $
  we are now looking for $hat(beta) = argmin_beta cal(V) (beta)$

  Let's consider $y(t+k)$ alone, we can split it into two parts
  $
    y(t+k) &= sum_(i=0)^(+infinity) w_i e(t-i-k) \
    &=
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
  // FIXME: this formula is monstrous
  // FIXME: it goes out of the page margins
  $
    cal(V)(beta) &= EE[(y(t+k) - hat(y)(t+k|t, beta))^2] \
    &= EE[( sum_(j=0)^(k-1) w_j e(t+k-j) +
        sum_(i=0)^(+infinity) w_(k+i) e(t-i) -
        sum_(i=0)^(+infinity) beta_i e(t-i) )^2] \
    &= EE[ underbrace((sum_(j=0)^(k-1) w_j e(t+k-j))^2, "not a function of" beta) +
      (sum_(i=0)^(+infinity) (w_(k+i) - beta_i) e(t-i))^2 + underbrace(
        cancel(
          2 (sum_(j=0)^(k-1) w_j e(t+k-j)) (sum_(i=0)^(+infinity) (w_(k+i) - beta_i) e(t-i))
        ),
        "all uncorrelated because WN at different time instants"
        ) ]
  $

  So the optimal predictor is given by the minimizer of
  $ EE[(sum_(i=0)^(+infinity) (w_(k+i) - beta_i) e(t-i))^2] >= 0 $
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

  // TODO: move this somewhere more appropriate
  #definition(title: "Diophantine equation")[
    $ C(z) = E(z)A(z) + z^(-k) F(z) $
  ]

  If we substitute $W(z) = E(z) + (z^(-k) F(z)) / A(z)$ into $y(t+k)$ we can observe that
  $
    y(t+k) &= W(z) e(t+k) \
    &= (E(z) + (z^(-k) F(z)) / A(z)) e(t+k) \
    &= E(z) e(t+k) + F(z) / A(z) z^(-k) e(t+k) \
    &=
    underbrace(E(z) e(t+k), "(a)") +
    underbrace(F(z) / A(z) e(t), "(b)")
  $

  // TODO: actual links
  #[
    #set enum(numbering: "(a)")
    + depends on future samples of $e(t)$ and is unpredictable at time $t$
    + depends on past samples of $e(t)$ and is computable at time $t$ and will be our predictor
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
  the _whitening filter_ is
  // TODO: i called it X because we did not give it a name
  $ X(z) = A(z) / C(z) $

  The whitening filter "de-colors" $y(t)$ and turns it into white noise
  #{
    import fletcher: diagram, node, edge
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

    $hat(y)$ is a SSP. // TODO: is it ARMA process?
  ]
]

#remark[
  How can we use the predictor from data in practice?
  Since our Stochastic process is stationary, we have that: 
    $ hat(y)(t+k|t) = F(z) / C(z) y(t) arrow.r.double.long  hat(y)(t|t-k) = F(z) / C(z) y(t-k) $ 
  This means we can get a predicition for the next possible value of our time-series. and we can compute it with data that we have, $y(t-k)$.
]

#remark[
  The quality of the prediction gets worse with increasing $k$ since going forward will mean considering less samples $angle.l y(0), dots, y(t) angle.r arrow.r.long^(k arrow.t) angle.l y(k), y(t) angle.r$ since we don't have any data after $y(t)$
]
== Analysis of the prediction error

Let's consider the prediction error of the optimal predictor
$
  epsilon(t+k|t) &= y(t+k) - hat(y) (t+k|t) \
  &= W(z)e(t+k) - hat(y) (t+k|t) \
  &= (E(z) + F(z) / A(z) z^(-k)) e(t+k) - F(z) / A(z) e(t) \
  &= E(z)e(t+k) + F(z) / A(z) z^(-k) e(t+k) - F(z) / A(z) e(t) \
  &= E(z)e(t+k) + cancel(F(z) / A(z) e(t)) - cancel(F(z) / A(z) e(t)) \
  &= E(z)e(t+k) \
  &= cal(e)_o e(t) + cal(e)_1 e(t-1) + dots + cal(e)_(k-1) e(t-k+1)
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
    // TODO: check if this sequence of e_i coefficients is realistic
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

It follows that // TODO: explain why it follows
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
    tilde(y)(t) &= y(t) - m_y = W(z) tilde(e)(t) \
    tilde(e)(t) &= e(t) - mu
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
    hat(y)(t+k|t) &= hat(tilde(y))(t+k|t) + m_y \
    &= F(z) / C(z) tilde(y)(t) + m_y \
    &= F(z) / C(z) (y(t) - m_y) + m_y \
    &= F(z) / C(z) y(t) - F(z) / C(z) m_y + m_y
  $

  Since $m_y$ is a constant // TODO: cite the theorem
  $ F(z) / C(z) m_y = F(1) / C(1) m_y $

  And we can finally write
  $ hat(y)(t+k|t) = F(z) / C(z) y(t) + "bias" $
  where $"bias" = (1 - F(1) / C(1)) m_y$
]
