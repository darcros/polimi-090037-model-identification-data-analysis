#import "/prelude.typ": *

= Stochastic processes

Stochastic processes extend the notion of random variable to signals.

#definition(title: "Stochstic process")[
  A stochastic process is an infinite sequence of random variables, all defined on the same probabilistic space

  $ v(1, s), v(2, s), dots, v(t, s) $

  where $t$ is the time index and $s in S$ denotes the outcome in the sample space $S$.
]

#figure(
  cetz.canvas({
    import suiji: gen-rng
    import cetz.draw: *
    import cetz-plot: *
    import "../util.typ": random-series
    let rng = gen-rng(1)
    let (rng, realization) = random-series(rng, 20, from: -1, to: 20)

    plot.plot(
      size: (10, 2),
      axis-style: "school-book",
      x-min: -1,
      x-max: +20,
      x-tick-step: none,
      x-label: [$t$],
      y-tick-step: none,
      y-label: [$v(t, overline(s))$],
      {
        plot.add(label: "realization", realization, line: "spline")
      },
    )
  }),
  caption: [If $s$ is fixed we get a deterministic function],
)

#figure(
  cetz.canvas({
    import suiji: gen-rng, normal
    import cetz.draw: *
    import cetz-plot: *
    import "../util.typ": random-series

    let rng = gen-rng(1)
    let (rng, v) = normal(rng, size: 15)
    let pts = v.map(value => (5.0, value))

    let gaussian(x, mu: 0.0, sigma: 1.0) = {
      let a = 1 / (sigma * calc.sqrt(2 * calc.pi))

      let sigma2 = calc.pow(sigma, 2)
      let exponent = -calc.pow(x - mu, 2) / (2 * sigma2)

      return a * calc.exp(exponent)
    }

    plot.plot(
      size: (10, 2),
      axis-style: "school-book",
      x-min: -1,
      x-max: +20,
      x-tick-step: none,
      x-label: [$t$],
      y-tick-step: none,
      y-label: [$v(overline(t), s)$],
      {
        plot.add(label: "realizations", pts, mark: "+", style: (stroke: none))
        plot.add(label: "distribution", domain: (-3, +3), x => 5.0 + 5 * gaussian(x), axes: ("y", "x"))
      },
    )
  }),
  caption: [If $t$ is fixed we get a random variable],
)

#definition(title: "Stochstic equivalence")[
  Two signals $v_1(t), v_2(t)$ are stochastically equivalent if they are realizations of the same SP
]

== Wide-Sense characterization

To fully describe a SP we would need to know how the probability distribution of $v(s)$ evolves over time.
This is too heavy for modeling purposes.

For this reason we use the wide-sense characterization, which means that we describe a SP with two functions:
- the mean value $m(t)$
- the covariance function $gamma(t_1, t_2)$

#definition(title: "Wide-Sense equivalence")[
  Two SPs $v_1$ and $v_2$ are wide-sense equivalent iff
  - $m_v_1 = m_v_2$ their mean value is the same
  - $gamma_v_1 = gamma_v_2$ their covariance function is the same
]

#definition(title: "Mean value of a SP")[
  $m(t) &= EE[v(t, s)] = integral_PP v(t,s) "pdf"(s) \ds$
]

#figure(
  cetz.canvas({
    import suiji: gen-rng
    import cetz.draw: *
    import cetz-plot: *
    import "../util.typ": avg-series, random-series
    let rng = gen-rng(1)

    let realizations = ()
    let i = 0
    while (i < 5) {
      let (r, realization) = random-series(rng, 20, from: -1, to: 20)
      realizations.push(realization)

      rng = r
      i += 1
    }

    let mean = avg-series(realizations)

    plot.plot(
      size: (10, 2),
      axis-style: "school-book",
      x-min: -1,
      x-max: +20,
      x-tick-step: none,
      x-label: [$t$],
      y-tick-step: none,
      y-label: none,
      {
        for r in realizations {
          plot.add(r, line: "spline", style: (stroke: blue.lighten(25%)))
        }
        plot.add-legend("realizations")

        plot.add(label: "mean", mean, line: "spline", style: (stroke: (paint: red, thickness: 2pt)))
      },
    )
  }),
  caption: [Mean value and a bunch of realizations],
)

#definition(title: "Covariance function of a SP")[
  $gamma(t_1, t_2) = EE[(v(t_1)-m(t_1))(v(t_2)-m(t_2))]$

  The #underline[variance function] $gamma(t) = gamma(t_1, t_2)|_(t_1 = t_2) = EE[(v(t)-m(t))^2]$ is a particular case of the covariance function, it measures the distance from the mean.
]

#remark[
  Intuitively, when $gamma(tau)$ decays *slowly*, realizations appear *smooth* (slowly varying). When $gamma(tau)$ decays *quickly* (small for $tau >= 1$), realizations appear *nervous* (rapidly fluctuating), approaching white noise behavior.
]

#figure(
  grid(
    columns: 2,
    gutter: 1em,
    cetz.canvas({
      import suiji: gen-rng
      import cetz.draw: *
      import cetz-plot: *

      // Slow decay: AR(1) with a close to 1 → smooth realization
      let rng = gen-rng(42)
      let a = 0.95
      let n = 60
      let points = ()
      let y = 0.0
      for i in range(n) {
        let (r, e) = suiji.uniform(rng, low: -0.3, high: 0.3)
        rng = r
        y = a * y + e
        points.push((i, y))
      }

      plot.plot(
        size: (6, 2.5),
        axis-style: "school-book",
        x-tick-step: none,
        x-label: [$t$],
        y-tick-step: none,
        y-label: [$y(t)$],
        {
          plot.add(points, line: "linear", style: (stroke: blue))
        },
      )
      content((3, -2), [Slow decay ($a = 0.95$): smooth])
    }),
    cetz.canvas({
      import suiji: gen-rng
      import cetz.draw: *
      import cetz-plot: *

      // Fast decay: AR(1) with a close to 0 → nervous realization
      let rng = gen-rng(42)
      let a = 0.1
      let n = 60
      let points = ()
      let y = 0.0
      for i in range(n) {
        let (r, e) = suiji.uniform(rng, low: -1.0, high: 1.0)
        rng = r
        y = a * y + e
        points.push((i, y))
      }

      plot.plot(
        size: (6, 2.5),
        axis-style: "school-book",
        x-tick-step: none,
        x-label: [$t$],
        y-tick-step: none,
        y-label: [$y(t)$],
        {
          plot.add(points, line: "linear", style: (stroke: red))
        },
      )
      content((3, -2), [Fast decay ($a = 0.1$): nervous])
    }),
  ),
  caption: [Effect of covariance decay rate on realization appearance. Left: $gamma(tau)$ decays slowly (high autocorrelation) producing smooth trajectories. Right: $gamma(tau)$ decays quickly producing nervous, rapidly fluctuating trajectories.],
)

#definition(title: "Correlation function")[
  $tilde(gamma)(t_1, t_2) = EE[v(t_1) v(t_2)]$
]

#caution-box[
  The correlation function is similar but *NOT THE SAME* as the covariance function.

  They are related by:
  $ gamma(t_1, t_2) = tilde(gamma)(t_1, t_2) - mu(t_1) mu(t_2) $

  They coincide iff the mean value of the SP is 0
  $ gamma_v(t_1, t_2) eq.triple tilde(gamma)_v(t_1, t_2) <==> m_v = 0 $

  For stationary processes: $tilde(gamma)(tau) = gamma(tau) + mu^2$
]

#definition(title: "Cross-covariance function")[
  The cross-covariance relates two different processes:
  $ gamma_(12)(t_1, t_2) = EE[(v_1(t_1) - mu_1(t_1))(v_2(t_2) - mu_2(t_2))] $
]

#definition(title: "Normalized covariance (Pearson correlation coefficient)")[
  $ rho(t_1, t_2) = gamma(t_1, t_2) / sqrt(gamma(t_1, t_1) dot gamma(t_2, t_2)) $

  Properties:
  - $rho(t, t) = 1$
  - $|rho(t_1, t_2)| <= 1$
  - $|rho| = 1$: maximal correlation
  - $rho = 0$: incorrelation
]

#theorem(title: "Cauchy-Schwarz inequality for covariance")[
  $ |"Cov"[v(t_1), v(t_2)]| <= sqrt("Var"[v(t_1)] dot "Var"[v(t_2)]) $
]

#proof[
  Consider the random vector $bold(v) = vec(v(t_1), v(t_2))$.
  Its variance matrix is:
  $ "Var"[bold(v)] = mat("Var"[v(t_1)], "Cov"[v(t_1), v(t_2)]; "Cov"[v(t_1), v(t_2)], "Var"[v(t_2)]) gt.eq 0 $

  Since this matrix is positive semidefinite, its determinant is non-negative:
  $ "Var"[v(t_1)] dot "Var"[v(t_2)] - "Cov"[v(t_1), v(t_2)]^2 >= 0 $
]

== Stationary Stochastic Processes

#definition(title: "Strongly stationary process")[
  A process is strongly stationary if, for any $n$ instants $t_1, dots, t_n$ and any time shift $T$:
  $ F_(t_1, t_2, dots, t_n)(q_1, dots, q_n) = F_(t_1+T, t_2+T, dots, t_n+T)(q_1, dots, q_n) quad forall T $

  This is a complete but extremely complex characterization.
]

#definition(title: "Weakly stationary stochastic process (SSP)")[
  A SP with
  - $m(t) = m quad forall t$
  - $gamma(t_1, t_2) = gamma(t_1, t_1 + tau) = gamma(tau)$
  Both the expected value and the covariance function do not depend on the specific instant, the covariance function depends on the distance of two instants we call *$tau$* or *lag*.
]

#remark[
  If a process is strongly stationary and its expected value and covariance function exist, then it is also weakly stationary.
]

#example(title: "Non stationary stochastic processes")[
  An example of non stationary stochastic process is a Random Walk:
  - $EE[v(t)]$ is not constant $forall t$
  - $gamma(t)$ is an increasing function
]

#note-box[
  $gamma (t_1, t_2) = gamma(t_3, t_4) "if" t_1 -t_2 = t_3 - t_4 "even if" t_1 != t_2 != t_3 != t_4$
]

#properties(title: "Covariance of SSPs")[
  / positive: $gamma(0) = E[(v(t,s) - m)^2] >= 0$
  / non-increasing: $|gamma(tau)| <= gamma(0) quad forall tau$
  / even function: $gamma(tau) = gamma(-tau) quad forall tau$
]

#proof(title: "Proof that γ(τ) = γ(−τ)")[
  $gamma(tau) = EE[(v(t)-mu)(v(t+tau)-mu)] = EE[(v(t+tau)-mu)(v(t)-mu)] = gamma(-tau)$
]

#properties(title: "Toeplitz matrix property")[
  The covariance matrix built from $gamma(tau)$ is positive semidefinite:
  $
    mat(gamma(0), gamma(1), dots, gamma(N-1); gamma(1), gamma(0), dots, gamma(N-2); dots.v, , dots.down, dots.v; gamma(N-1), dots, , gamma(0)) gt.eq 0
  $

  This can be verified by observing that it equals $EE[bold(V) bold(V)^TT]$ where $bold(V) = vec(v(0)-mu, v(1)-mu, dots.v, v(N-1)-mu)$.

  An even function $gamma(tau)$ is a valid covariance function iff the associated Toeplitz matrix is positive semidefinite for all $N$.
]

=== Gaussian processes

#definition(title: "Gaussian process")[
  A process is Gaussian if, for any $n$ instants $t_1, dots, t_n$, the random variables $v(t_1), dots, v(t_n)$ are jointly Gaussian.
]

#theorem[
  Weakly stationary + Gaussian $arrow.r.double$ Strongly stationary.
]

=== Ergodic processes

#definition(title: "Ergodic process")[
  A stationary process is *ergodic* if statistical properties can be derived (with probability 1) from a single realization.
  $ lim_(N -> infinity) 1 / N sum_(i=1)^N (dot) --> EE[dot] $

  For ergodic processes, time averages converge to ensemble averages.
]



=== White Noise

#definition(title: "White Noise (WN)")[\
  A stationary stochastic process $e(t)$ is called *white noise* with variance $lambda^2$ if:
  - $EE[e(t)] = 0 quad forall t$
  - $gamma_e (tau) = cases(lambda^2 &"if" tau = 0, 0 &"if" tau != 0)$

  We denote this as $e(t) tilde WN(0, lambda^2)$.
]

#remark[
  White noise is *completely unpredictable*: knowing the past of $e(t)$ provides no information about future values. This makes it the fundamental building block for stochastic process models.
]

=== Wold's Decomposition

#theorem(title: "Wold's decomposition theorem")[
  Any stationary process $v(t)$ can be uniquely decomposed as the sum of two uncorrelated components:
  $ v(t) = w(t) + d(t) $
  where:
  - $w(t)$ is a *purely non-deterministic* (PND) component, representable as:
    $ w(t) = sum_(k=0)^(infinity) c_k e(t-k) quad e(t) tilde WN(0, lambda^2) $ with $c_0 = 1$ and $sum_(k=0)^infinity c_k^2 < infinity$
  - $d(t)$ is a *purely deterministic* component that can be perfectly predicted from its past
  - $gamma_(w d)(tau) = 0 quad forall tau$ (the two components are uncorrelated)
]

#remark[
  - The PND component has a representation as an $MA(infinity)$ driven by white noise
  - For most practical applications, the deterministic component is absent ($d(t) = 0$), and $v(t)$ is purely non-deterministic
  - A sinusoidal process $v(t) = A sin(omega_0 t + phi)$ with $phi$ random is an example of a purely deterministic stationary process
]

=== Operatorial representation

#note-box[
  This stuff is explained later in the lectures.
  I moved it here to keep things thematically grouped.
]

We introduce the $z$ operator which will aid us in writing the transfer function of a SP.
This is needed in order to use @thm:stationarity

#definition(title: "Shift operator")[
  The $z^k$ with $k in CC$ operator allows moves expressions from the time-domain to the z-domain and vice versa thereby "shifting" time.

  It is defined as:
  / Backward-shift operator: $ z^(-1) dot x(t) = x(t-1) $
  / Forward-shift operator: $ z^1 dot x(t) = x(t+1) $
]

#note-box[
  In exercises do not mix positive with negative operators.
]

#properties(title: "Shift operator properties")[
  / Linearity: $ z^(-1)(a x(t) + b y(t)) = a x(t-1) + b y(t-1) quad a,b in RR $
  / Recursive application: $ z^(-1)(z^(-1) (z^(-1) x(t))) = z^(-3) x(t) $
  / Linear composition: $ (a z^(-1) + b z + c z^(-3) + d z^(+2)) = a x(t-1) + b x(t+1) + c x(t-3) + d x(t+2) $
]

We can use the $z$ operator to calculate the *transfer function* of various models, we will see it with $AR(m)$, $MA(n)$, $ARMA(m, n)$ and $ARMAX$ models.

Most notably:
- all $AR(m)$, $MA(n)$ and $ARMA(m, n)$ models can be written as
  $ W(z) e(t) = C(z) / A(z) e(t) quad z in CC $
- all $ARMAX$ models can be written as
  $ W(z) e(t) = B(z) / A(z) z^(-k) u(t) + C(z) / A(z) e(t) quad z in CC $

#definition(title: "Discreet time Transfer function of a SP")[
  The transfer function of a SP is the function $W(z)$ such that $y(t) = W(z) e(t)$ where $e(t)$ is a white noise process.
]



#note-box[
  All the properties of generic transfer functions apply:
  - *Series* of filters: $y(t) = M(z)W(z)u(t)$
    #{
      import fletcher: diagram, edge, node
      figure(
        diagram(
          node-shape: "rect",
          node-stroke: 1pt,
          edge((0, 0), "r", "-|>")[$u(t)$],
          node((1, 0))[W(z)],
          edge("r", "-|>")[$s(t)$],
          node((2, 0))[M(z)],
          edge("r", "-|>")[$y(t)$],
        ),
      )
    }

  - *Parallel* of filters: $y(t) = M(z)u(t) plus.minus W(z)u(t)$
  #{
    import fletcher: diagram, edge, node
    figure(
      diagram(
        node-shape: "rect",
        node-stroke: 1pt,
        spacing: (3em, 1em),
        edge((0, 0), "r,u,r", "-|>", label-pos: 15%)[$u(t)$],
        node((2, -1))[W(z)],
        edge("r,d", "-|>", label-pos: 25%)[$w(t)$],

        edge((0, 0), "r,d,r", "-|>"),
        node((2, 1))[M(z)],
        edge("r,u", "-|>", label-pos: 25%)[$m(t)$],

        node((3, 0), shape: "circle", inset: 2pt, sym.plus.minus),
        edge("r", "-|>")[$y(t)$],
      ),
    )
  }

  - Filters in *Feedback configurations*:
    $y(t) = M(z) / (1 minus.plus W(z)M(z)) u(t)$
    #{
      import fletcher: diagram, edge, node
      figure(
        diagram(
          node-shape: "rect",
          node-stroke: 1pt,
          edge((0, 0), "r", "-|>")[$u(t)$],
          node((1, 0), shape: "circle", inset: 2pt, sym.plus.minus),
          edge("r", "-|>")[$epsilon(t)$],
          node((2, 0))[M(z)],
          edge("rr", "-|>", label-pos: 75%)[$y(t)$],
          edge("r,u,l", "-|>"),
          node((2, -1))[W(z)],
          edge("l,d", "-|>"),
        ),
      )
    }
]


=== Stability and stationarity of digital filters

#note-box[
  This stuff is explained later in the lectures.
  I moved it here to keep things thematically grouped.
]

#definition(title: "Singularities of a transfer function")[
  Let $W(z) = C(Z) / A(z)$ be a transfer function.

  We define:
  / zero: $z in CC$ such that $C(z) = 0$
  / pole: $p in CC$ such that $A(p) = 0$
  / singularity: a pole or a zero
]

#theorem(title: "Asymptotic stability of digital filters")[
  A linear digital filter with transfer function $W(z)$ is asymptotically stable iff
  all of its poles are strictly inside the unit circle (in the complex plane).

  #figure(
    cetz.canvas({
      import suiji: gen-rng
      import cetz.draw: *
      import cetz-plot: *
      import "../util.typ": random-pts-in-circle

      let rng = gen-rng(100000)
      let (rng, poles) = random-pts-in-circle(rng, 4, radius: 1.0)


      plot.plot(
        size: (4, 4),
        axis-style: "school-book",
        x-tick-step: none,
        x-label: $Re$,
        x-ticks: (1,),
        x-equal: "y",
        y-min: -1.25,
        y-max: +1.25,
        y-tick-step: none,
        y-label: $Im$,
        {
          // unit circle
          plot.add(domain: (0, 2 * calc.pi), t => (calc.cos(t), calc.sin(t)))


          plot.add(label: "poles", poles, mark: "x", style: (stroke: none))

          plot.add-contour(
            x-domain: (-1, 1),
            y-domain: (-1, 1),
            style: (fill: blue.transparentize(95%), stroke: none),
            fill: true,
            op: "<", // Find contours where data < z
            z: (1, 1, 1), // Z values to find contours for
            (x, y) => calc.sqrt(x * x + y * y),
          )
        },
      )
    }),
  )
]

#theorem(title: "Stationarity of steady-state output")[
  #{
    import fletcher: diagram, edge, node
    figure(
      diagram(
        node-shape: "rect",
        node-stroke: 1pt,
        edge((-1, 0), "r", "-|>")[$e(t)$],
        node((0, 0))[$W(z)$],
        edge((0, 0), "r", "-|>")[$y(t)$],
      ),
    )
  }

  The *steady-state* output (i.e. the output after all initial-condition transients have decayed) $y(t)$ is stationary iff
  - $e(t)$ is SSP
  - $W(z)$ is asymptotically stable
] <thm:stationarity>


#theorem(title: "Stationarity of ARMA(n,m) models")[
  The ARMA(n,m) model represents a SSP if and only if W(z) is asymptotically stable.
]
