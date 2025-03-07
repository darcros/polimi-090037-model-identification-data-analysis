#import "/prelude.typ": *

= Stochastic processes

Stochastic processes extend the notion of random variable to signals.

#definition(title: "Stochstic process")[
  A stochastic process is an infinite sequence of random variables, all defined on the same probabilistic space

  $ v(1, s), v(2, s), dots, v(t, s) $

  where $t$ is the time index.
  // FIXME: is $s$ the probabilisti space or something else?
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
      size: (10,2),
      axis-style: "school-book",
      x-min: -1,
      x-max: +20,
      x-tick-step: none,
      x-label: [$t$],
      y-tick-step: none,
      y-label: [$v(t, overline(s))$],
      {
        plot.add(label: "realization", realization, line: "spline")
      }
    )
  }),
  caption: [If $s$ is fixed we get a deterministic function]
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
      let a = 1 / (sigma * calc.sqrt(2  * calc.pi))
      
      let sigma2 = calc.pow(sigma, 2)
      let exponent = -calc.pow(x - mu , 2) / (2 * sigma2)
      
      return a * calc.exp(exponent)
    }
    
    plot.plot(
      size: (10,2),
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
      }
    )
  }),
  caption: [If $t$ is fixed we get a random variable]
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
    import "../util.typ": random-series, avg-series
    let rng = gen-rng(1)

    let realizations = ()
    let i = 0;
    while (i < 5) {
      let (r, realization) = random-series(rng, 20, from: -1, to: 20)
      realizations.push(realization)

      rng = r
      i += 1
    }

    let mean = avg-series(realizations)
    
    plot.plot(
      size: (10,2),
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
      }
    )
  }),
  caption: [Mean value and a bunch of realizations]
)

#definition(title: "Covariance function of a SP")[
  $gamma(t_1, t_2) = EE[(v(t_1)-m(t_1))(v(t_2)-m(t_2))] $
  
  The #underline[variance function] $gamma(t) = gamma(t_1, t_2)|_(t_1 = t_2) = EE[(v(t)-m(t))^2]$ is a particular case of the covariance function, it measures the distance from the mean.
]

// TODO: graphical representation of the covariance function and how it relates to the realizations (slow vs. nervous realization)

#definition(title: "Correlation function")[
  // FIXME: I don't know which simbol to use so i called it "corr"
  $"corr"(t_1, t_2) = EE[v(t_1) v(t_2)]$
]

#caution-box[
  The correlation function is similar but *NOT THE SAME* as the covariance function.

  They coincide iff the mean value of the SP is 0
  $ gamma_v(t_1, t_2) eq.triple "corr"_v(t_1, t_2) <==> m_v = 0 $
]

== Stationary Stochastic Processes

#definition(title: "Stationary stochastic process")[
  A SP with
  - $m(t) = m quad forall t$
  - $gamma(t_1, t_2) = gamma(t_1, t_1 + tau) = gamma(tau)$
  Both the expected value and the covariance function do not depend on the specific instant, the covariance function depends on the distance of two instants we call *$tau$* or *lag*.
]

#example(title: "Non stationary stochastic processes")[
  An example of non stationary stochastic process is a Random Walk:
  - $EE[v(t)]$ is not constant $forall t$
  - $gamma(t)$ is an increasing function
]

#note-box[
  $gamma (t_1, t_2) = gamma(t_3, t_4) "if" t_1 -t_2 = t_3 - t_4 "even if" t_1 != t_2  != t_3 != t_4$
]

#properties(title: "Covariance of SSPs")[
  / positive: $gamma(0) = E[(v(t,s) - m)^2] >= 0$
  / non-increasing: $|gamma(tau)| <= gamma(0) quad forall tau$
  // FIXME: lo ho anche io così negli appunti ma non dovrebbe essere per t e t + tau?
  / even function: $gamma(tau) = gamma(-tau) quad forall tau$
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

#properties(title: "Shift operator properties")[
  / Linearity: $ z^(-1)(a x(t) + b y(t)) = a x(t-1) + b y(t-1) quad a,b in RR $
  / Recursive application: $ z^(-1)(z^(-1) (z^(-1) x(t))) = z^(-3) x(t) $
  / Linear composition: $ (a z^(-1) + b z + c z^(-3) + d z^(+2)) = a x(t-1) + b x(t+1) + c x(t-3) + d x(t+2) $
]

We can use the $z$ operator to calculate the *transfer function* of various models, we will see it with $"AR"(m)$, $"MA"(n)$, $"ARMA"(m, n)$ and $"ARMAX"$ models.

Most notably:
- all $"AR"(m)$, $"MA"(n)$ and $"ARMA"(m, n)$ models can be written as
  $ W(z) e(t) = C(z) / A(z) e(t) quad z in CC $
- all $"ARMAX"$ models can be written as
  $ W(z) e(t) = B(z) / A(z) z^(-k) u(t) + C(z) / A(z) e(t) quad z in CC $

#note-box[
  All the properties of generic transfer functions apply:
  - *Series* of filters: $y(t) = M(z)W(z)u(t)$
    #{
      import fletcher: diagram, node, edge
      figure(
        diagram(
          node-shape: "rect",
          node-stroke: 1pt,
          edge((0, 0), "r", "-|>")[$u(t)$],
          node((1,0))[W(z)],
          edge("r", "-|>")[$s(t)$],
          node((2,0))[M(z)],
          edge("r", "-|>")[$y(t)$],
        )
      )
    }

  - *Parallel* of filters: $y(t) = M(z)u(t) plus.minus W(z)u(t) $
  #{
    import fletcher: diagram, node, edge
    figure(
      diagram(
        node-shape: "rect",
        node-stroke: 1pt,
        spacing: (3em, 1em),
        edge((0, 0), "r,u,r", "-|>", label-pos: 15%)[$u(t)$],  
        node((2,-1))[W(z)],
        edge("r,d", "-|>", label-pos: 25%)[$w(t)$],
        
        edge((0, 0), "r,d,r", "-|>"),
        node((2,1))[M(z)],
        edge("r,u", "-|>", label-pos: 25%)[$m(t)$],
        
        node((3,0), shape: "circle", inset: 2pt, sym.plus.minus),
        edge("r", "-|>")[$y(t)$],
      )
    )
  }

- Filters in *Feedback configurations*:
  $y(t) = M(z) / (1 minus.plus W(z)M(z)) u(t) $
  #{
    import fletcher: diagram, node, edge
    figure(
      diagram(
        node-shape: "rect",
        node-stroke: 1pt,
        edge((0, 0), "r", "-|>")[$u(t)$],
        node((1, 0), shape: "circle", inset: 2pt, sym.plus.minus),
        edge("r", "-|>")[$epsilon(t)$],
        node((2,0))[M(z)],
        edge("rr", "-|>", label-pos: 75%)[$y(t)$],
        edge("r,u,l", "-|>"),
        node((2, -1))[W(z)],
        edge("l,d", "-|>"),
      )
    )
  }
]

  
=== Stationarity theorem // FIXME: I don't know what it's actually called

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

// TODO: what was the name of the theorem?
#theorem(title: "Asymptotic stability for linear digital filters")[
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
        size: (4,4),
        axis-style: "school-book",
        x-tick-step: none,
        x-label: $Re$,
        x-ticks: (1, ),
        x-equal: "y",
        y-min: -1.25,
        y-max: +1.25,
        y-tick-step: none,
        y-label: $Im$,
        {
          // unit circle
          plot.add(domain: (0, 2 * calc.pi), (t) => (calc.cos(t),calc.sin(t)))

          
          plot.add(label:"poles", poles, mark: "x", style: (stroke: none))

          plot.add-contour(x-domain: (- 1,  1), y-domain: (- 1, 1 ),
          style: (fill: blue.transparentize(95%), stroke: none),
          fill: true,
          op: "<", // Find contours where data < z
          z: (1,1,1), // Z values to find contours for
          (x, y) => calc.sqrt(x * x + y * y))       
        }
      )
    })
  )
]

// TODO: what was the name of the theorem?
#theorem(title: "Stationarity for steady state outputs of digital filters")[
  #{
    import fletcher: diagram, node, edge
    figure(
      diagram(
        node-shape: "rect",
        node-stroke: 1pt,
        edge((-1, 0), "r", "-|>")[$e(t)$],
        node((0, 0))[$W(z)$],
        edge((0, 0), "r", "-|>")[$y(t)$],
      )
    )
  }

  // TODO: define steady-state
  The steady-state output $y(t)$ is stationary iff
  - $e(t)$ is SSP
  - $W(z)$ is asymptotically stable
] <thm:stationarity>