#import "/prelude.typ": *

= Model validation
So far we have introduced an algorithm to find the best parameters given a model. We're minimizing the empirical variance of the prediction error.
We need a quality assessment of the model. We need to check if the model is good enough to be used for prediction

$D_N = { (u(1), y(1)), (u(2), y(2)), \ldots, (u(N), y(N)) }$ our dataset of dimension $N$

$M_theta = {M(Theta), theta in Theta subset RR^(n_theta)}$ the model we want to validate.

$hat(Theta)_N = argmin_theta J_N (theta)$ with $J_N (theta) = 1/N sum_{i=1}^N (y(i) - M_theta(u(i)))^2$.

#remark()[
  $hat(theta)$ is going to change with the dataset since it depends on the realization of our loss function $J_N$ and therefore our process $y(t)$ and noise $e(t)$.
]

#figure(
  cetz.canvas({
    import cetz.draw: *

    // Show multiple J_N curves (different realizations) each with its own minimum
    let w = 8
    let h = 4

    // Axes
    line((-0.5, 0), (w + 0.5, 0), stroke: gray + 0.5pt, mark: (end: ">"))
    content((w + 0.8, -0.15), $theta$)
    line((0, -0.3), (0, h + 0.3), stroke: gray + 0.5pt, mark: (end: ">"))
    content((0.5, h + 0.4), $cal(J)_N$)

    // Draw several parabolas with different minima (different realizations)
    let curves = (
      (center: 3.2, scale: 0.15, color: blue, label: $cal(J)_N^((1))$),
      (center: 4.0, scale: 0.12, color: red, label: $cal(J)_N^((2))$),
      (center: 3.6, scale: 0.18, color: green.darken(30%), label: $cal(J)_N^((3))$),
    )

    for curve in curves {
      let pts = ()
      for i in range(0, 41) {
        let x = i * w / 40
        let y = curve.scale * calc.pow(x - curve.center, 2) + 0.3
        pts.push((x, y))
      }
      // Draw curve
      for i in range(1, pts.len()) {
        line(pts.at(i - 1), pts.at(i), stroke: curve.color + 0.8pt)
      }
      // Mark minimum
      circle((curve.center, 0.3), radius: 0.06, fill: curve.color)
      line(
        (curve.center, 0),
        (curve.center, 0.25),
        stroke: (paint: curve.color, dash: "dotted", thickness: 0.6pt),
      )
    }

    // Labels for minima
    content((3.2, -0.3), text(fill: blue, size: 7pt, $hat(theta)_N^((1))$))
    content((4.0, -0.3), text(fill: red, size: 7pt, $hat(theta)_N^((2))$))
    content((3.6, -0.55), text(fill: green.darken(30%), size: 7pt, $hat(theta)_N^((3))$))

    // Theoretical J
    line((0.3, 0.6), (w - 0.3, 0.6), stroke: (paint: black, dash: "dashed", thickness: 0.5pt))
    content((w + 0.1, 0.6), text(size: 7pt, $overline(cal(J))$))
  }),
  caption: [Different realizations of the data produce different empirical cost functions $cal(J)_N^((i))$, each with a different minimizer $hat(theta)_N^((i))$. As $N -> infinity$, they all converge to the theoretical $overline(cal(J))(theta)$.],
)

#theorem()[
  Under current assumptions, as $N -> infinity$:
  $
    J_N (hat(theta), #text(fill: red, "s")) -->_(N -> infinity) dash(J) (theta) = EE_(#text(fill: red, "s"))[epsilon (t|t-1, theta, #text(fill: red, "s"))^2]
  $
  Moreover, by letting
  $ Delta = {theta^*, J(theta^*) <= dash(J)(theta^*) forall theta} $
  be the set of global minima points of $dash(J)(theta^*)$ we have
  $hat(theta_N) (#text(fill: red, "s")) -->_(N -> infinity) Delta$ with $PP(dot) = 1$

]

#corollary()[
  if $Delta= {theta^*}$ is a singleton we have that $ hat(theta_N) (s) -->_(N -> infinity) theta^* "with "PP(dot) = 1 $
]

#remark[
  If $Delta$ is a singleton (the true system belongs to the model class with a unique parameterization), convergence is to a single point. If $cal(S) in.not cal(M)_theta$, the estimate converges to the *best approximation* within the class.
]

#remark(title: "What we really care about")[
  The convergence result above is general but not fully satisfying. We ask ourselves: are we happy with $hat(theta)_N -->_(N -> infinity) Delta$ (the set of minima)?

  What we are *really interested in* is guaranteeing that:
  $
    text("If") quad cal(S) in cal(M)_theta quad (exists theta^0: cal(S) equiv cal(M)(theta^0)), quad text("then") quad hat(theta)_N -->_(N -> infinity) theta^0
  $

  In other words, under the assumption that the true system belongs to the model class (exact correct specification), PEM will recover the true parameters. This is a much stronger and more practical guarantee.
]


#remark(title: "Convergence of the minimum")[
  Furthermore:
  $ min{cal(J)_N (theta)} -->_(N -> infinity) min{overline(cal(J))(theta)} $
  The minimum of the empirical cost converges to the minimum of the theoretical cost.
]

=== PEM converges to the true system

#theorem(title: "PEM convergence when S ∈ M(θ)")[
  If the true system $cal(S) in cal(M)(theta)$ (i.e., there exists $theta^0$ such that $cal(M)(theta^0) = cal(S)$), then PEM guarantees:
  $ hat(theta)_N -->_(N -> infinity) theta^0 $
]

#proof[
  Decompose the prediction error at a generic $theta$:
  $
    epsilon(t|t-1, theta) & = y(t) - hat(y)(t|t-1, theta) \
                          & = y(t) - underbrace(hat(y)(t|t-1, theta^0) + hat(y)(t|t-1, theta^0), 0) - hat(y)(t|t-1, theta) \
                          & = e(t) + [hat(y)(t|t-1, theta^0) - hat(y)(t|t-1, theta)]
  $

  Computing the expected squared error:
  $
    EE[epsilon(t|t-1, theta)^2] &= EE[(e(t) + hat(y)(t|t-1, theta^0) - hat(y)(t|t-1, theta))^2] \
    &= underbrace(EE[e(t)^2], lambda^2) + underbrace(EE[(hat(y)(t|t-1, theta^0) - hat(y)(t|t-1, theta))^2], text("function of θ ≥ 0, min at θ = θ⁰")) + underbrace(2 EE[e(t)(hat(y)(t|t-1, theta^0) - hat(y)(t|t-1, theta))], #text(fill: red, "0: independent of θ")))
  $

  The cross-term vanishes because $e(t)$ is uncorrelated with quantities depending on past data. Thus:
  $ EE[epsilon(t|t-1, theta)^2] = lambda^2 + EE[(hat(y)(t|t-1, theta^0) - hat(y)(t|t-1, theta))^2] >= lambda^2 $

  Equality holds iff $theta = theta^0$. Since PEM converges to the global minimum of $overline(cal(J))(theta) = EE[epsilon^2]$, we get $hat(theta)_N -> theta^0$.
]

#remark[
  This theorem is very significant: if the true system lies within the model set, PEM will converge asymptotically (with enough data) to exactly that system.
]

=== Four cases of PEM convergence

#figure(
  table(
    columns: 3,
    align: (center, left, left),
    table.header([], [$Delta$ singleton], [$Delta$ not singleton]),
    [$cal(S) in cal(M)(theta)$],
    [
      $hat(theta)_N -> theta^0$\
      Unique correct solution.
    ],
    [
      $hat(theta)_N$ converges to one value in $Delta$.\
      Multiple models equally represent $cal(S)$ (over-parameterization).
    ],

    [$cal(S) in.not cal(M)(theta)$],
    [
      $hat(theta)_N -> theta^*$ (best proxy).\
      Model class insufficient but unique best approximation.
    ],
    [
      No guarantees. $hat(theta)_N$ converges to one of the best proxies.\
      Model class too limited and/or data not representative enough.
    ],
  ),
  caption: [Summary of PEM convergence cases depending on whether the true system belongs to the model class and whether the set of global minima $Delta$ is a singleton.],
)

== Model order selection
Let's find the best dimension of $cal(M)_theta$ such that our system $cal(S) in cal(M)_theta$.

The procedure is to sweep over increasing model orders $n = 1, 2, dots, n_max$:
+ Define the model class $cal(M)_theta^((n))$ of order $n$
+ Estimate $hat(theta)_N^((n)) = argmin_theta cal(J)_N^((n))(theta)$
+ Evaluate $cal(J)_N^((n))(hat(theta)_N^((n)))$
+ Select the order that best balances fit quality and model complexity


#figure()[


  #cetz.canvas({
    import suiji: gen-rng
    import cetz.draw: *
    import cetz-plot: *
    import "../util.typ": random-series

    let n_max = 20
    let actual_model_oder = 10

    let lambda2 = 1.5 // asymptotic value of J_N
    let points = range(0, n_max).map(x => (x + 3, 20 / (1 + x) + lambda2))

    plot.plot(
      size: (5, 5),
      axis-style: "school-book",
      x-label: [$n$ model order],
      x-max: n_max,
      y-max: n_max,
      x-tick-step: 50,
      y-tick-step: 50,
      y-label: [$overline(cal(J))_N^((n)) (hat(theta))$],
      x-ticks: ((actual_model_oder, [$n_theta$]), 0),
      y-ticks: ((lambda2, $lambda^2$),),
      {
        plot.add(
          points,
          line: "spline",
        )

        // λ² horizontal asymptote
        plot.add-hline(lambda2, style: (stroke: (dash: "dashed", paint: gray)))

        plot.add-vline(actual_model_oder, style: (stroke: (dash: "dashed")))
        plot.add-fill-between(
          range(0, actual_model_oder + 1).map(x => (x, 0)),
          range(0, actual_model_oder + 1).map(x => (x, 20)),
          style: (fill: rgb("#ff00006c"), stroke: none),
          label: ["underfitting"], // Highlight underfitting in red
        )
        plot.add-fill-between(
          range(actual_model_oder, actual_model_oder + 11).map(x => (x, 0)),
          range(actual_model_oder, actual_model_oder + 11).map(x => (x, 20)),
          style: (fill: rgb("#0000ff46"), stroke: none),
          label: ["overfitting"],
        )
      },
    )
  })
]

#note-box[We see that our loss function $J_N (theta)$ is inversely proportional to the number of parameters in the model. The more parameters we have, the better we can fit the data. But this is not always a good thing since we may be fitting the noise of the data too. To avoid *overfitting*, we need to find a balance between the number of parameters and the goodness of fit.]


Three criteria for model order selections are
+ Whiteness test on the residuals (Anderson's test)
+ Cross validation
+ Identification of the model order penalties.


==== Whiteness test on the residuals (Anderson's test)

#definition(title: "Anderson's test")[
  Given prediction errors $epsilon(t) = y(t) - hat(y)(t|t-1, hat(theta)_N)$, compute its sample auto-covariance:
  $ hat(gamma)_epsilon (tau) = 1/N sum_(t=1)^(N-tau) epsilon(t) epsilon(t+tau) $

  Under the null hypothesis (what we want to reject) $H_0$: $epsilon(t)$ is white noise, for large $N$:
  $
    hat(rho)_epsilon (tau) = hat(gamma)_epsilon (tau) / hat(gamma)_epsilon (0) tilde.dot cal(N)(0, 1/N) quad "for" tau != 0
  $

  The test checks whether $hat(rho)_epsilon (tau)$ falls within the *confidence band*:
  $ |hat(rho)_epsilon (tau)| <= z_(alpha/2) / sqrt(N) quad forall tau = 1, 2, dots, T_max $

  If any value exceeds the band, the model is inadequate at significance level $alpha$.
]

#remark(title: "Practical implementation")[
  Typically $T_max approx N/4$ and $alpha = 0.05$ (95% confidence), giving bands $plus.minus 1.96/sqrt(N)$.

  The test is applied to the *normalized* auto-covariance $hat(rho)(tau)$, not $hat(gamma)(tau)$.
]

#figure(
  grid(
    columns: 2,
    gutter: 2em,
    cetz.canvas({
      import cetz.draw: *

      // Ideal situation
      let w = 5
      let h = 3
      line((-0.3, 0), (w + 0.3, 0), stroke: gray + 0.5pt, mark: (end: ">"))
      content((w + 0.6, -0.15), text(size: 7pt, [$n$ (order)]))
      line((0, -0.3), (0, h + 0.3), stroke: gray + 0.5pt, mark: (end: ">"))

      // Model order labels and test results
      for (i, result) in ((0.8, "KO"), (1.6, "KO"), (2.4, "KO"), (3.2, "OK"), (4.0, "OK")) {
        let color = if result == "KO" { red } else { green.darken(20%) }
        content((i, 0.5), text(fill: color, size: 8pt, weight: "bold", result))
      }

      // Arrow pointing to the best candidate
      line((3.2, 1.0), (3.2, 0.75), stroke: blue + 1pt, mark: (end: ">"))
      content((3.2, 1.3), text(fill: blue, size: 7pt)[best candidate])

      content((w / 2, -0.7), text(size: 8pt, weight: "bold")[Ideal situation])
    }),
    cetz.canvas({
      import cetz.draw: *

      // Real situation
      let w = 5
      let h = 3
      line((-0.3, 0), (w + 0.3, 0), stroke: gray + 0.5pt, mark: (end: ">"))
      content((w + 0.6, -0.15), text(size: 7pt, [$n$ (order)]))
      line((0, -0.3), (0, h + 0.3), stroke: gray + 0.5pt, mark: (end: ">"))

      // All KO
      for i in (0.8, 1.6, 2.4, 3.2, 4.0) {
        content((i, 0.5), text(fill: red, size: 8pt, weight: "bold", "KO"))
      }

      // Question mark
      content((w / 2, 1.3), text(fill: orange, size: 16pt, weight: "bold", "?"))

      content((w / 2, -0.7), text(size: 8pt, weight: "bold")[Real situation])
    }),
  ),
  caption: [Ideal: the whiteness test fails (KO) for under-specified orders and passes (OK) starting at the true order. In practice, with finite data the test may keep failing for all tested orders, making the decision less clear-cut.],
)

===== Why the whiteness test fails in practice

Even when the model is correctly specified, the Anderson test often fails to pass in practice due to several fundamental issues:

#definition(title: "The paradox")[
  According to theory, if $cal(S) in cal(M)(theta)$ and $hat(theta)_N -->_(N -> infinity) theta^0$, then:
  $ epsilon(t|t-1, hat(theta)_N) = y(t) - hat(y)(t|t-1, hat(theta)_N) -->_(N -> infinity) e(t) $
  which should be *asymptotically* white noise. So why does the test fail?
]

#caution-box[The Anderson test is checking *asymptotic* whiteness, but we work with *finite* data. Several practical issues arise:

  1. *Finite-sample correlation*: Even with a correctly specified model, residuals from finite samples contain spurious correlations that would vanish only in the limit $N -> infinity$. The sample autocorrelation $hat(rho)_epsilon(tau)$ has nonzero variance, and randomly some lags will exceed the confidence band just by chance (multiple comparisons problem).

  2. *Estimation error in parameters*: Since $hat(theta)_N != theta^0$ exactly, the residuals $epsilon(t|t-1, hat(theta)_N)$ retain some structure from the true system that cannot be fully removed. This is particularly severe for short time series or high-dimensional models.

  3. *Underspecification*: In many real problems, $cal(S) in.not cal(M)(theta)$ — the true system is not exactly representable by the chosen model class. The residuals will then contain the information lost by the model class, leading to systematic correlations. This is the *common* case in practice.

  4. *Non-stationarity and structural breaks*: Real data often exhibits time-varying properties or abrupt changes that violate the stationarity assumption underlying white noise. The test assumes constant statistics over the entire data window.

  5. *Multiple comparisons*: We test $T_max approx N/4$ different lags against the same confidence band. With $alpha = 0.05$, we expect $approx 0.05 T_max$ false positives just by chance, making the test overly conservative.
]

#remark(
  title: "Practical consequence",
)[The whiteness test is often viewed as a *necessary but not sufficient* condition. A "failure" of the test may indicate missing dynamics, but a "pass" does not guarantee an adequate model. Many practitioners use it as a qualitative guidance tool alongside visual inspection of residuals and other tests, rather than as a definitive decision rule.
]

===== Visualizing whiteness: when does the test pass?

For a correct model (when $cal(S) in cal(M)(theta)$), the residuals should behave as white noise. Let's visualize what this means.

#figure(
  cetz.canvas({
    import cetz.draw: *
    import cetz-plot: *

    let time = range(0, 150)
    // Simulated white noise: truly random without structure
    let wn = time.map(t => (t, 0.8 * calc.cos(t * 0.1 + t * 0.001 * t)))

    plot.plot(
      size: (7, 3),
      axis-style: "school-book",
      x-label: [$t$ (time)],
      y-label: [$epsilon(t) approx WN(0, lambda^2)$],
      x-max: 150,
      y-min: -1.2,
      y-max: 1.2,
      {
        plot.add(wn, line: "linear", fill: false)
      },
    )
  }),
  caption: [If $epsilon(t|t-1, hat(theta)_N) arrow epsilon(t) tilde.op WN(0, lambda^2)$: residuals appear as random fluctuations around zero, indicating an adequate model with no systematic structure left unexplained.],
)

===== Why the whiteness test fails: theory vs practice

#figure(
  grid(
    columns: 2,
    rows: 2,
    gutter: 1.5em,
    // Top-left: Theory - Time domain
    rect(
      fill: green.lighten(95%),
      stroke: green.lighten(50%) + 1.5pt,
      inset: 8pt,
      cetz.canvas({
        import cetz.draw: *

        let w = 4
        let h = 2

        content((w / 2, h + 0.3), text(size: 8pt, weight: "bold")[Theory: Time domain])

        // Time domain section
        line((-0.3, 1.0), (w + 0.3, 1.0), stroke: gray + 0.5pt, mark: (end: ">"))
        content((w + 0.4, 0.95), text(size: 6pt, [$t$]))
        line((0, 0.4), (0, 1.6), stroke: gray + 0.5pt, mark: (end: ">"))
        content((-0.4, 1.6), text(size: 6pt, [$epsilon(t)$]))

        // Random scattered points (white noise in time domain)
        for i in (10, 15, 20, 25, 30, 35, 40, 45, 50, 55) {
          let y = 1.0 + 0.25 * calc.cos(i * 0.3) * calc.exp(-i / 40)
          circle((i / 15, y), radius: 0.04, fill: green.lighten(40%))
        }
        content((w / 2, 0.1), text(size: 7pt, weight: "bold", fill: green.darken(20%), "Random"))
      }),
    ),

    // Top-right: Practice - Time domain
    rect(
      fill: orange.lighten(95%),
      stroke: orange.lighten(50%) + 1.5pt,
      inset: 8pt,
      cetz.canvas({
        import cetz.draw: *

        let w = 4
        let h = 2

        content((w / 2, h + 0.3), text(size: 8pt, weight: "bold")[Practice: Time domain])

        // Time domain section
        line((-0.3, 1.0), (w + 0.3, 1.0), stroke: gray + 0.5pt, mark: (end: ">"))
        content((w + 0.4, 0.95), text(size: 6pt, [$t$]))
        line((0, 0.4), (0, 1.6), stroke: gray + 0.5pt, mark: (end: ">"))
        content((-0.4, 1.6), text(size: 6pt, [$epsilon(t)$]))

        // Points that look scattered but have hidden structure
        for i in (10, 15, 20, 25, 30, 35, 40, 45, 50, 55) {
          let y = 1.0 + 0.2 * calc.sin(i * 0.4) * calc.exp(-i / 35) + 0.1 * calc.cos(i * 0.15)
          circle((i / 15, y), radius: 0.04, fill: orange.lighten(40%))
        }

        content((w / 2, 0.1), text(size: 7pt, weight: "bold", fill: orange.darken(20%), "Appears random"))
      }),
    ),

    // Bottom-left: Theory - Spectrum
    rect(
      fill: green.lighten(95%),
      stroke: green.lighten(50%) + 1.5pt,
      inset: 8pt,
      cetz.canvas({
        import cetz.draw: *

        let w = 4
        let h = 2

        content((w / 2, h + 0.3), text(size: 8pt, weight: "bold")[Theory: Spectrum])

        // Spectrum section
        line((-0.3, 0.5), (w + 0.3, 0.5), stroke: gray + 0.5pt, mark: (end: ">"))
        content((w + 0.4, 0.45), text(size: 6pt, [$omega$]))
        line((0, 0.5), (0, 1.5), stroke: gray + 0.5pt, mark: (end: ">"))
        content((-0.5, 1.5), text(size: 6pt, [$Phi(omega)$]))

        // Spectrum axis marks
        line((-0.06, 0.5), (0.06, 0.5), stroke: gray + 0.3pt)
        content((-0.25, 0.5), text(size: 5pt, $0$))
        line((-0.06, 1.1), (0.06, 1.1), stroke: gray + 0.3pt)
        content((-0.35, 1.1), text(size: 5pt, $lambda^2$))

        // Draw flat spectrum line
        line((0.3, 0.9), (w - 0.3, 0.9), stroke: green.darken(20%) + 1.2pt)

        content((w / 2, 0.1), text(size: 7pt, weight: "bold", fill: green.darken(20%), "Flat"))
      }),
    ),

    // Bottom-right: Practice - Spectrum
    rect(
      fill: orange.lighten(95%),
      stroke: orange.lighten(50%) + 1.5pt,
      inset: 8pt,
      cetz.canvas({
        import cetz.draw: *

        let w = 4
        let h = 2

        content((w / 2, h + 0.3), text(size: 8pt, weight: "bold")[Practice: Spectrum])

        // Spectrum section
        line((-0.3, 0.5), (w + 0.3, 0.5), stroke: gray + 0.5pt, mark: (end: ">"))
        content((w + 0.4, 0.45), text(size: 6pt, [$omega$]))
        line((0, 0.5), (0, 1.5), stroke: gray + 0.5pt, mark: (end: ">"))
        content((-0.5, 1.5), text(size: 6pt, [$Phi(omega)$]))

        // Spectrum axis marks
        line((-0.06, 0.5), (0.06, 0.5), stroke: gray + 0.3pt)
        content((-0.25, 0.5), text(size: 5pt, $0$))
        line((-0.06, 1.3), (0.06, 1.3), stroke: gray + 0.3pt)
        content((-0.35, 1.3), text(size: 5pt, "high"))

        // Draw jagged spectrum with peaks
        let spectrum_pts = (
          (0.3, 0.58),
          (0.5, 0.68),
          (0.7, 0.62),
          (0.9, 0.78),
          (1.1, 0.95),
          (1.3, 0.88),
          (1.5, 0.65),
          (1.7, 0.82),
          (1.9, 0.92),
          (2.1, 0.77),
          (2.3, 0.88),
          (2.5, 1.1),
          (2.7, 0.78),
          (2.9, 0.68),
          (3.1, 0.82),
          (3.3, 0.95),
          (3.5, 0.68),
          (3.7, 0.82),
          (3.9, 0.77),
          (4.1, 0.65),
        )

        for i in range(1, spectrum_pts.len()) {
          line(spectrum_pts.at(i - 1), spectrum_pts.at(i), stroke: orange.darken(20%) + 1pt)
        }

        // Highlight peaks
        circle((1.1, 0.95), radius: 0.06, fill: red.lighten(50%), stroke: red + 0.8pt)
        circle((2.5, 1.1), radius: 0.06, fill: red.lighten(50%), stroke: red + 0.8pt)
        circle((3.3, 0.95), radius: 0.06, fill: red.lighten(50%), stroke: red + 0.8pt)

        content((w / 2, 0.1), text(size: 7pt, weight: "bold", fill: red.darken(20%), "Not flat, not white."))
      }),
    ),
  ),
  caption: [Whiteness test paradox visualized in four subfigures. *Top row*: time-domain residuals appear random in both theory (left, ideal white noise) and practice (right, with hidden structure). *Bottom row*: frequency-domain analysis reveals the difference—theoretical spectrum is flat (left), while practical spectrum shows peaks at specific frequencies (right, circled in red). This demonstrates why spectral analysis is superior to time-domain visual inspection alone for detecting model inadequacy.],
)



==== Cross-correlation test

#definition(title: "Cross-correlation test (for input-output models)")[
  For models with exogenous input, we also test that the residuals are *uncorrelated with the input*:
  $ hat(gamma)_(epsilon u)(tau) = 1/N sum_(t=1)^N epsilon(t) u(t-tau) $

  Under $H_0$ (correct model):
  $ hat(rho)_(epsilon u)(tau) tilde.dot cal(N)(0, 1/N) quad forall tau $

  If the cross-correlation exceeds $plus.minus 1.96/sqrt(N)$, it indicates:
  - Incorrect transfer function $B(z)/A(z)$ (if correlation at specific lags)
  - Wrong delay $d$
  - Missing dynamics
]

#caution-box[
  The Anderson test (auto-correlation of residuals) checks the *noise model* $C(z)/A(z)$.
  The cross-correlation test checks the *input-output model* $B(z)/A(z)$.
  Both tests should be performed for input-output models.
]

==== Simple cross-validation (train/validation split)

The simplest form of cross-validation: partition the dataset into an *identification set* (training) used to build the model and a *validation set* used to assess its performance.

#remark[
  The validation data is "reserved" and not used during model building. This leads to potential *data wastage* --- with less training data, the model may be less accurate, and with less validation data, the performance estimate is noisier.
]

#figure(
  cetz.canvas({
    import cetz.draw: *

    let total_width = 10
    let train_width = 6.5
    let val_width = total_width - train_width

    // Axes
    line((-0.5, 0), (total_width + 0.5, 0), stroke: gray + 0.5pt, mark: (end: ">"))
    content((total_width + 0.8, -0.15), text(size: 8pt, weight: "bold", [$t$ (time)]))

    // Training segment (blue shaded)
    rect((0, -0.3), (train_width, 0.3), fill: blue.lighten(80%), stroke: blue + 1pt)
    content((train_width / 2, 0.8), text(size: 8pt, weight: "bold", fill: blue, "Training set"))
    content((train_width / 2, 0.55), text(size: 7pt, fill: blue, "used to build the model"))

    // Validation segment (red shaded)
    rect((train_width, -0.3), (total_width, 0.3), fill: red.lighten(80%), stroke: red + 1pt)
    content((train_width + val_width / 2, 0.8), text(size: 8pt, weight: "bold", fill: red, "Validation set"))
    content((train_width + val_width / 2, 0.55), text(size: 7pt, fill: red, "performance estimation"))

    // Markers
    line((0, -0.5), (0, -0.35), stroke: black + 0.8pt)
    content((0, -0.75), text(size: 7pt, $1$))

    line((train_width, -0.5), (train_width, -0.35), stroke: black + 0.8pt)
    content((train_width, -0.75), text(size: 7pt, $k$))

    line((total_width, -0.5), (total_width, -0.35), stroke: black + 0.8pt)
    content((total_width, -0.75), text(size: 7pt, $N$))
  }),
  caption: [Simple train/validation split: The dataset is partitioned into a training set (used to estimate parameters) and a held-out validation set (used to estimate generalization performance).],
)


==== Cross validation with K-Fold

We split the dataset into k folds. We train the model on k-1 folds and test it on the last fold. We repeat this process k times, each time using a different fold as the test set. We then compute the average error over all k folds.

#caution-box[

  This is not adequate for time series data since the data is not independent and we would be splicing together different time windows and we would be inserting in the dataset fake temporal correlations.
  We need to use a different approach for time series data.
]

#figure(
  cetz.canvas({
    import cetz.draw: *

    let total = 10
    let fold_width = total / 5
    let row_height = 0.8

    // Title
    content((total / 2, 3.8), text(size: 9pt, weight: "bold", "k-Fold Cross-Validation"))

    // Draw 5 iterations (k=5)
    for iteration in range(0, 5) {
      let y_pos = 3.2 - iteration * row_height

      // Iteration label
      content((-1.0, y_pos), text(size: 7pt, "Iter " + str(iteration + 1)))

      // Draw 5 folds
      for fold in range(0, 5) {
        let x_start = fold * fold_width
        let x_end = x_start + fold_width

        if fold == iteration {
          // Validation fold (red)
          rect((x_start, y_pos - 0.3), (x_end, y_pos + 0.3), fill: red.lighten(70%), stroke: red + 0.8pt)
          content((x_start + fold_width / 2, y_pos), text(size: 6pt, weight: "bold", fill: red, "val"))
        } else {
          // Training fold (blue)
          rect((x_start, y_pos - 0.3), (x_end, y_pos + 0.3), fill: blue.lighten(80%), stroke: blue + 0.5pt)
          content((x_start + fold_width / 2, y_pos), text(size: 6pt, fill: blue, "train"))
        }
      }
    }

    // Legend
    content((total / 2, -0.5), text(size: 7pt, fill: gray.darken(30%), "Average error over all 5 iterations"))
  }),
  caption: [k-Fold cross-validation (k=5): Each iteration uses a different fold as validation (red) while the remaining k-1 folds serve as training data (blue). Performance metrics are averaged over all k iterations.],
)

#caution-box[
  For time series data, k-Fold is problematic because it shuffles the temporal order. Splicing together non-contiguous time windows creates artificial correlations and violates the temporal dependency structure. For dynamical systems, use *time-series-aware* cross-validation methods that preserve temporal continuity (e.g., forward-chaining or rolling windows).
]


==== Cross validation with model order penalties
Instead of minimizing $J_N(theta)$, we can minimize a penalized version of it:

#figure()[


  #cetz.canvas({
    import suiji: gen-rng
    import cetz.draw: *
    import cetz-plot: *
    import "../util.typ": random-series

    let n_max = 20
    let actual_model_oder = 10

    let points = range(0, n_max).map(x => (x + 3, 20 / (1 + x) + 1))

    let AIC_points = range(0, n_max).map(x => (x + 3, calc.ln((50 / (1 + x) + 10) + 2 * x / (n_max))))
    let FPE_points = range(0, n_max).map(x => (x + 3, 20 / (1 + x) + 1 + (n_max + x) / (n_max - x)))
    let MDL_points = range(0, n_max).map(x => (x + 3, calc.ln(100 / (1 + x) + 10) + x * calc.ln(x + 1) / (2 * n_max)))

    plot.plot(
      size: (5, 5),
      axis-style: "school-book",
      x-label: [$n$ model order],
      x-max: n_max,
      y-max: n_max,
      x-tick-step: 50,
      y-tick-step: 50,
      y-label: [$J_N (theta)$],
      x-ticks: ((actual_model_oder, [$n_theta$]), 0),
      {
        plot.add(
          points,
          line: "spline",
          label: ["$J_N(theta)$"],
        )

        plot.add-vline(actual_model_oder, style: (stroke: (dash: "dashed")))
        plot.add(
          AIC_points,
          line: "spline",
          label: ["AIC"],
        )
        plot.add(
          FPE_points,
          line: "spline",
          label: ["FPE"],
        )
        plot.add(
          MDL_points,
          line: "spline",
          label: ["MDL"],
        )
      },
    )
  })
]

+ *FPE*: Final Prediction Error $"FPE"(n) = (N+ n)/(N-n) J_N(hat(theta)^((n)))$
+ *AIC*: Akaike Information Criterion $"AIC"(n) = ln(J_N(hat(theta)^((n)))) + 2n/N$
+ *MDL*: Minimum Description Length $"MDL"(n) = ln(J_N(hat(theta)^((n)))) + ln(N) dot n/(2N)$

#remark(title: "Derivation and interpretation of FPE")[
  The *ideal objective* is to minimize:
  $ min_theta EE[J_N(theta)] = min_theta EE[(1/N sum_(t=1)^N epsilon(t|t-1, theta)^2)] $
  This is not affordable since we don't know the true distribution.

  *FPE* is a numerical approximation:
  $ "FPE"(n) = frac(N + n, N - n) J_N(hat(theta)^((n))) $

  The factor $frac(N + n, N - n)$ corrects for the bias introduced by using the same data to both estimate parameters and evaluate performance. As $N -> infinity$, this factor approaches 1 and FPE approaches $J_N$.
]

#remark(title: "Relationship: FPE vs AIC")[
  The logarithm of FPE is:
  $ ln("FPE") = ln(frac(N + n, N - n)) + ln(J_N(hat(theta)^((n)))) $

  Using the approximation $ln(1 + x) approx x$ for small $x$:
  $ ln(frac(N + n, N - n)) = ln(1 + frac(2n, N - n)) approx frac(2n, N) (quad "for large" N) $

  Therefore:
  $ ln("FPE") approx ln(J_N(hat(theta)^((n)))) + 2n/N = "AIC"(n) $

  This shows that *AIC and FPE return the same optimal model order* asymptotically, making them equivalent choices for large datasets.
]

#remark(title: "Comparison: AIC vs MDL")[
  Both criteria have the form: $[ln(J_N(hat(theta)^((n))))] + [text("penalty term")]$

  - *AIC penalty*: $2n / N$ (constant with $N$)
  - *MDL penalty*: $ln(N) dot n / (2N)$ (grows logarithmically with $N$)

  When $ln(N) > 2$, i.e., $N > e^2 approx 7.39$, the MDL penalty exceeds the AIC penalty, so *MDL penalizes model complexity more strongly*.

  $ "MDL"(n) > "AIC"(n) quad arrow.r.double quad text("MDL tends to select simpler (lower-order) models") $
]

#caution-box[
  *Practical guidance on criteria selection:*
  + If the true system $cal(S) in cal(M)_theta$ and both low-order and higher-order ARX models contain the system, *MDL is theoretically optimal* for identifying the true order.
  + In practice, when $cal(S) in.not cal(M)_theta$ (the true system is not in the model class) or we have multiple competing models, *AIC is often preferred* because it allows for slight over-parameterization, providing additional flexibility and better generalization.
  + For time series with strong dependencies and limited data, *overfitting with AIC is a concern*; MDL is more conservative in these cases.
]

#properties(title: "Theoretical properties")[
  - All three criteria are *consistent* for selecting the correct model order (as $N -> infinity$)
  - *FPE* is derived from minimizing the expected prediction error on future data
  - *AIC* is derived from minimizing the Kullback-Leibler divergence between the true and estimated distributions
  - *MDL* is derived from information theory (minimum description length principle)
  - *FPE* and *AIC* are asymptotically equivalent (return the same optimal order for large $N$)
  - *MDL* penalizes complexity more heavily than *AIC* $arrow.r$ tends to select simpler models
]

#note-box[
  If $cal(S) in cal(M)_theta$, and $cal(M)_theta$ is in the set of *ARX* model, MDL is right.
  If $cal(S) in cal(M)_theta$, and $cal(M)_theta$ is in the set of *ARMAX* model we prefer slightly to overfit so we tend to use AIC.
]
