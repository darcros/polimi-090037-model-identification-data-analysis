#import "/prelude.typ": *

= Model validation
So far we have introduced an algorrithm to find the best parameters given a model. We're minimizing the empirical variance of the prediction error.
We need a quality assesment of the model. We need to check if the model is good enough to be used for prediction

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
  Under current assumptionss, as $N -> infinity$:
  $ J_N (hat(theta), s) -->_(N -> infinity) dash(J) (theta) = EE[epsilon (t|t-1, theta, s)^2] $
  Moreover, by letting
  $ Delta = {theta^*, J(theta^* <= dash(J)(theta^*) forall theta} $
  be the set of global minima points of $dash(J)(theta^*)$ we have
  $hat(theta_N) (s) -->_(N -> infinity) Delta$ with $PP(dot) = 1$
]

#corollary()[
  if $Delta= {theta^*}$ is a singleton we have that $ hat(theta_N) (s) -->_(N -> infinity) theta^* "with "PP(dot) = 1 $
]

#remark[
  If $Delta$ is a singleton (the true system belongs to the model class with a unique parameterization), convergence is to a single point. If $cal(S) in.not cal(M)_theta$, the estimate converges to the *best approximation* within the class.
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
                          & = y(t) - hat(y)(t|t-1, theta^0) + hat(y)(t|t-1, theta^0) - hat(y)(t|t-1, theta) \
                          & = e(t) + [hat(y)(t|t-1, theta^0) - hat(y)(t|t-1, theta)]
  $

  Computing the expected squared error:
  $
    EE[epsilon(t|t-1, theta)^2] &= EE[(e(t) + hat(y)(t|t-1, theta^0) - hat(y)(t|t-1, theta))^2] \
    &= underbrace(EE[e(t)^2], lambda^2) + EE[(hat(y)(t|t-1, theta^0) - hat(y)(t|t-1, theta))^2] + underbrace(2 EE[e(t)(hat(y)(t|t-1, theta^0) - hat(y)(t|t-1, theta))], 0)
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
  Given prediction errors $epsilon(t) = y(t) - hat(y)(t|t-1, hat(theta)_N)$, compute the sample auto-covariance:
  $ hat(gamma)_epsilon (tau) = 1/N sum_(t=1)^(N-tau) epsilon(t) epsilon(t+tau) $

  Under the null hypothesis $H_0$: $epsilon(t)$ is white noise, for large $N$:
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

==== Cross validation with k fold
We split the dataset into k folds. We train the model on k-1 folds and test it on the last fold. We repeat this process k times, each time using a different fold as the test set. We then compute the average error over all k folds.

This is not adequeate for time series data since the data is not independent and we would be splicing together different time windows and we would be inserting in the dataset fake temporal correlations.
We need to use a different approach for time series data.

==== Cross validation with model order penalties
Instead of minimizing J_N(theta), we can minimize a penalized version of it:

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
          label: ["J_N(theta)"],
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

+ *FPE*: Final Prediction Error $"FPE" = (N+ n)/(N-n) J_N(theta)$
+ *AIC*: Akaike Information Criterion $"AIC"(n) = ln(J_N(theta)) + 2n/(N)$
+ *MDL*: Minimum Description Length $"MDL"(n) = ln(J_N(theta)) + n ln(N)/(2N)$

#properties(title: "Theoretical properties")[
  - All three criteria are *consistent* for selecting the correct model order (as $N -> infinity$)
  - *FPE* is derived from minimizing the expected prediction error on future data
  - *AIC* is derived from minimizing the Kullback-Leibler divergence between the true and estimated distributions
  - *MDL* is derived from information theory (minimum description length principle)
  - *MDL* penalizes complexity more heavily than AIC $arrow.r$ tends to select simpler models
  - For $N -> infinity$: MDL penalty $>$ AIC penalty, so MDL selects lower-order models
]

#note-box[
  If $cal(S) in cal(M)_theta$, and $cal(M)_theta$ is in the set of *ARX* model, MDL is right.
  If $cal(S) in cal(M)_theta$, and $cal(M)_theta$ is in the set of *ARMAX* model we prefer slightly to overfit so we tend to use AIC.
]
