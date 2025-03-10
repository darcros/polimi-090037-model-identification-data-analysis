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
  $hat(y) (t+k|t)$ is the is the predictor of SSP $y(t)$.

  // TODO: explain
]

#definition(title: "Linear predictor")[
  $
    hat(y)(t+k|t) &= alpha_0 y(t) + alpha_1 y(t - 1) + dots \
    &= sum_(i=0)^(+infinity) (alpha_i) y(t-i)
  $

  with $alpha_i in RR$ such that $sum(i=0)^(+infinity) alpha_i < +infinity$

  // TODO: note that alpha_i is info from the model and y(t-i) is info from the data
]

#definition(title: "Prediction error")[
  The prediction error of $hat(y)$ is
  $ epsilon(t+k|t) = t(t+k) - hat(y) (t+k|t) $
]

#definition(title: "Mean squared prediction error")[
  The MSP of $hat(y)$ is
  $
    cal(J) (hat(y)) &= EE[epsilon(t+k|t)^2] \
    &= EE[(t(t+k) - hat(y) (t+k|t))^2]
  $
]

#problem(title: "Optimal predictor problem")[
  We would like to find the predictor $hat(y) (t+k|t)$ that minimizes $cal(J)(hat(y))$.

  Since we work with linear (ARMA/ARMAX) models, we will restrict ourselves to linear predictors.
  So we want to find $hat(alpha) = argmin_alpha cal(J) (alpha)$
]

#solution[
  // TODO: add skipped steps

  We reformulate the optimization problem wrt $alpha$ in a new optimization problem wrt $beta$.

  The optimal predictor is given by the minimizer of
  $ EE[(sum_(i=0)^(+infinity) (w_(k+i) - beta_i) e(t-i))^2] >= 0 $
  This means that $hat(beta)_i = w_(k+i), i = 0, 1, dots$

  So our optimal predictor is
  $ hat(y)(t+k|t) = sum_(i=0)^(+infinity) w_(k+i) e(t-i) $
]

#remark[
  This solution is not completely satisfactory as it
  - requires samples of $e(t)$ (we can only sample $y(t)$)
  - requires the use of an infinite number of coefficients
]

// TODO: the part about the optimal predictor analysis
