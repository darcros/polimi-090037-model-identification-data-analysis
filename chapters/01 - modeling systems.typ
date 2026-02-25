#import "/prelude.typ": *

= Modeling systems

#{
  import fletcher: diagram, edge, node
  figure(
    diagram(
      node-shape: "rect",
      node-stroke: 1pt,
      edge((0, 0), "rr", "-|>")[$u(t)$],
      node((2, 0))[system],
      edge("r", "-|>")[$y(t)$],

      edge((0, 0), "r,d,r", "-|>"),
      node((2, 1))[model],
      edge("r", "-|>")[$hat(y)(t)$],
    ),
  )
}

We create models to predict the output $y(t)$ of a system.
The model takes the same input $u(t)$ as the system and produces an approximation of the output $hat(y)(t)$.

Models can be:
/ white box: derived from first principles and known physical laws
/ black box: derived from historical data + statistical analysis
/ gray box: derived from a mix of historical data (+ analysis) and simplified physical laws

Some challenges in model identification are:
- *Uncertainty*: #underline[structural] if the physical laws that govern the system are unknown or unavailable and #underline[parametric] if we are oblivious to its parameters.
- *Complexity and model purpose*: the model must be complex enough to describe the model for the intended use but simple enough that computations and predictions are feasible

== Modeling error

#{
  import fletcher: diagram, edge, node
  figure(
    diagram(
      node-shape: "rect",
      node-stroke: 1pt,
      edge((0, 0), "rr", "-|>")[$u(t)$],
      node((2, 0))[system],
      edge("r,d", "-|>", label-pos: 75%)[$+$],

      edge((0, 0), "r,d,r", "-|>"),
      node((2, 1))[model],
      edge("r", "-|>")[$-$],
      node((3, 1), shape: circle, radius: 4pt, stroke: 1pt),
      edge("r", "-|>")[$epsilon(t)$],
    ),
  )
}

== Static and dynamic systems

#definition(title: "Static system")[
  The output $y(t)$ of the system is determined _only_ by its input $u(t)$.

  The system does not change over time.
]

#definition(title: "Dynamic system")[
  The output $y(t)$ of the system is determined by
  - its input $u(t)$
  - the previous state of the system.

  The system changes over time. The change can be
  - *Discrete*: described by difference equations: $y(t) = a dot y(t-1) + e(t)$
  - *Continuous*: described by differential equations $dot(y)(t) = a dot y(t) + e(t)$
]

== Dealing with uncertainty

The natural assumption when building models is that the output can be determined exactly from the input.
This is unrealistic because there are always signals beyond our control (measurement noise, uncontrollable inputs, etc.).
The simplest way to deal with this is to introduce an *additive noise term* into the model:
$ y(t) = f(u(t)) + v(t) $
where $v(t)$ is a disturbance, whose value is _not known a priori_, but past values may allow reasonable estimates of future values using *probability theory*.

== The estimation problem

#definition(title: "Estimation problem")[
  An estimation problem arises whenever there is an unknown parameter $theta$ whose value must be obtained from experimental observations.

  - *Unknown quantity*: $theta$ (scalar or vector, constant or time-varying)
  - *Observations*: $d = {d(t), t in T}$, with $T = {1, 2, dots, N}$
]

#definition(title: "Estimator and estimate")[
  An *estimator* is a function $hat(theta) = f(d)$ that associates to the data a value of the parameter.

  An *estimate* is the numerical value taken by the estimator for specific observed data.
]

#remark[
  Two types of estimation:
  - *Constant parameter* ($theta = "const"$) $arrow.r$ parametric estimation / identification
  - *Time-varying parameter* ($theta(t)$) $arrow.r$
    - $t > t_N$: *prediction*
    - $t = t_N$: *filtering*
    - $t < t_N$: *smoothing* (regularization / interpolation)
]

== The prediction problem

Given a time series $v(1), v(2), dots, v(t-1)$, we want to predict $v(t)$.

#definition(title: "Predictor")[
  $ hat(v)(t|t-1) = f(v(t-1), v(t-2), dots, v(1)) $
  Typically a *linear*, *finite-memory*, *time-invariant* function:
  $ hat(v)(t|t-1) = a_1 v(t-1) + a_2 v(t-2) + dots + a_n v(t-n) $
  with parameter vector $theta = mat(a_1, a_2, dots, a_n)^T$.
]

#definition(title: "Prediction error")[
  $ epsilon(t) = v(t) - hat(v)(t|t-1) $
  The prediction problem becomes an *identification* problem: finding $theta$ that minimizes the cost function
  $ cal(J)(theta) = sum_(i=n+1)^(t-1) epsilon(i)^2 $
]

#remark[
  Minimizing $cal(J)(theta)$ is *necessary but not sufficient*. We must also analyze the structure of the residual $epsilon(t)$:
  - If $epsilon(t)$ has *non-zero mean* $arrow.r$ predictor systematically under/over-estimates
  - If $epsilon(t)$ has *regular patterns* $arrow.r$ further information can be extracted to improve the predictor

  If no regularity remains in $epsilon(t)$, the residual is a *white noise* and the predictor is optimal.
]

A predictor whose residual is white noise leads to the decomposition:
$ v(t) = underbrace(hat(v)(t|t-1), "predictable part") + underbrace(epsilon(t), "white noise") $
which, after substitution, yields the stochastic dynamical system:
$ v(t) = a_1 v(t-1) + a_2 v(t-2) + dots + a_n v(t-n) + epsilon(t) $

== Fundamental elements of identification

The four fundamental elements of a system identification problem are:
+ *System ($cal(S)$)*: the data-generating mechanism (generally a stochastic dynamical system)
+ *Model ($cal(M)$)*: a parametric model with parameter vector $theta$; the choice of model family is a crucial decision
+ *Identification method ($cal(I)$)*: the tool for selecting the best model (typically an optimization procedure)
+ *Identification experiment ($cal(E)$)*: generates data for identification; the experimental conditions should maximize information content in the data
