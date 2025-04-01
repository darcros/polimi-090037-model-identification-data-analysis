#import "/prelude.typ": *

= Non parametric identification

== Parametric vs non parametric identification
So far, we have explored model based identification. From the true measurement we have created a dataset $D_n = {(u(1), y(1)), (u(2), y(2)), \ldots, (u(N), y(N))}$ and we have used it to identify our Model $W(z)$.

Let's suppose we only need a set of statistics reguarding the system. Computing all the parameter $theta$ of the model is not necessary. We can directly use the dataset to compute a set of statistics that represent the realization. And hopefully for a large enough dataset, the statistics will be enough to represent the whole system.

Let's first define some properties of the estimators we are going to use.
#definition(title: "Correctness")[
  $hat(T)_N$ is a correct estimator of $T^*$ if the expectation of the estimator is equal to the true value of the parameter, i.e. $ EE[hat(T)_N] = T^* $.
]

#definition(title: "Consistency")[
  $hat(T)_N$ is a consistent estimator of $T^*$ if the variance of the estimator goes to zero as the number of samples goes to infinity, i.e. $ "VAR"[hat(T)_N] = EE[(hat(T)_N - T^*)^2] -->_(N -> infinity) 0 $.
]

