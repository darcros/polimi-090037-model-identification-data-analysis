#import "/prelude.typ": *

= Non parametric identification

== Parametric vs non parametric identification
So far, we have explored model based identification. From the true measurement we have created a dataset $D_n = {(u(1), y(1)), (u(2), y(2)), \ldots, (u(N), y(N))}$ and we have used it to identify our Model $W(z)$.

Let's suppose we only need a set of statistics reguarding the system. Computing all the parameter $theta$ of the model is not necessary. We can directly use the dataset to compute a set of statistics that represent the realization. And hopefully for a large enough dataset, the statistics will be enough to represent the whole system.

Let's first define some properties of the estimators we are going to use.
#definition(title: "Correctness")[
  $hat(T)_N$ is a *correct estimator* of $T^*$ if the expectation of the estimator is equal to the true value of the parameter, i.e. $ EE[hat(T)_N] = T^* $
  *Asymptotic correctness* is defined as $ EE[hat(T)_N] -->_(N -> infinity) T^* $
]

#definition(title: "Consistency")[
  $hat(T)_N$ is a *consistent estimator* of $T^*$ if the variance of the estimator goes to zero as the number of samples goes to infinity, i.e. $ "VAR"[hat(T)_N] = EE[(hat(T)_N - T^*)^2] -->_(N -> infinity) 0 $.
]

Some of the most common non parametric estimators are:
- *Sample mean*
- *Sameple covariance function*
- *Sample spectral density*


#properties-box(title: "Sample mean")[
  $ hat(mu)_N = 1 / N sum_(i=1)^N y(i) $
  - $hat(mu)_N$ is a *correct* estimator of $mu^*$.
  - $hat(mu)_N$ is a *consistent* estimator of $mu^*$.
]
#properties-box(title: "Sample covariance function")[
  $ hat(gamma)_N (tau) = 1 / (N - |tau|) sum_(i=1)^(N - |tau|) y(i) y(i + tau) $
  - $hat(gamma)_N (tau)$ is a *correct estimator* of $gamma^* (tau)$.
  - $hat(gamma)_N (tau)$ is a *consistent estimator* of $gamma^* (tau)$.
  We're going to use N- |tau| samples to compute the covariance function. This is a consequence of the fact that we need to have $y(t)$ and $y(t + tau)$ in the same dataset. The number of samples is going to decrease as $|tau| $increases.
]
#properties-box(title: "Sample spectral density")[
  $ hat(Gamma)_N (omega) = 1 / (2 pi) sum_(tau = -(N - 1))^(N - 1) hat(R)_N (tau) e^(-j omega tau) $
  Here we have two sources of approximation:
  + the approximation of the estimatator covariance function $hat(gamma)_N (tau)$. This makes the sample spectral density an *indirect estimator* of the spectrum $Gamma^* (omega)$.
  + the approximation of the Fourier transform considering only a finite number of samples.
  Some properties still hold:
  - $hat(Gamma)_N (omega)$ is an *asymptotically correct* estimator of $Gamma^* (omega)$.
  - $ EE(hat(Gamma)_N (omega)) -->_(N -> infinity) Gamma^* (omega) $.
]

#note-box(title: "Other estimators of the spectral denstilty")[
  The sample spectral density is an indirect estimator of the spectrum and usually is computationally expensive to compute as it is an indirect estimator and uses another estimator (the covariance function) to compute the spectrum.

  We can solve this problem by using another estimator of the spectrum. The *discrete Fourier transform* is a direct estimator of the spectrum. It is defined as:
  $
    hat(Gamma)_N (omega) = 1 / N sum_(t=0)^(N - 1) |y(t) e^(-j omega t)|^2
  $

  - We don't use the covariance function to compute the spectrum. We use directly the dataset.
  - The DFT is a *direct estimator* of the spectrum.
  - The DFT is *no longer a correct estimator* of the spectrum. The expectation of the DFT is not equal to the spectrum. But for large enough datasets, the DFT is a *asymptotically consistent estimator* of the spectrum.
]
