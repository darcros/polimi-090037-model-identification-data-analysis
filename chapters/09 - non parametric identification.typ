#import "/prelude.typ": *

= Non parametric identification

== Parametric vs non parametric identification
So far, we have explored model based identification. From the true measurement we have created a dataset $D_n = {(u(1), y(1)), (u(2), y(2)), \ldots, (u(N), y(N))}$ and we have used it to identify our Model $W(z)$.

Let's suppose we only need a set of statistics (e.g. $m_y, gamma_y(tau), Gamma_y(omega)$, ... ) regarding the system. Computing all the parameter $theta$ of the model is not necessary. We can directly use the dataset to compute a set of statistics that represent the realization. And hopefully for a large enough dataset, the statistics will be enough to represent the whole system.

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

#proof(title: "Unbiasedness of sample mean")[
  $ EE[hat(mu)_N] = EE[1/N sum_(t=1)^N y(t)] = 1/N sum_(t=1)^N EE[y(t)] = 1/N dot N mu = mu $
]

#proof(title: "Consistency of sample mean")[
  $ "Var"[hat(mu)_N] = "Var"[1/N sum_(t=1)^N y(t)] = 1/N^2 sum_(t_1=1)^N sum_(t_2=1)^N gamma(t_1 - t_2) $
  Since $gamma(tau) -> 0$ as $tau -> infinity$ (for ergodic processes), this variance $-> 0$ as $N -> infinity$.
]

#properties-box(title: "Sample covariance function")[
  $ hat(gamma)_N (tau) = 1 / (N - |tau|) sum_(i=1)^(N - |tau|) y(i) y(i + tau) $
  - $hat(gamma)_N (tau)$ is a *correct estimator* of $gamma^* (tau)$.
  - $hat(gamma)_N (tau)$ is a *consistent estimator* of $gamma^* (tau)$.
  We're going to use $N- |tau|$ samples to compute the covariance function. This is a consequence of the fact that we need to have $y(t)$ and $y(t + tau)$ in the same dataset. The number of samples is going to decrease as $|tau|$increases.
]

#remark(title: "Key insights on sample covariance")[

  *Practical limitation:* Variable sample sizes across lags

  Notice that the sample covariance at different lags is computed over different numbers of samples:
  $
      hat(gamma)_u(0) & = frac(1, N) sum_(t=1)^N y(t)^2 \
      hat(gamma)_u(1) & = frac(1, N-1) sum_(t=1)^(N-1) y(t) y(t+1) \
                      & dots.v \
    hat(gamma)_u(N-1) & = frac(1, 1) y(1) dot y(N)
  $

  We can only estimate $hat(gamma)_u(tau)$ reliably for $tau << N-1$. For large lags, the estimate becomes increasingly noisy due to fewer samples.

  *Symmetry preservation with absolute values*

  By definition: $ gamma(tau) = gamma(-tau) $. This is preserved using $|tau|$ in the denominator:
  $
    hat(gamma)_u(tau) = frac(1, N - |tau|) sum_(t=1)^(N-|tau|) y(t) y(t+|tau|)
  $

  #proof(title: "Unbiasedness of sample covariance")[
    $
      EE[hat(gamma)_u(tau)] &= EE[frac(1, N-tau) sum_(t=1)^(N-tau) y(t) y(t+tau)] \
      &= frac(1, N-tau) sum_(t=1)^(N-tau) EE[y(t) y(t+tau)] = frac(1, N-tau) dot (N-tau) gamma(tau) = gamma(tau)
    $
  ]

  #proof(title: "Consistency of sample covariance")[
    For ergodic processes where $gamma(tau) -->_(tau -> infinity) 0$:
    $
      "Var"[hat(gamma)_u(tau)] = frac(1, (N-tau)^2) sum_(t_1, t_2) text("cov")[y(t_1) y(t_1+tau), y(t_2) y(t_2+tau)] -->_(N -> infinity) 0
    $
  ]
]
#properties-box(title: "Sample spectral density")[
  $ hat(Gamma)_N (omega) = 1 / (2 pi) sum_(tau = -(N - 1))^(N - 1) hat(R)_N (tau) e^(-j omega tau) $
  Here we have two sources of approximation:
  + the approximation of the estimatator covariance function $hat(gamma)_N (tau)$. This makes the sample spectral density an *indirect estimator* of the spectrum $Gamma^* (omega)$.
  + the approximation of the Fourier transform considering only a finite number of samples.
  Some properties still hold:
  - $hat(Gamma)_N (omega)$ is an *asymptotically correct* estimator of $Gamma^* (omega)$.
  $ EE(hat(Gamma)_N (omega)) -->_(N -> infinity) Gamma^* (omega) $
]

#remark(title: "Key insights on sample spectral density")[

  #proof(title: "Asymptotic correctness of sample spectral density")[
    $
      EE[hat(Gamma)_N(omega)] & = frac(1, 2pi) sum_(tau=-(N-1))^(N-1) EE[hat(gamma)_u(tau)] dot e^(-j omega tau) \
                              & = frac(1, 2pi) sum_(tau=-(N-1))^(N-1) gamma(tau) dot e^(-j omega tau) \neq Gamma(omega)
    $
    The true spectrum uses an infinite sum, but we only compute a finite sum with available data. As $N -> infinity$, the truncation error vanishes: $EE[hat(Gamma)_N(omega)] -->_(N -> infinity) Gamma(omega)$.
  ]

  *Consistency and decorrelation across frequencies:*

  As $N -> infinity$, variance decreases and fluctuations at different frequencies become independent:
  $
    EE[(hat(Gamma)_N(omega_1) - Gamma(omega_1))(hat(Gamma)_N(omega_2) - Gamma(omega_2))] -->_(N -> infinity) 0
  $
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

=== Alternative covariance estimator

#definition(title: "Alternative (non-normalized) covariance estimator")[
  $ hat(gamma)'_N (tau) = 1/N sum_(t=1)^(N-|tau|) (y(t) - hat(mu)) (y(t+|tau|) - hat(mu)) $

  This estimator divides by $N$ instead of $N - |tau|$.
  - It is *biased* for finite $N$ (but asymptotically unbiased)
  - However, the resulting Toeplitz matrix is guaranteed to be *positive semidefinite*, which is required for valid spectral estimation
]

=== Periodogram

#definition(title: "Periodogram")[
  $ hat(Gamma)_N^("per")(omega) = 1/N |sum_(t=0)^(N-1) y(t) e^(-j omega t)|^2 $

  The periodogram is the squared magnitude of the DFT of the data, normalized by $N$.
]

#caution-box(title: "Inconsistency of the periodogram")[
  The periodogram is *not a consistent* estimator of the spectrum:
  $ "Var"[hat(Gamma)_N^("per")(omega)] arrow.r.not 0 "as" N -> infinity $

  The variance does not decrease with $N$ — it remains proportional to $Gamma(omega)^2$.
]

=== Improved spectral estimators

#definition(title: "Bartlett's method")[
  Split the $N$ data points into $K$ non-overlapping segments of length $M = N/K$.
  Compute the periodogram for each segment and average:
  $ hat(Gamma)_N^("Bart")(omega) = 1/K sum_(k=1)^K hat(Gamma)_M^("per",k)(omega) $

  This reduces variance by a factor of $K$, but also reduces frequency resolution.
]

#definition(title: "Windowing method")[
  Apply a *window function* $w(tau)$ to the sample covariance before Fourier transforming:
  $ hat(Gamma)_N^("win")(omega) = sum_(tau=-(M-1))^(M-1) w(tau) hat(gamma)_N(tau) e^(-j omega tau) $

  Common windows: Bartlett (triangular), Hann, Hamming, Blackman.

  The window should satisfy:
  - $w(0) = 1$
  - $w(tau) = w(-tau)$ (symmetric)
  - $w(tau) = 0$ for $|tau| >= M$ (finite support)
  - The resulting spectrum estimate should be non-negative
]
