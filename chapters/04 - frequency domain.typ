#import "/prelude.typ": *

= Frequency domain

Since there are infinitely many realizations for a SP we cannot simply analyze the signal in the frequency domain.
So we analyze the covariance function $gamma$ instead.

#definition(title: "Spectral density of a SSP")[
  The spectral density (or spectrum) of a SSP $y$ is

  $ Gamma(omega) = sum_(-infinity)^(+infinity) gamma_y(t) e^(-j omega tau) $

  where
  - $j$ imaginary unit
  - $omega$ frequency
]

#remark[
  $Gamma(omega) eq.triple cal(F){gamma_y(t)}$ is the Discrete Fourier Transform (DFT) of the covariance function $gamma_y(t)$.
]

#properties[
  The function $Gamma(omega)$ is

  / real: $Im(Gamma(omega)) = 0, forall omega in RR$

  / positive: $Gamma(omega) >= 0, forall omega in RR$

  / even: $Gamma(omega) = Gamma(-omega), forall omega in RR$

  / $2 pi$ periodic: $Gamma(omega) = Gamma(omega + 2 pi k), forall omega in RR, k in ZZ$
]

// FIXME: is it a Lemma, Corollary or proposition?
#lemma(title: "Application of the Shannon Theorem")[
  The minimum time period that can be represented by $Gamma(omega)$ in is two samples.

  In fact
  $
    omega &<= omega_max \
    omega = (2 pi) / T &<= pi = omega_max \
    T &>= (2 pi) / omega_max = (w pi) / pi = 2 = T_min
  $
]

// FIXME: is this really a definition? or does this just follow from the definition of Gamma?
#definition(title: "Inverse DFT")[
  $ gamma_y(tau) = cal(F)^(-1){Gamma(omega)} = 1 / (2 pi) integral_(-pi)^(+pi) Gamma(omega) e^(+j omega tau) d w $

  Note the exponent of $e$, it is *positive*.
]

#theorem(title: "spectral factorization")[
  #{
    import fletcher: diagram, node, edge
    figure(
      diagram(
        node-shape: "rect",
        node-stroke: 1pt,
        edge((-1, 0), "r", "-|>")[$v(t)$],
        node((0, 0))[$W(z)$],
        edge((0, 0), "r", "-|>")[$y(t)$],
      ),
    )
  }

  Where
  - $v(t)$ a SSP
  - $W(z)$ an asymptotically stable filter

  Then $ Gamma_y(omega) = |W(e^(j omega))|^2 Gamma_v(omega) $
]

#remark[
  $ W(z)|_(z = e^(j omega)) = W(e^(j omega)) $ is called "Frequency response" of $W(z)$
]

// TODO: examples

