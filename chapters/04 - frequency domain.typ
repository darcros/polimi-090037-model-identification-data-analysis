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
  $ gamma_y (tau) = cal(F)^(-1){Gamma(omega)} = 1 / (2 pi) integral_(-pi)^(+pi) Gamma(omega) e^(+j omega tau) d w $

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

  Then we can reconduce the spectral density of the output to the spectral density of the input.
  $ Gamma_y(omega) = |W(e^(j omega))|^2 Gamma_v(omega) $
]
#remark(title: "Remark: Gamma function of a white noise")[
  From the definition of *inverse DFT*, if $u(t)$ is a *white noise* process,  the only non-zero term for $gamma_y (tau)$ is at $tau = 0$.
  
   $ Gamma_u (omega) = gamma_u (0) dot e^(- j omega 0) = gamma_u (0)= lambda^2 forall omega $
]

#remark[
  $ W(z)|_(z = e^(j omega)) = W(e^(j omega)) $ is called "Frequency response" of $W(z)$
]

#example(
  title: "Spectral density of a SSP",
)[

$y(t) = (z+c)/(z+a) e(t), |a| < 1$

$e(t) ~ WN(0,lambda^2)$

$ Gamma_e (omega) = gamma_e (0) e^(-j omega 0) = lambda^2$

$Gamma_y (omega) &= |W(e^(j omega))|^2 Gamma_e(omega) \ 
 &= |(e^(j omega) + c) / (e^(j omega) + a)|^2 lambda^2 "we solve the comlex conjugate of the denominator and numerator" \
 &= |((e^(j omega) + c) (e^(-j omega)+c))/ ((e^(j omega) + a)((e^(-j omega) + a)))| lambda^2 " and knowing that" cos(omega) = (e^(j omega) + e^(-j omega))/2, sin(omega) = (e^(j omega) - e^(-j omega))/2j \
  &=( 1 + 2c cos(omega) + c^2) / (1 + 2a cos(omega) + a^2) lambda^2  $
]


#remark(title: "Remark: poles and zeroes on the unit disk and the spectral density")[
  If $e^(j omega)$ is a pole of $W(z)$ then $Gamma_y(omega) = +infinity$

  If $e^(j omega)$ is a zero of $W(z)$ then $Gamma_y(omega) = 0$.

  If $e^(j omega)$ is nearing a *pole* or a *zero* of $W(z)$ then $Gamma_y(omega) -> infinity "or" 0 $ 
]


#note-box( title: "Note: computing the spectrum of a SSP"
)[
  A fast way to compute an approximate value of $Gamma_y(omega)$ is to compute the poles and zeroes of $W(z)$ and use the distances from the poles and zeroes to $e^(j omega)$ to compute the values of $Gamma_y(omega)$ in the most common $Omega$ values. (e.g. $0, pi/2, pi, pi$)

  Assuming our transfer function has one zero and one pole $W(z) = (z-z_("zero")) / (z-z_("pole")) $

  Let
  $ overline(a) = e^(j omega)-z_("zero") $ and $ overline(b) = e^(j omega)-z_("pole") $

  Then
  $ Gamma_y(omega) = |W(e^(j omega))|^2 Gamma_e(omega) = |(e^(j omega)-z_("zero")) / (e^(j omega)-z_("pole"))|^2 lambda^2 = |overline(a) / overline(b)|^2 lambda^2 = |a / b|^2 lambda^2 $
  
  Finally we compute the value of $a$ and $b$ for each common value of $omega$.
  //TODO add figure
]

