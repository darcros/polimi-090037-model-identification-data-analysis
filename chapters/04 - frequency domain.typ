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

#remark(title: "Existence condition")[
  The Fourier series converges (and $Gamma(omega)$ exists) iff:
  $ sum_(tau = -infinity)^(+infinity) |gamma(tau)| < infinity $
  i.e. $gamma(tau)$ is *absolutely summable*.
]

#properties[
  The function $Gamma(omega)$ is

  / real: $Im(Gamma(omega)) = 0, forall omega in RR$

  / positive: $Gamma(omega) >= 0, forall omega in RR$

  / even: $Gamma(omega) = Gamma(-omega), forall omega in RR$

  / $2 pi$ periodic: $Gamma(omega) = Gamma(omega + 2 pi k), forall omega in RR, k in ZZ$
]

#proof(title: "Proof that Γ(ω) ≥ 0")[
  For processes with rational spectrum: $Gamma(omega) = |W(e^(j omega))|^2 lambda^2 >= 0$, since it is a product of non-negative terms.
]

#corollary(title: "Variance as area under the spectrum")[
  Setting $tau = 0$ in the inverse DFT formula:
  $ gamma(0) = "Var"[v(t)] = 1/(2pi) integral_(-pi)^(+pi) Gamma(omega) d omega $
  The variance equals the total area under the spectrum (divided by $2 pi$).
]

#lemma(title: "Shannon sampling bound")[
  The minimum time period that can be represented by $Gamma(omega)$ in is two samples.

  In fact
  $
                 omega & <= omega_max \
    omega = (2 pi) / T & <= pi = omega_max \
                     T & >= (2 pi) / omega_max = (w pi) / pi = 2 = T_min
  $
]

#properties(title: "Inverse DFT")[
  $ gamma_y (tau) = cal(F)^(-1){Gamma(omega)} = 1 / (2 pi) integral_(-pi)^(+pi) Gamma(omega) e^(+j omega tau) d w $

  Note the exponent of $e$, it is *positive*.
]

#theorem(title: "spectral factorization")[
  #{
    import fletcher: diagram, edge, node
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
  From the definition of *inverse DFT*, if $u(t)$ is a *white noise* process, the only non-zero term for $gamma_y (tau)$ is at $tau = 0$. Therefore, the spectral density of a white noise process is constant.

  $ Gamma_u (omega) = gamma_u (0) dot e^(- j omega 0) = gamma_u (0)= lambda^2 forall omega $
]

#remark[
  $ W(z)|_(z = e^(j omega)) = W(e^(j omega)) $ is called "Frequency response" of $W(z)$
]

#example(title: "Spectral density of a SSP")[

  $y(t) = (z+c) / (z+a) e(t), |a| < 1$

  $e(t) ~ WN(0, lambda^2)$

  $Gamma_e (omega) = gamma_e (0) e^(-j omega 0) = lambda^2$

  $Gamma_y (omega) &= |W(e^(j omega))|^2 Gamma_e(omega) \
  &= |(e^(j omega) + c) / (e^(j omega) + a)|^2 lambda^2 "we solve the comlex conjugate of the denominator and numerator" \
  &= |((e^(j omega) + c) (e^(-j omega)+c)) / ((e^(j omega) + a)((e^(-j omega) + a)))| lambda^2 " and knowing that" cos(omega) = (e^(j omega) + e^(-j omega)) / 2, sin(omega) = (e^(j omega) - e^(-j omega)) / 2j \
  &=( 1 + 2c cos(omega) + c^2) / (1 + 2a cos(omega) + a^2) lambda^2$
]


#remark(title: "Remark: poles and zeroes on the unit disk and the spectral density")[
  If $e^(j omega)$ is a pole of $W(z)$ then $Gamma_y(omega) = +infinity$

  If $e^(j omega)$ is a zero of $W(z)$ then $Gamma_y(omega) = 0$.

  If $e^(j omega)$ is nearing a *pole* or a *zero* of $W(z)$ then $Gamma_y(omega) -> infinity "or" 0$
]


#note-box(title: "Note: computing the spectrum of a SSP")[
  A fast way to compute an approximate value of $Gamma_y(omega)$ is to compute the poles and zeroes of $W(z)$ and use the distances from the poles and zeroes to $e^(j omega)$ to compute the values of $Gamma_y(omega)$ in the most common $Omega$ values. (e.g. $0, pi/2, pi, pi$)

  Assuming our transfer function has one zero and one pole $W(z) = (z-z_("zero")) / (z-z_("pole"))$

  Let
  $ overline(a) = e^(j omega)-z_("zero") $ and $ overline(b) = e^(j omega)-z_("pole") $

  Then
  $
    Gamma_y(omega) = |W(e^(j omega))|^2 Gamma_e(omega) = |(e^(j omega)-z_("zero")) / (e^(j omega)-z_("pole"))|^2 lambda^2 = |overline(a) / overline(b)|^2 lambda^2 = |a / b|^2 lambda^2
  $

  Finally we compute the value of $a$ and $b$ for each common value of $omega$.

  #figure(
    cetz.canvas({
      import cetz.draw: *

      // Unit circle
      circle((0, 0), radius: 2.0, stroke: gray + 0.8pt)

      // Axes
      line((-2.5, 0), (2.5, 0), stroke: gray + 0.5pt)
      line((0, -2.5), (0, 2.5), stroke: gray + 0.5pt)
      content((2.7, -0.2), $Re$)
      content((0.3, 2.6), $Im$)

      // Pole at 0.7 on real axis
      let pole = (1.4, 0)
      content(pole, $times$)
      content((1.4, -0.35), text(size: 8pt, $z_"pole"$))

      // Zero at -0.5 on real axis
      let zero = (-1.0, 0)
      circle(zero, radius: 0.1, stroke: 1pt)
      content((-1.0, -0.35), text(size: 8pt, $z_"zero"$))

      // Point on unit circle at omega = pi/3
      let omega = calc.pi / 3
      let pt = (2.0 * calc.cos(omega), 2.0 * calc.sin(omega))
      circle(pt, radius: 0.06, fill: black)
      content((pt.at(0) + 0.4, pt.at(1) + 0.2), $e^(j omega)$)

      // Distance vectors
      line(zero, pt, stroke: (paint: blue, thickness: 1.2pt, dash: "dashed"))
      content((-0.1, 1.2), text(fill: blue, size: 8pt, $overline(a)$))

      line(pole, pt, stroke: (paint: red, thickness: 1.2pt, dash: "dashed"))
      content((1.5, 1.2), text(fill: red, size: 8pt, $overline(b)$))

      // Arc for omega
      arc((0, 0), start: 0deg, stop: 60deg, radius: 0.6, stroke: 0.5pt)
      content((0.55, 0.25), text(size: 7pt, $omega$))
    }),
    caption: [Graphical method: $Gamma_y(omega) prop |overline(a) \/ overline(b)|^2$ where $overline(a)$ and $overline(b)$ are distances from the zero and pole to $e^(j omega)$ on the unit circle.],
  )
]

== Complex spectrum and input-output relations

#definition(title: "Complex spectrum (bilateral Z-transform of covariance)")[
  $ Phi(z) = sum_(tau = -infinity)^(+infinity) gamma(tau) z^(-tau) $

  The power spectral density is obtained by evaluating on the unit circle: $Gamma(omega) = Phi(e^(j omega))$.
]

#definition(title: "Cross-covariance and cross-spectrum")[
  Given two processes $u(t)$ and $y(t)$ with $y(t) = W(z) u(t)$:

  $
    gamma_(u y)(tau) = sum_(i=0)^infinity w(i) gamma_(u u)(tau - i) quad arrow.r.double quad Phi_(u y)(z) = W(z) Phi_(u u)(z)
  $

  $
    gamma_(y y)(tau) = sum_(i=0)^infinity w(i) gamma_(y u)(tau - i) quad arrow.r.double quad Phi_(y y)(z) = W(z) Phi_(y u)(z)
  $

  Since $Phi_(y u)(z) = Phi_(u y)(z^(-1))$, we get:
  $ Phi_(y y)(z) = W(z) W(z^(-1)) Phi_(u u)(z) $
]

#remark[
  The cross-spectrum $Gamma_(u y)(omega)$ is in general *complex* (unlike auto-spectra which are real).
]

#remark(title: "MIMO systems")[
  For multiple-input multiple-output systems:
  $
    Phi_(y y)(z) = W(z) Phi_(u u)(z) W(z^(-1))^TT, quad Gamma_(y y)(omega) = W(e^(j omega)) Gamma_(u u)(omega) W(e^(-j omega))^TT
  $
]

== Spectrum examples

#example(title: "Spectrum of MA(1)")[
  Let $v(t) = e(t) + c e(t-1)$, $e(t) tilde WN(0, lambda^2)$.

  Using the definition directly (only $gamma(0), gamma(plus.minus 1) != 0$):
  $
    Gamma(omega) = gamma(0) + gamma(1) e^(-j omega) + gamma(-1) e^(j omega) = (1 + c^2) lambda^2 + 2c lambda^2 cos(omega)
  $
]

#example(title: "Spectrum of MA(2)")[
  Let $v(t) = e(t) + c_1 e(t-1) + c_2 e(t-2)$.

  $ Gamma(omega) = (1 + c_1^2 + c_2^2) lambda^2 + 2(c_1 + c_1 c_2) lambda^2 cos(omega) + 2 c_2 lambda^2 cos(2 omega) $
]

#example(title: "Spectrum of AR(1)")[
  Let $v(t) = a v(t-1) + e(t)$, $|a| < 1$.

  Using the spectral factorization with $W(z) = 1/(1 - a z^(-1))$:
  $ Gamma(omega) = |W(e^(j omega))|^2 lambda^2 = lambda^2 / (1 + a^2 - 2a cos(omega)) $
]

== Spectrum of deterministic processes

#remark(title: "Wold decomposition in the frequency domain")[
  From Wold's decomposition $v(t) = w(t) + d(t)$:
  - The PND component $w(t)$ has a *continuous* spectrum distributed over all frequencies
  - The purely deterministic component $d(t)$ has a spectrum consisting of *impulses* (Dirac deltas) at specific frequencies
]

#example(title: "Constant process")[
  $v(t) = v(t-1) = v$, where $v tilde cal(N)(0, lambda^2)$.

  $tilde(gamma)(tau) = lambda^2 forall tau$, so $Gamma(omega) = lambda^2 delta(omega)$ (impulse at frequency 0).
]

#example(title: "Sinusoidal process")[
  $v(t) = v cos(omega_0 t + theta)$, where $v tilde cal(N)(mu_v, lambda^2)$, $theta tilde cal(U)(-pi, pi)$ independent.

  This process is stationary with $EE[v(t)] = 0$ and $tilde(gamma)(tau) = lambda^2/2 cos(omega_0 tau)$.

  $ Gamma(omega) = lambda^2/4 [delta(omega - omega_0) + delta(omega + omega_0)] $
]

#example(title: "Mixed process")[
  If $v(t) = tilde(v)(t) + overline(v)(t)$ with $tilde(v)$ deterministic (spectrum $Gamma_(tilde(v))$) and $overline(v)$ PND (spectrum $Gamma_(overline(v))$), and the two components are independent:
  $ Gamma_v (omega) = Gamma_(tilde(v))(omega) + Gamma_(overline(v))(omega) $
]

