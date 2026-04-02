#import "/prelude.typ": *

= Frequency domain

Since there are infinitely many realizations for a SP we cannot simply analyze the signal in the frequency domain.
So we analyze the covariance function $gamma$ instead.

#figure(
  rect(
    stroke: 0.5pt,
    grid(
      columns: 3,
      align: center + horizon,
      gutter: 1.5em,
      [
        *Covariance function*\
        time domain ($tau$)\
        $gamma(tau)$
      ],
      [
        #text(fill: red)[Fourier Transform\
          (Discrete)]\
        $arrow.double.l.r$
      ],
      [
        *Spectrum*\
        frequency domain ($omega$)\
        $Gamma(omega)$
      ],
    ),
  ),
  caption: [The covariance function $gamma(tau)$ and the spectral density $Gamma(omega)$ are a discrete Fourier transform pair.],
)

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

#figure(
  cetz.canvas({
    import cetz.draw: *

    // Axes
    line((-5, 0), (5, 0), stroke: black + 0.8pt)
    line((0, -0.5), (0, 4), stroke: black + 0.8pt)

    // Axis labels
    content((5.3, -0.3), $omega$, anchor: "west")
    content((-0.3, 4.2), $Gamma(omega)$, anchor: "south")

    // Axis ticks and labels
    line((-calc.pi, -0.2), (-calc.pi, 0.2), stroke: black + 0.5pt)
    content((-calc.pi, -0.6), $-pi$, anchor: "north", font-size: 9pt)

    line((calc.pi, -0.2), (calc.pi, 0.2), stroke: black + 0.5pt)
    content((calc.pi, -0.6), $pi$, anchor: "north", font-size: 9pt)

    line((0, -0.2), (0, 0.2), stroke: black + 0.5pt)

    // Draw smooth spectrum curve using many segments
    let spectrum(omega) = {
      let a = 0.6
      let lambda2 = 1.0
      return lambda2 / (1.0 + a * a - 2.0 * a * calc.cos(omega)) * 1.0
    }

    let points = ()
    for i in range(0, 201) {
      let omega = -calc.pi + i * 2.0 * calc.pi / 200.0
      let val = spectrum(omega)
      points.push((omega * 0.8, val + 0.3))
    }

    // Draw hatching under the curve
    for i in range(0, 200) {
      let x = -calc.pi * 0.8 + i * calc.pi * 1.6 / 200.0
      let y = spectrum(x / 0.8) * 1.0 + 0.3
      line((x, 0.3), (x, y), stroke: purple + 0.5pt)
    }

    // Draw curve
    let prev = none
    for pt in points {
      if prev != none {
        line(prev, pt, stroke: purple + 1.5pt)
      }
      prev = pt
    }

    // Annotation for the area
    content(
      (0, 1.8),
      $gamma(0) dot 2pi$,
      anchor: "south",
      fill: white,
      padding: 2pt,
      stroke: purple + 0.5pt,
      font-size: 10pt,
    )
    line((0, 1.7), (0, 1.1), stroke: purple + 1pt, mark: (end: ">"))
  }),
  caption: [The variance $gamma(0) = "Var"[y(t)]$ equals the area under the spectrum curve: the shaded region equals $gamma(0) dot 2pi = integral_(-pi)^(pi) Gamma(omega) dif omega$.],
)


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

#properties(title: "Linearity of the spectrum")[
  / Scalar multiple: If $y(t) = a x(t)$, then $ Gamma_y (omega) = a^2 Gamma_x (omega) $
  / Sum of uncorrelated processes: If $z(t) = x(t) + y(t)$ and $x, y$ are uncorrelated, then $ Gamma_z (omega) = Gamma_x (omega) + Gamma_y (omega) $
]

#example(title: "Linear combination of uncorrelated processes")[
  If $v(t) = a x(t) - b y(t)$ with $x(t)$ and $y(t)$ uncorrelated:
  $ Gamma_v (omega) = a^2 Gamma_x (omega) + b^2 Gamma_y (omega) $
]

#note-box(title: "Euler formula recall")[
  Useful identities for spectrum computation:
  - $e^(j omega) = cos omega + j sin omega$
  - $e^(-j omega) = cos omega - j sin omega$
  - $e^(-j omega) + e^(j omega) = 2 cos(omega)$
  - $e^(j omega) - e^(-j omega) = 2 j sin(omega)$

  On the unit circle: $e^(j omega)$ and $e^(-j omega)$ are complex conjugates. The first rotates counterclockwise by angle $omega$, decomposing into $cos omega$ (real) and $sin omega$ (imaginary). Its conjugate is reflected across the real axis.
]

#figure(
  cetz.canvas(length: 1cm, {
    import cetz.draw: *

    // Axes
    line((-2.3, 0), (2.3, 0), stroke: 0.6pt)
    line((0, -2.3), (0, 2.3), stroke: 0.6pt)

    // Axis labels
    content((2.45, -0.12), [Re], anchor: "west", size: 7pt)
    content((-0.12, 2.45), [Im], anchor: "south", size: 7pt)

    // Unit circle
    let steps = 60
    for i in range(steps) {
      let angle1 = i * 360deg / steps
      let angle2 = (i + 1) * 360deg / steps
      let x1 = 2.0 * calc.cos(angle1)
      let y1 = 2.0 * calc.sin(angle1)
      let x2 = 2.0 * calc.cos(angle2)
      let y2 = 2.0 * calc.sin(angle2)
      line((x1, y1), (x2, y2), stroke: 0.5pt)
    }

    // e^(jω) at 50 degrees
    let omegaa = 50deg
    let x_pos = 2.0 * calc.cos(omegaa)
    let y_pos = 2.0 * calc.sin(omegaa)

    // e^(-jω) reflected
    let x_neg = 2.0 * calc.cos(omegaa)
    let y_neg = -2.0 * calc.sin(omegaa)

    // Lines from origin
    line((0, 0), (x_pos, y_pos), stroke: red + 0.8pt)
    line((0, 0), (x_neg, y_neg), stroke: green + 0.8pt)

    // Projections on axes
    line((x_pos, 0), (x_pos, y_pos), stroke: (paint: gray, thickness: 0.4pt, dash: "dashed"))
    line((0, y_pos), (x_pos, y_pos), stroke: (paint: gray, thickness: 0.4pt, dash: "dashed"))

    // Points
    circle((x_pos, y_pos), radius: 0.06, fill: red, stroke: none)
    circle((x_neg, y_neg), radius: 0.06, fill: green, stroke: none)

    // Labels
    content((x_pos + 0.2, y_pos + 0.15), $e^(j omega)$, anchor: "south-west", size: 6pt)
    content((x_neg + 0.2, y_neg - 0.15), $e^(-j omega)$, anchor: "north-west", size: 6pt)
    content((0.3, 0.1), $omega$, anchor: "south-west", size: 6pt)

    // Axis projections labels
    content((x_pos, -0.2), $cos omega$, anchor: "north", size: 6pt)
    content((-0.25, y_pos), $sin omega$, anchor: "east", size: 6pt)
  }),
  caption: [Unit circle: $e^(j omega)$ and its conjugate $e^(-j omega)$.],
)





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
  $ Gamma_y(omega) = |W(e^(j omega))|^2 Gamma_v(omega) = W(e^(j omega)) W(e^(-j omega)) Gamma_v(omega) $

  Note: $W(e^(j omega))$ is called the frequency response of the filter W(z).
]
#remark(title: "Remark: Gamma function of a white noise")[
  From the definition of *inverse DFT*, if $u(t)$ is a *white noise* process, the only non-zero term for $gamma_y (tau)$ is at $tau = 0$. Therefore, the spectral density of a white noise process is constant.

  $ Gamma_u (omega) = gamma_u (0) dot e^(- j omega 0) = gamma_u (0)= lambda^2 forall omega $
]

#figure(
  cetz.canvas({
    import cetz.draw: *

    // Axes
    line((-5, 0), (5, 0), stroke: black + 0.8pt)
    line((0, -0.3), (0, 2), stroke: black + 0.8pt)

    // Axis labels
    content((5.3, -0.2), $omega$, anchor: "west")
    content((-0.3, 2.1), $Gamma(omega)$, anchor: "south")

    // Axis ticks and labels
    line((-calc.pi, -0.15), (-calc.pi, 0.15), stroke: black + 0.5pt)
    content((-calc.pi, -0.4), $-pi$, anchor: "north", font-size: 9pt)

    line((calc.pi, -0.15), (calc.pi, 0.15), stroke: black + 0.5pt)
    content((calc.pi, -0.4), $pi$, anchor: "north", font-size: 9pt)

    line((0, -0.15), (0, 0.15), stroke: black + 0.5pt)

    // Constant spectrum for white noise
    let lambda2 = 1.0
    let wn-spectrum = 0.3 + lambda2 * 0.5

    // Draw constant horizontal line
    line((-calc.pi * 0.8, wn-spectrum), (calc.pi * 0.8, wn-spectrum), stroke: purple + 1.5pt)

    // Draw hatching under the (constant) line
    for i in range(0, 200) {
      let x = -calc.pi * 0.8 + i * calc.pi * 1.6 / 200.0
      line((x, 0.3), (x, wn-spectrum), stroke: purple + 0.5pt)
    }

    // Annotation for constant spectrum
    content(
      (calc.pi * 0.9, wn-spectrum),
      $lambda^2$,
      anchor: "west",
      fill: white,
      padding: 2pt,
      stroke: purple + 0.5pt,
      font-size: 10pt,
    )
  }),
  caption: [The spectrum of white noise is *constant* across all frequencies: $Gamma(omega) = lambda^2$ for all $omega$. This indicates that white noise contains all frequencies with equal power.],
)




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

== Final value theorem

#theorem(title: "Theorem of the gain (final value theorem)")[
  For a stable transfer function $W(z)$:
  $ EE[y(t)] = W(1) dot EE[u(t)] $

  This holds even if $W(z)$ is not in canonical form.
]

#remark(title: "De-biasing procedure")[
  The final value theorem is useful when dealing with non-zero-mean processes:

  + Find the mean using the gain theorem: $EE[y(t)] = W(1) dot mu = overline(y)$
  + Define de-biased processes:
    $ tilde(y)(t) = y(t) - overline(y), quad tilde(e)(t) = e(t) - mu $
  + Compute covariance etc. on the de-biased process: $gamma_(tilde(y))(tau) = gamma_y (tau) quad forall tau$
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

