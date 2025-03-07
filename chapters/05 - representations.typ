#import "/prelude.typ": *

= Representations
We have seen how we can represent SP with at least 4 representation types:
- *Time domain*
$y(t) = a_1 y(t) + a_2 y(t-1) + ... + a_n y(t-m) + c_0 e(t) + b_1 e(t-1) + ... + c_n u(t-n)$

- *Operatorial* (Transfer function)
$y(t) = C(z) / A(z) e(t)$

- *Probabilistic*
$m_y(t) = E[y(t)]$

$gamma_y(tau) = E[(y(t) - m_y(t)) (y(t-tau) - m_y(t-tau))]$

- *Frequency domain*
$m_y(t) = E[y(t)]$

$Gamma(omega) = sum_(-infinity)^(+infinity) gamma_y(t) e^(-j omega tau)$

#remark[
  Reliable forecasts are based on unique representations.
  It is possible to represent a class of equivalent SP in different representations, but for each of them there are infinitely many representations $W(z)$, $e(t)$ that produce a specific given $y(t)$ as outcome.

  Some cases follow:

  Cases:

  - $y(t) = W(z) e(t)$
  1. $y(t) = W(z) alpha / alpha e(t)
      cases(
    tilde(W)(z) = W(z) / alpha \
    tilde(e)(t) = alpha e(t)
  )$

  2. $y(t) = W(z) z^n / z^n e(t)
      cases(
    tilde(W)(z) = W(z) z^n \
    tilde(e)(t) = z^{-n} e(t) = e(t-n)
  )$

  3. $y(t) = W(z) (z-p) / (z-p) e(t) , p in CC, abs(p)<1
      cases(
    tilde(W)(z) = W(z) (z-p) / (z-p) \
    tilde(e)(t) = e(t)
  )$
  has a different difference equation than the original one but same dynamics.

  4. $y(t) = W(z) 1 / q (z-q) / (z-1 / q) e(t) , p in CC, abs(p)>1$

  Which is dimostrably equivalent $Gamma_y(t) = abs(Gamma(e^{j omega}))^2 = Gamma_e(omega)$

]

For each of these possible cases, we can define properties so to define a SSP in a unique way, selecting a *Canonical representation*.

#theorem(title: "Canonical representation")[
  A SSP is uniquely defined by a *canonical representation*.

  Let $y(t)$ be a SSP with a *rational* $Gamma_y(omega) in QQ$.

  Then there exists a unique representation $W(z)$, $e(t)$ such that $y(t) = W(z) e(t)$ and \ $W(z)= C(z)/A(z), z in CC$ if the following conditions hold:
  + $$C(z), A(z)$, $ are *monic*: the leading coefficient is 1 \ $cases(
    C(z) = 1 + c_1z^{-1} + ... + c_n z^{-n} \
    A(z) = 1 + a_1z^{-1} + ... + a_m z^{-m}
    )$
  + $C(z), A(z)$ have *null relative degree*: $nu = n- m$
  + $C(z), A(z)$ are *coprime*: the greatest common divisor is 1
  + $C(z), A)(z)$ have *roots in the unit circle*: $abs(z_i) = 1, i = 1, ..., n$

]
