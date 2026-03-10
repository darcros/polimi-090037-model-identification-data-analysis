#import "/prelude.typ": *

= Model identification

#problem(title: "Identification problem")[
  Given a system $cal(S)$, we collect samples of input $u(t)$ and output $y(t)$ in order to derive a model of the system.

  #{
    import fletcher: diagram, edge, node
    figure(
      diagram(
        spacing: (1em, 1em),
        // debug: true,
        {
          edge((-2, 0), "rr", "-|>")[$u(t)$]
          node((0, 0), stroke: 1pt)[system $cal(S)$]
          edge("rr", "-|>")[$y(t)$]

          node((-1, 1))[$u(1) \ u(2) \ dots.v \ u(N)$]
          edge("d,r", "=>")
          node((+1, 1))[$y(1) \ y(2) \ dots.v \ y(N)$]
          edge("d,l", "=>")

          node((0, 2))[Machine learning]
          edge("d", "=>")
          node((0, 3))[Mathematical \ model of $cal(S)$]
        },
      ),
    )
  }
]

== Parametric model identification

=== Static modeling framework

Before considering dynamic models, it is instructive to consider the static case.

#definition(title: "Static model")[
  $ y(t) = f(x(t), theta) + e(t) $
  where $x(t)$ are known regressors (basis functions), $theta$ are parameters, $e(t)$ is noise.
]

#remark(title: "Basis function expansion")[
  For linear-in-parameters models: $y(t) = phi(t)^TT theta + e(t)$, where $phi(t) = vec(f_1(x(t)), f_2(x(t)), dots.v, f_n(x(t)))$.

  The basis functions $f_i$ can be polynomials, radial basis functions, wavelets, etc. Neural networks provide a form of *universal approximation*.
]

#remark(title: "Occam's razor")[
  Model complexity should be as low as possible while still capturing the essential dynamics of the system. Overly complex models overfit the data and generalize poorly.
]

The parametric approach to model identification consists of the following steps:
+ Experiment design (and data collection)
+ Selection of a parametric model class $cal(M)(theta) = {cal(M)_theta : theta in Theta}$
+ Choice of the identification criterion $cal(J)_N (theta) >= 0$
+ Minimization of $cal(J)_N (theta)$ wrt $theta$ and computation of $hat(theta)_N$
  $ hat(theta)_N = argmin_theta cal(J)_N(theta) $
+ Model validation

=== Experiment design

When designing the experiment we need to:
- Design the input $u(t)$.
  The information content of the output $y(t)$ depends on this.

- Select a suitable length for the experiment (number of collected samples $N$).
  With longer experiments we can be more confident in our data.

=== Choice of the model class

#definition[
  A model class
  $ cal(M)(theta) = {cal(M)_theta : theta in Theta} $
  is a set of models $cal(M)_theta$ which are parametrized over $theta$.
  $Theta$ is the set of admissible values for $theta$.
]

#remark[
  There are many options to choose from:
  - #highlight[discrete time] vs. continuous time
  - #highlight[linear] vs. non-linear
  - #highlight[time-invariant] vs. time-varying
  - static vs. #highlight[dynamic]

  In this course we focus in the #highlight[highlighted] ones.
]

#example[
  Consider the class of $ARMAX(m, n, k, p)$ models
  $
    mat(
      delim: #none, align: #left, row-gap: #1em, column-gap: #0.5em,
      y(t) = #h(1em), , #hide[+] a_1 y(y-1), + a_2 y(t-2), dots, + a_m y(t-m);
      , +b_0 u(t-k), +b_1 u(t-k-1), + b_2 u(t-k-2), dots, + b_m u(t-k-p);
      , +c_0 e(t), +c_1 e(t-1), + c_2 e(t-2), dots, + c_m e(t-n);
    )
  $
  where $e(t) ~ WN(0, lambda^2)$.

  Then
  $ theta = mat(a_1, a_2, dots, b_0, b_1, b_2, dots, c_0, c_1, c_2, dots)^T in Theta subset.eq RR^(m + n + p) $

  #note-box[
    $theta$ is not sufficient to fully specify the model.
    We also need to know
    - $lambda^2$: the variance of the white noise
    - $k$: the pure delay of the system
    - m, n, p: the degrees of the $AR$, $MA$ and $"X"$ parts

    The first two are not critical in the sense that they can be obtained "for free" through the minimaztion of our identification criterion (see @model-ident:step-3).
    While the latter must be specified ahead a priori.
  ]

  #note-box[
    $Theta != RR^(m + n + p)$ because some values in $RR^(m + n + p)$ are not admissible since we admit only canonical representations (?).

    Consider $y(t) = a_1 y(t-1) + e(t) ~ AR(1) = ARMAX(1, 0, 0, 0)$

    In order for $y(t)$ to be static, we need $|a_1| < 1$ (see @thm:stationarity).

    So $Theta = {a : |a| < 1} subset RR^1$ and our model class would be
    $cal(M)_theta = {AR(1), theta = a, |a| < 1}$
  ]
]

=== Choice of the identification criterion <model-ident:step-3>

We use a _predictive approach_: a model is considered good if it can accurately predict future values. The key idea is that if we found the "true" model, the prediction error would be white noise, and no further improvement would be possible.

Our ideal objective would be
$ cal(J)(theta) = EE[(y(t+1) - hat(y)(t+1|t, theta))^2] $
but that is not computable, because we only have *one* realization.

So we use a sample-based version instead
$ cal(J)(theta) = 1 / N sum_(t=1)^N (y(t+1) - hat(y)(t|t-1, theta))^2 $

#remark(title: "Why 1-step ahead prediction")[
  We "look" 1-step ahead (instead of $2, 3, dots, k$) because

  $ epsilon(t+k|t) = y(t+k) - hat(y)(t+k|t) $
  $ "var"[e(t)] <= "var"[epsilon(t+k|t)] <= "var"[y(t)] $

  if the model is exact and we compute the optimal 1-step ahead error
  $ epsilon(t+k|t) = e(t) $
  and
  $ "var"[epsilon(t+1|t)] = "var"[e(t)] = lambda^2 $

  Then
  $ cal(J)(hat(theta)) = EE[e(t)^2] = lambda^2 $
  where $hat(theta)$ is the parameter vector of optimal model.

  So we effectively get $lambda^2$ for free.
] <remark:free-lambda>

=== Minimization of the identification criterion

We will consider two situations:
- $cal(J)(theta)$ is quadratic, this is the case for $AR$ and $"ARX"$ models
- $cal(J)(theta)$ is non-quadratic, this is the case for $ARMA$ and $ARMAX$ models

#figure(
  grid(
    columns: 2,
    gutter: 2em,
    cetz.canvas({
      import cetz.draw: *

      // Quadratic (ARX): paraboloid contours, single minimum
      let cx = 0
      let cy = 0
      line((-3, -1.5), (3, -1.5), stroke: gray + 0.5pt, mark: (end: ">"))
      content((3.3, -1.5), $theta_1$)
      line((0, -2.5), (0, 2.5), stroke: gray + 0.5pt, mark: (end: ">"))
      content((0.4, 2.6), $theta_2$)

      for r in (0.4, 0.9, 1.4, 2.0) {
        circle((cx, cy), radius: (r * 1.3, r * 0.8), stroke: blue.lighten(30%) + 0.6pt)
      }
      circle((cx, cy), radius: 0.08, fill: red)
      content((0.9, 0.3), text(fill: red, size: 8pt)[unique min])
      content((0, -3.0), text(size: 8pt)[*Quadratic* (ARX): unique global minimum])
    }),
    cetz.canvas({
      import cetz.draw: *

      // Non-quadratic (ARMAX): multiple local minima
      line((-3, -1.5), (3, -1.5), stroke: gray + 0.5pt, mark: (end: ">"))
      content((3.3, -1.5), $theta_1$)
      line((0, -2.5), (0, 2.5), stroke: gray + 0.5pt, mark: (end: ">"))
      content((0.4, 2.6), $theta_2$)

      // Multiple "basins"
      for c in ((-1.5, -0.5), (1.2, 0.8), (-0.3, 1.3)) {
        for r in (0.3, 0.6, 0.9) {
          circle(c, radius: (r * 0.8, r * 0.5), stroke: orange.lighten(30%) + 0.5pt)
        }
        circle(c, radius: 0.06, fill: red)
      }
      // Global minimum label
      content((-1.5, -1.0), text(fill: red, size: 7pt)[global])
      content((1.2, 0.3), text(fill: red, size: 7pt)[local])
      content((-0.3, 0.8), text(fill: red, size: 7pt)[local])
      content((0, -3.0), text(size: 8pt)[*Non-quadratic* (ARMAX): multiple local minima])
    }),
  ),
  caption: [Quadratic cost (ARX) has a unique global minimum found by OLS. Non-quadratic cost (ARMAX) may have multiple local minima, requiring iterative methods.],
)

=== Model validation

During the process we made assumptions that must be verified:
- The system $cal(S)$ lies within the model class $cal(M)_theta$ (correct structure and order). If this fails, the identified model is the *best approximation* within the class, but may exhibit systematic bias.
- The model orders $m$, $n$, $p$ have been "correctly" selected. If orders are too low, the model cannot capture the true dynamics; if too high, overfitting occurs (good fit on training data, poor generalization).
- The input is persistently exciting of sufficient order. If this fails, the parameter estimate may be non-unique.

== Model families and prediction forms

#definition(title: "Model families")[
  Common model families, classified by structure:
  / Time series models: No exogenous input. AR, ARMA.
  / Input-output models: With exogenous input. ARX, ARMAX, Output Error (OE), FIR, Box-Jenkins.

  The key distinction is *how noise enters the model*:
  - *Equation error* (ARX): noise enters additively after $A(z)$
  - *Output error* (OE): noise is on the output only
  - *ARMAX*: noise is filtered by $C(z)/A(z)$
]

#definition(title: "FIR model")[
  The *Finite Impulse Response* model is a special case of ARX with $A(z) = 1$:
  $ y(t) = B(z) u(t-d) + e(t) $
  It has no feedback and is always stable, but may require many parameters.
]

#definition(title: "Prediction form")[
  Every model can be written in *prediction form* for 1-step ahead prediction:
  $ hat(y)(t|t-1, theta) = phi(t, theta)^TT theta $

  - For ARX: $phi(t)$ depends only on past data $arrow.r$ *linear in data and parameters*
  - For ARMAX: $phi(t, theta)$ depends on $theta$ through $C(z)$ $arrow.r$ *nonlinear in parameters*
]

== Prediction Error Methods (PEM)

#definition(title: "PEM framework")[
  The PEM approach minimizes:
  $ cal(J)_N(theta) = 1/N sum_(t=1)^N epsilon(t, theta)^2, quad epsilon(t, theta) = y(t) - hat(y)(t|t-1, theta) $

  The prediction error $epsilon(t, theta)$ depends on the model family:
  - *ARX*: $epsilon(t, theta) = A(z) y(t) - B(z) u(t-d)$ $arrow.r$ linear in $theta$
  - *ARMAX*: $epsilon(t, theta) = A(z)/C(z) y(t) - B(z)/C(z) u(t-d)$ $arrow.r$ nonlinear in $theta$

  #{
    import fletcher: diagram, edge, node
    figure(
      diagram(
        spacing: (2em, 2em),
        node-stroke: 1pt,
        node-shape: "rect",
        {
          // System S
          node((0, 0))[System $cal(S)$]
          edge((-1.5, 0), (0, 0), "-|>", label: $u(t)$, label-side: center, label-sep: 0.5em)
          edge((0, 0), (2, 0), "-|>", label: $y(t)$, label-side: center, label-sep: 0.5em)

          // Model
          node((0, 1))[Model $hat(cal(M))(theta)$]
          edge((-1.5, 0), (-1.5, 1), (0, 1), "-|>")
          edge((0, 1), (1, 1), "-|>", label: $hat(y)(t|t-1, theta)$, label-side: center, label-sep: 0.5em)

          // Difference node
          node((2, 0.5), stroke: 0.8pt, shape: "circle", width: 1.5em, height: 1.5em)[$-$]
          edge((2, 0), (2, 0.5), "-|>")
          edge((1, 1), (2, 1), (2, 0.5), "-|>")
          edge((2, 0.5), (3.5, 0.5), "-|>", label: $epsilon(t, theta)$, label-side: center, label-sep: 0.5em)

          // Minimization
          node((3.5, 0.5), stroke: none, shape: "rect")[$min_theta cal(J)_N$]
        },
      ),
      caption: [PEM identification loop: the prediction error $epsilon(t, theta) = y(t) - hat(y)(t|t-1, theta)$ is formed and minimized over $theta$.],
    )
  }
]

#theorem(title: "Asymptotic PEM properties")[
  Under standard assumptions (stationary data, persistent excitation, true system in the model class):

  As $N -> infinity$:
  $ hat(theta)_N ->^p theta^* $
  where $theta^*$ minimizes the theoretical (expected) cost $overline(cal(J))(theta) = EE[epsilon(t, theta)^2]$.

  If the true system belongs to the model class, then $theta^* = theta_0$ (true parameters).
]

#remark(title: "Innovation form")[
  When the true system belongs to the model class, the prediction error becomes white noise:
  $ epsilon(t, theta_0) = e(t) tilde WN(0, lambda^2) $
  This is the basis for model validation (testing whiteness of residuals).
]

=== LS properties for ARX models

#theorem(title: "Unbiasedness of LS estimate")[
  If the true system is an ARX model with the correct orders:
  $ EE[hat(theta)_N] = theta_0 $
  The LS estimate is *unbiased*.
]

#theorem(title: "Variance of LS estimate")[
  $ "Var"[hat(theta)_N] = lambda^2 [sum_(t=1)^N phi(t) phi(t)^TT]^(-1) $
  The variance decreases as $1/N$ with more data.
]

#remark(title: "Noise variance estimate")[
  $ hat(lambda)^2 = cal(J)_N(hat(theta)_N) = 1/N sum_(t=1)^N epsilon(t, hat(theta)_N)^2 $
  is a consistent estimate of $lambda^2$ (since $cal(J)(hat(theta)) = lambda^2$ for the optimal model).
]

=== Persistent excitation

#definition(title: "Persistent Excitation (PE)")[
  An input signal $u(t)$ is *persistently exciting of order $n$* if the information matrix
  $ R_N = 1/N sum_(t=1)^N phi(t) phi(t)^TT $
  is positive definite for sufficiently large $N$.

  Equivalently, the Toeplitz matrix of input auto-covariances must be positive definite.
]

#remark[
  - A *white noise* input is PE of any order
  - A *sinusoidal* input at frequency $omega_0$ is PE of order 2 only
  - A *PRBS* (Pseudo-Random Binary Sequence) is approximately PE of high order
  - PE is necessary for the information matrix to be invertible, hence for identifiability
]

#definition(title: "Identifiability conditions")[
  For the LS estimate $hat(theta)_N$ to be unique and well-defined:
  + *Structural identifiability*: the model class $cal(M)(theta)$ is such that different $theta$ values give different models (no over-parameterization)
  + *Experimental identifiability*: the input is persistently exciting of sufficient order (the information matrix is non-singular)
]

== Identification of ARX models

$ cal(M)(theta): quad y(t) = B(z) / A(z) u(t-d) + bold(1) / A(z) e(t), quad e(t) ~ WN(0, lambda^2) $

#note-box[
  The ARX model is a special case of the ARMAX model with $C(z) = 1$.
]

Let's consider $A(z)$ and $B(z)$ in monic form:
$ A(z) = 1 + a_1 z^(-1) + a_2 z^(-2) + dots + a_n z^(-m) $
$ B(z) = b_0 + b_1 z^(-1) + b_2 z^(-2) + dots + b_p-1 z^(-p+1) $

Then the parameter vector is
$ theta = mat(a_1, a_2, dots, a_m, b_0, b_1, dots, b_(p-1))^T in Theta subset.eq RR^(m + p) $

Let's call $m_theta$ the size of the parameter vector $theta$. In this case,
$ m_theta = m + p $

We can rewrite the model as
$
  cal(M)(theta): quad & A(z)y(t) = B(z)u(t-d) + e(t) \
                      & bold(y(t) - y(t)) + A(z)y(t) = B(z)u(t-d) + e(t) \
                      & y(t) = (1- A(z))y(t) + B(z)u(t-d) + e(t)
$

$
  y(t) = &underbrace((a_1 z^(-1) + dots + a_m z^(-m))y(t) + (b_0 + b_1 z^(-1) + dots + b_(p-1) z^(-p+1))u(t-d), "Available at time" t-1"," \ "predictable") + &underbrace(e(t), "Not available," \ "prediction error")
$

Isolating the predictable part of the model, we obtain the *predictor* of the class of models $hat(cal(M))(theta)$:
$ hat(y)(t|t-1) = (a_1 z^(-1) + dots + a_m z^(-m))y(t) + (b_0 + b_1 z^(-1) + dots + b_(p-1) z^(-p+1))u(t-d)) $

We can use this predictor to define our optimization task in the *cost function*:

$ cal(J)_N (theta) = 1 / N sum_(t=1)^N (y(t+1) - hat(y)(t|t-1, theta))^2 $

#definition(title: "Regressor")[
  Given the class of predictor $hat(cal(M))(theta)$, it's called *regressor* the vector:
  $ phi(t) = [y(t-1), y(t-2), dots, u(t-d), dots, u(t-d-p+1)] $ namely, the vector of data points we're going to multiply by the chosen parameters to compute our our predictor.
]
Using that structure, we can rewrite the predictor in a compact vectorial form:
$ hat(y)(t|t-1) = theta^T phi(t) $

That reflects on the cost function:

$
  cal(J)_N (theta) & = 1 / N sum_(t=1)^N (y(t+1) - hat(y)(t|t-1, theta))^2 \
                   & = 1 / N sum_(t=1)^N (y(t+1) - theta^T phi(t) )^2 \
                   & = 1 / N sum_(t=1)^N (y(t+1) - phi(t)^T theta )^2
$

#remark[
  Since $hat(y)$ is linear in $theta$, $cal(J)_N$ is *quadratic* in $theta$. Hence, the global minimum $hat(theta)_N$ is unique.\
  For example, with $m_theta = 2$, $cal(J)_N$ is a paraboloid in $RR^3$.
]

#figure(
  cetz.canvas({
    import cetz.draw: *

    // Draw a simple 2D contour-style representation of a paraboloid
    // (isometric projection of level curves)
    let cx = 0
    let cy = 0

    // Axes
    line((-3.5, -1), (3.5, -1), stroke: gray + 0.5pt, mark: (end: ">"))
    content((3.8, -1), $theta_1$)
    line((-0.0, -2.5), (0.0, 2.5), stroke: gray + 0.5pt, mark: (end: ">"))
    content((0.4, 2.6), $theta_2$)

    // Elliptical contour lines (level sets of the paraboloid)
    for r in (0.4, 0.9, 1.4, 2.0, 2.7) {
      let rx = r * 1.3
      let ry = r * 0.8
      circle(
        (cx, cy),
        radius: (rx, ry),
        stroke: if r == 0.4 { blue + 1pt } else { blue.lighten(30%) + 0.6pt },
        fill: if r == 0.4 { blue.lighten(90%) } else { none },
      )
    }

    // Optimal point
    circle((cx, cy), radius: 0.06, fill: red)
    content((0.6, 0.3), text(fill: red, size: 9pt, $hat(theta)_N$))

    // Gradient arrow from some point to center
    let p = (1.8, 1.0)
    circle(p, radius: 0.05, fill: black)
    line(p, (0.4, 0.22), stroke: (paint: red, thickness: 1pt), mark: (end: ">"))
    content((2.3, 1.0), text(size: 8pt, $-nabla cal(J)_N$))

    // Labels
    content((2.5, -2.0), text(size: 8pt)[Level curves of $cal(J)_N (theta)$])
  }),
  caption: [For $m_theta = 2$, $cal(J)_N$ is a paraboloid. The contour lines are ellipses and the gradient points toward the unique global minimum $hat(theta)_N$.],
)

To find the global optimal $hat(theta)_N$ in general, the following two conditions must hold:
- $hat(theta)_N$ is a *stationary point* of $cal(J)_N$, therefore the gradient must be null:

$ evaluated(pdv(cal(J)_N (theta), theta))_(theta = hat(theta)_N) = 0 $



- $hat(theta)_N$ is a *minimum point* of $cal(J)_N$, therefore the hessian matrix must be positive definite:

$ evaluated(pdv(cal(J)_N (theta), theta, 2))_(theta = hat(theta)_N) succ 0 $

Let us analyze the gradient:

$
  pdv(cal(J)_N (theta), theta) &= vec(pdv(cal(J)_N (theta), a_1), dots, pdv(cal(J)_N (theta), b_(p-1))) \
  &= dv(, theta)[cal(J)_N (theta)] \
  &= dv(, theta)[1 / N sum_(t=1)^N (y(t) - phi(t)^TT theta )^2]\
  // invert summation and derivative
  &= 1 / N sum_(t=1)^N dv(, theta)[(y(t) - phi(t)^TT theta )^2]\
  &= 1 / N sum_(t=1)^N 2(y(t) - phi(t)^TT theta ) underbrace((dv(, theta)[y(t) - phi(t)^TT theta ]), bold(-phi(t)) \ "Not transposed, to make" \ "the gradient a column vector")\
  &= - 2 / N sum_(t=1)^N phi(t)(y(t) - phi(t)^TT theta )\
$

Now, let's set the gradient to zero and find the optimal $hat(theta)_N$:

$
  evaluated(pdv(cal(J)_N (theta), theta))_(theta = hat(theta)_N) &= 0\
  - 2 / N sum_(t=1)^N phi(t)(y(t) - phi(t)^TT hat(theta)_N ) &= 0\
  cancel(- 2 / N) sum_(t=1)^N phi(t)y(t) &= cancel(- 2 / N) sum_(t=1)^N phi(t)phi(t)^TT hat(theta)_N\
  [sum_(t=1)^N phi(t)phi(t)^TT] hat(theta)_N &= sum_(t=1)^N phi(t)y(t)\
$

If $sum_(t=1)^N phi(t)phi(t)^TT$, called *information matrix*, is non-singular, then we can solve for $hat(theta)_N$:

$
  hat(theta)_N & = [sum_(t=1)^N phi(t)phi(t)^TT]^(-1) sum_(t=1)^N phi(t)y(t) \
$

This result is called *ordinary least squares formulas*. It's the only case of explicit solution for identification problems adn depends on the data only.

Let's analyze the hessian matrix to check if the solution is a minimum:

$
  // start from the gradient and then derive it (divide in two summations)
  pdv(cal(J)_N (theta), theta, 2) &= dv(, theta)[pdv(cal(J)_N (theta), theta)] \
  &= dv(, theta)[underbrace(- 2 / N sum_(t=1)^N phi(t)y(t), "Does not depends on" theta) + 2 / N sum_(t=1)^N phi(t) phi(t)^TT theta] \
  &= 0 + 2 / N sum_(t=1)^N phi(t)phi(t)^TT quad forall theta\
$

We need to check if that result is positive definite.

#remark[
  A quadratic matrix $M$ is *positive definite* if
  $ x^T M x > 0 quad forall x != 0 $
]

So, we need to verify that:

$
  x^T pdv(cal(J)_N (theta), theta, 2) x &> 0 \
  x^T 2 / N sum_(t=1)^N phi(t)phi(t)^TT x &= 2 / N sum_(t=1)^N (underbrace(x^T phi(t), "scalar") quad underbrace(phi(t)^TT x, "scalar"))\
  &= 2 / N sum_(t=1)^N (phi(t)^TT x)^T phi(t)^TT x\
  &= 2 / N sum_(t=1)^N (phi(t)^TT x)^2 >= 0 quad forall x != 0\
$

We actually proved that
$ x^T pdv(cal(J)_N (theta), theta, 2) x >= 0 quad forall x != 0 $
that means that two cases can occour:
- if $pdv(cal(J)_N (theta), theta, 2) succ 0$ the hessian matrix is positive definite and $hat(theta)_N$ is a minimum;
- if $pdv(cal(J)_N (theta), theta, 2) = 0$ (degenerate case) the hessian matrix is singular and we have an infinite number of minimum points.

There could be two reasons behind the degenerate case:
- *experimantal identifiability issue*: data are not representive of the full underlying physical process; This means the information matrix $phi(t)phi(t)^T$ is not all informative and some data is the linear combination of some other data.
- *structural identifiability issue*: the selected model class is too large, there is more than one model that can fit the data.



== Identification of ARMA(X) models

$ cal(M)(theta): quad y(t) = B(z) / A(z) u(t-d) + C(z) / A(z) e(t), quad e(t) ~ WN(0, lambda^2) $
With
- $C(z) = 1+ c_1 z^(-1) + c_2 z^(-2) + dots + c_n z^(-n)$
- $A(z) = 1 - ( a_1 z^(-1) + a_2 z^(-2) + dots + a_m z^(-m))$
- $B(z) = b_0 + b_1 z^(-1) + b_2 z^(-2) + dots + b_p-1 z^(-p+1)$

Then the parameter vector is
$ theta = mat(a_1, a_2, dots, a_m, b_0, b_1, dots, b_(p-1), c_1, c_2, dots, c_n)^T in Theta subset.eq RR^(m + p + n) $
Let's call $n_theta$ the size of the parameter vector $theta$.

We have the same loss function:
$
  J = 1 / N sum_(t=1)^N (y(t) - y(t-1|t, theta))^2
$

We're interested in the one step ahead predictor $k=1$.


#let longdiv(all-columns, ..cells) = {
  let longdiv = grid
  let cols = if type(all-columns) == array {
    all-columns.len()
  } else if type(all-columns) == int {
    all-columns
  } else {
    1
  }
  set grid(
    columns: cols,
    inset: 5pt,
    align: right,
    stroke: (x, y) => (
      // Add left stroke to the last column
      left: if x == cols - 1 { black },
      bottom: if (
        // Add bottom stroke to the top right cell
        y == 0 and x == cols - 1 // Add bottom stroke every two rows (calc.odd check),
          // but for one less column each time
          or x < cols - 1 and calc.odd(y) and x + 1 >= y / 2
      ) {
        black
      },
    ),
  )
  grid(..cells)
}

$
  longdiv(
    #2, C(z), A(z),
    -A(z), E(z) = 1,
    C(z)- A(z)
  )
$

Let's call $F(z)z^(-k) := C(z)- A(z)$ the rest of our long division.
Our predictor is therefore:
$ hat(y)(t|t-1) = (C(z)-A(z)) / C(z) y(t) +( B(z) E(z)) / C(z) u(t-d) $

The stationarity of the process gives us that, for a generic k
$ hat(y)(t + k |t) = (C(z)-A(z)) / C(z) y(t + k) +( B(z) E(z)) / C(z) u(t+k-d) $

Returning to $k=1$, we can compute the prediction error:
$
  epsilon(t|t-1, theta) = A(z) / C(z) y(t) - B(z) / C(z)u(t-d)
$

#note-box()[We cannot apply OLS since $hat(theta)$ appears in the denominator of the prediction appears in C(z).
  We cannot write $hat(y)(t|t-1, theta)$ as a linear function $phi(t)^T theta$.
  We have to use non linear function of $theta$ so to optimize the prediction error we're going to use an *heuristic / iterative numerical appproach*]

#note-box()[
  Iterative algorithms are guaranteed to converge to local minimima not global minima.
]

To update our parameters, there are various strategies. We could use *Newtons method* and find the tangent paraboloid around the current iteration of the estimate of , our parameter $Theta^((i))$, then take the minimum of the computed paraboloid as the new estimate of $J_N (theta)$

$
  nu (Theta) = J_N (Theta^((i))) + pdv(J_N (Theta), Theta)_(Theta = Theta^((i))) (Theta - Theta^((i))) + 1 / 2 (Theta - Theta^((i)))^T pdv(J_N (Theta), Theta, 2)_(Theta = Theta^((i))) (Theta - Theta^((i)))
$

=== Results
After computing the expressions of the Gradient and Hessian matrix, we can choose to update our parameters $theta^((i))$ always in the direction of the gradient, since it's going to give us the direction toward a stationary point of our loss function, and use the hessian matrix convex, hopefully or approximated definite positive to go towards a minima.

We can choose different strategies to change the magnitude of our jump:
- *Newton's rule*: to have the learning rate of the algorithm according to the curvature of the paraboloid.
$ theta^((i+1)) = theta^((i)) - ["hessian"]^(-1) dot "gradient" $
- *Quasi-Newton*: to have the learning rate of the algorithm according to the curvature of the paraboloid.
  $ theta^((i+1)) = theta^((i)) - ["positive terms of the hessian"]^(-1) dot "gradient" $
- *Gradient descent*: to have a scalar and fixed learning rate. Fastest and simplest method.
$ theta^((i+1)) = theta^((i)) - eta ["gradient"] $

=== Quasi-Newton derivation

The idea of Newton's method is to approximate $cal(J)_N (theta)$ with a quadratic (Taylor expansion around the current iterate $theta^((i))$), then take the minimum of the paraboloid as the next estimate.

#figure(
  cetz.canvas({
    import cetz.draw: *

    let w = 10
    let h = 5

    // Axes
    line((-0.5, 0), (w + 0.5, 0), stroke: gray + 0.5pt, mark: (end: ">"))
    content((w + 0.8, -0.2), $theta$)
    line((0.5, -0.3), (0.5, h + 0.3), stroke: gray + 0.5pt, mark: (end: ">"))

    // J_N(theta) -- non-quadratic curve with local structure
    let pts = ()
    for i in range(0, 81) {
      let x = i * w / 80
      let y = (
        0.3
          + 0.8 * calc.pow(x - 4, 2) / 10
          + 0.5 * calc.sin(x * 1.2) * calc.exp(-0.15 * calc.pow(x - 2, 2))
          + 0.2 * calc.exp(-0.5 * calc.pow(x - 2, 2))
      )
      pts.push((x, y))
    }
    for i in range(1, pts.len()) {
      line(pts.at(i - 1), pts.at(i), stroke: black + 1.2pt)
    }
    content((w - 0.5, h - 0.8), $cal(J)_N (theta)$)

    // Tangent parabola at theta^(i) = 7
    let ti = 7.0
    let yi = (
      0.3
        + 0.8 * calc.pow(ti - 4, 2) / 10
        + 0.5 * calc.sin(ti * 1.2) * calc.exp(-0.15 * calc.pow(ti - 2, 2))
        + 0.2 * calc.exp(-0.5 * calc.pow(ti - 2, 2))
    )
    let parab-pts = ()
    for i in range(0, 41) {
      let x = ti - 3 + i * 6 / 40
      let dx = x - ti
      let y = yi - 0.8 * dx + 0.25 * dx * dx // tangent parabola
      if y > 0 and x > 0 and x < w {
        parab-pts.push((x, y))
      }
    }
    for i in range(1, parab-pts.len()) {
      line(parab-pts.at(i - 1), parab-pts.at(i), stroke: blue + 1pt)
    }
    content((8.5, h - 0.3), text(fill: blue, $nu_i (theta)$))

    // Iterates
    let theta-i = ti
    let theta-i1 = 5.4
    let theta-i2 = 4.5

    line((theta-i, 0), (theta-i, yi), stroke: (dash: "dashed", paint: blue))
    circle((theta-i, 0), radius: 0.06, fill: blue)
    content((theta-i, -0.35), text(fill: blue, size: 8pt, $theta^((i))$))

    line((theta-i1, 0), (theta-i1, 0.5), stroke: (dash: "dashed", paint: red))
    circle((theta-i1, 0), radius: 0.06, fill: red)
    content((theta-i1, -0.35), text(fill: red, size: 8pt, $theta^((i+1))$))

    line((theta-i2, 0), (theta-i2, 0.35), stroke: (dash: "dashed", paint: orange))
    circle((theta-i2, 0), radius: 0.06, fill: orange)
    content((theta-i2, -0.35), text(fill: orange, size: 8pt, $theta^((i+2))$))
  }),
  caption: [Newton's method: at each iterate $theta^((i))$, a tangent parabola $nu_i (theta)$ approximates $cal(J)_N (theta)$. The minimum of the parabola gives $theta^((i+1))$, converging toward the local minimum.],
)

$
  nu (theta) = cal(J)_N (theta^((i))) + pdv(cal(J)_N, theta)|_(theta^((i))) (theta - theta^((i))) + 1 / 2 (theta - theta^((i)))^TT pdv(cal(J)_N, theta, 2)|_(theta^((i))) (theta - theta^((i)))
$

Setting $pdv(nu, theta) = 0$ gives the update rule:
$ theta^((i+1)) = theta^((i)) - [pdv(cal(J)_N, theta, 2)|_(theta^((i)))]^(-1) pdv(cal(J)_N, theta)|_(theta^((i))) $

The *gradient* of the cost function is:
$ pdv(cal(J)_N (theta), theta) = 2 / N sum_(t=1)^N epsilon(t|t-1, theta) dot pdv(epsilon(t|t-1, theta), theta) $

The *Hessian* is (applying the product rule to the gradient):
$
  pdv(cal(J)_N (theta), theta, 2) = 2 / N sum_(t=1)^N pdv(epsilon, theta) (pdv(epsilon, theta))^TT + 2 / N sum_(t=1)^N epsilon(t) pdv(epsilon, theta, 2)
$

#remark(title: "Quasi-Newton approximation")[
  The second-order term involving $pdv(epsilon, theta, 2)$ is neglected:
  $ pdv(cal(J)_N, theta, 2) approx 2 / N sum_(t=1)^N pdv(epsilon, theta) (pdv(epsilon, theta))^TT $

  This *approximate Hessian* is always positive semi-definite (it is a sum of outer products), which guarantees the update direction points toward a minimum --- we always get a "bowl-shaped" paraboloid.
]

The Quasi-Newton update rule then becomes:
$
  theta^((i+1)) = theta^((i)) - [sum_(t=1)^N pdv(epsilon(t), theta) (pdv(epsilon(t), theta))^TT]^(-1) [sum_(t=1)^N epsilon(t) pdv(epsilon(t), theta)]
$

=== Computing $pdv(epsilon, theta)$ with auxiliary signals

The gradient $pdv(epsilon(t|t-1, theta), theta)$ is a vector of partial derivatives with respect to each parameter in $theta = mat(a_1, dots, a_m, b_0, dots, b_(p-1), c_1, dots, c_n)^TT$.

Three *auxiliary signals* are defined to simplify the computation:
$
  alpha(t) = -1 / C(z) y(t), quad
  beta(t) = -1 / C(z) u(t), quad
  gamma(t) = -1 / C(z) epsilon(t|t-1, theta)
$

#remark[
  These are all time signals: the derivative of a signal with respect to a parameter is itself a signal. Note that $gamma(t)$ depends on the current value of $theta$.
]

Using them, the gradient of the prediction error takes the compact form:
$
  pdv(epsilon(t|t-1, theta), theta) = vec(alpha(t-1), dots.v, alpha(t-m), beta(t-1), dots.v, beta(t-p+1), gamma(t-1), dots.v, gamma(t-n))
$

#figure(
  {
    import fletcher: diagram, edge, node
    diagram(
      spacing: (2em, 1.5em),
      node-stroke: 1pt,
      node-shape: "rect",
      {
        // Inputs
        node((-2, -0.5), stroke: none)[$u(t)$]
        node((-2, 0.5), stroke: none)[$y(t)$]
        node((-2, 1.5), stroke: none)[$epsilon$]

        // -1/C(z) filter block
        node((0, 0.5))[$display(-1 / C(z))$]
        edge((-2, -0.5), (-0.5, -0.5), (0, 0.5), "-|>")
        edge((-2, 0.5), (0, 0.5), "-|>")

        // epsilon filter
        node((0, 1.5))[$display(-1 / C(z))$]
        edge((-2, 1.5), (0, 1.5), "-|>")

        // Auxiliary signals
        node((2, -0.5), stroke: none)[$alpha(t)$]
        node((2, 0.5), stroke: none)[$beta(t)$]
        node((2, 1.5), stroke: none)[$gamma(t)$]
        edge((0, 0.5), (2, -0.5), "-|>")
        edge((0, 0.5), (2, 0.5), "-|>")
        edge((0, 1.5), (2, 1.5), "-|>")

        // Gradient vector
        node((4, 0.5))[$display(pdv(epsilon, theta)) = vec(alpha(t-1), dots.v, gamma(t-n))$]
        edge((2, -0.5), (4, 0.5), "-|>")
        edge((2, 0.5), (4, 0.5), "-|>")
        edge((2, 1.5), (4, 0.5), "-|>")

        // Update rule
        node((6, 0.5))[UPDATE \ RULE]
        edge((4, 0.5), (6, 0.5), "-|>", label: $pdv(epsilon, theta)$, label-side: center, label-sep: 0.8em)

        // theta^(i) input
        node((6, -0.5), stroke: none)[$theta^((i))$]
        edge((6, -0.5), (6, 0.5), "-|>")

        // theta^(i+1) output
        node((8, 0.5), stroke: none)[$theta^((i+1))$]
        edge((6, 0.5), (8, 0.5), "-|>")
      },
    )
  },
  caption: [Processing scheme for the Newton / Quasi-Newton update law. The signals $u(t), y(t), epsilon$ are filtered through $-1 slash C(z)$ to produce auxiliary signals $alpha(t), beta(t), gamma(t)$, which form $pdv(epsilon, theta)$ and drive the parameter update.],
)

=== ARMAX PEM details

#definition(title: "Pseudo-regressor for ARMAX")[
  The ARMAX 1-step predictor can be written as:
  $ hat(y)(t|t-1, theta) = psi(t, theta)^TT theta $
  where $psi(t, theta)$ is the *pseudo-regressor*:
  $
    psi(t, theta) = vec(y(t-1), dots, y(t-m), u(t-d), dots, u(t-d-p+1), epsilon(t-1, theta), dots, epsilon(t-n, theta))
  $

  Note: $psi$ depends on $theta$ through the past prediction errors $epsilon(t-i, theta)$.
]

#remark(title: "Initialization strategy")[
  The iterative algorithm requires an initial guess $theta^((0))$. Common strategies:
  - Start from the LS estimate of a simpler ARX model (ignoring $C(z)$)
  - Use instrumental variables
  - Try multiple random initializations and keep the best

  The algorithm should be run until convergence: $||theta^((i+1)) - theta^((i))|| < "tolerance"$.
]

== Uncertainty of parameter estimates

#theorem(title: "Variance of PEM estimates")[
  For ARX models, the covariance matrix of $hat(theta)_N$ is:
  $ "Var"[hat(theta)_N] = lambda^2 [sum_(t=1)^N phi(t) phi(t)^TT]^(-1) $

  Asymptotically (as $N -> infinity$):
  $ sqrt(N)(hat(theta)_N - theta_0) -> cal(N)(0, P) $
  where $P$ depends on the model structure.
]

#remark(title: "Confidence intervals")[
  For each parameter $theta_i$, the $(1-alpha)$ confidence interval is:
  $ hat(theta)_(N,i) plus.minus z_(alpha\/2) sqrt("Var"[hat(theta)_(N,i)]) $
  where $z_(alpha\/2)$ is the quantile of the standard normal distribution.
]

== Bias analysis

#caution-box(title: "Wrong model class")[
  If the true system does not belong to the selected model class, the LS estimate is generally *biased*:
  $ EE[hat(theta)_N] != theta_0 $

  Common sources of bias:
  - Wrong model orders ($m, n, p$)
  - Missing dynamics (e.g., using ARX when the true system is ARMAX)
  - Correlated noise (using ARX when noise is colored by $C(z) != 1$)
]

#remark(title: "Bias in ARX with colored noise")[
  If the true system is ARMAX but we fit an ARX model:
  $ y(t) = B(z)/A(z) u(t-d) + C(z)/A(z) e(t) $
  and we use $hat(y)(t|t-1) = phi(t)^TT theta$ (ARX predictor), the LS estimate will be *biased* because the noise term $C(z)/A(z) e(t)$ is not white and correlates with the regressor.
]

== Dummy identification exercise

#example(title: "Identifying an MA(1) process")[
  Suppose the true system is:
  $ y(t) = e(t) + 1 / 2 e(t-1), quad e(t) tilde WN(0, 1) $

  and we want to identify it with a model $cal(M)$:
  $ y(t) = eta(t) + b eta(t-1), quad eta(t) tilde WN(0, lambda^2) $

  The prediction error for 1-step prediction is:
  $ epsilon(t|t-1) = A(z) / C(z) y(t) = 1 / (1 + b z^(-1)) y(t) = (1 + 1 / 2 z^(-1)) / (1 + b z^(-1)) e(t) $

  For $epsilon(t)$ to be white noise (equal to $e(t)$ with $lambda^(2*) = 1$), we need the transfer function to simplify to 1, which requires:
  $ b^* = 1 / 2 $

  With $b^* = 1/2$, the prediction error becomes $epsilon(t) = e(t) tilde WN(0, 1)$, confirming the model matches the true system.

  The key insight is to recognize this by inspection: the model structure directly matches the system --- no optimization is needed.
]

