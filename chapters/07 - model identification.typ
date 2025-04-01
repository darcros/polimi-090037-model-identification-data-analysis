#import "/prelude.typ": *

= Model identification

#problem(title: "Identification problem")[
  Given a system $cal(S)$, we collect samples of input $u(t)$ and output $y(t)$ in order to derive a model of the system.

  #{
    import fletcher: diagram, node, edge
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
    mat(delim: #none, align: #left, row-gap: #1em, column-gap: #0.5em,
      y(t) = #h(1em)
      , , #hide[+] a_1 y(y-1), + a_2 y(t-2), dots, + a_m y(t-m);
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
    // FIXME: the notation for the model class is a bit weird, but it is what the prof. wrote so...
  ]
]

=== Choice of the identification criterion <model-ident:step-3>

We use a _predictive approach_. We assume our model is good if the prediction error given by the predictor we have chosen to estimate our future value of our time series is low.// TODO: explain a little bit better

Our ideal objective would be
$ cal(J)(theta) = EE[(y(t+1) - hat(y)(t+1|t, theta))^2] $
but that is not computable, because we only have *one* realization.

So we use a sample-based version instead
$ cal(J)(theta) = 1 / N sum_(t=1)^N (y(t+1) - hat(y)(t|t-1, theta))^2 $

// TODO: explain better
#remark[
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

=== Model validation

// TODO: explain better the assumptions, why they could be wrong, what happens if they are
During the process we made a few assumptions that must be verified in order to validate the quality of the model:
- the system $cal(S)$ lies within the model class $cal(M)_theta$, meaning that the model has the same number of parameters, in the same structure and elevated to the correct magnitude.
- $m$, $n$, $p$ are specified a priori to data have been "correctly" selected (they are informative)

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
  cal(M)(theta): quad& A(z)y(t) = B(z)u(t-d) + e(t) \
  & bold(y(t) - y(t)) + A(z)y(t) = B(z)u(t-d) + e(t) \
  & y(t) = (1- A(z))y(t) + B(z)u(t-d) + e(t)
$

$
  y(t) = &underbrace((a_1 z^(-1) + dots + a_m z^(-m))y(t) + (b_0 + b_1 z^(-1) + dots + b_(p-1) z^(-p+1))u(t-d), "Available at time" t-1"," \ "predictable") + &underbrace(e(t), "Not available," \ "prediction error")
$

Isolating the predictable part of the model, we obtain the *predictor* of the class of models $ hat(cal(M))(theta)$:
$ hat(y)(t|t-1) = (a_1 z^(-1) + dots + a_m z^(-m))y(t) + (b_0 + b_1 z^(-1) + dots + b_(p-1) z^(-p+1))u(t-d)) $

We can use this predictor to define our optimization task in the *cost function*:

$ cal(J)_N (theta) = 1 / N sum_(t=1)^N (y(t+1) - hat(y)(t|t-1, theta))^2 $

#definition(title: "Regressor")[
  Given the class of predictor $hat(cal(M))(theta)$, it's called *regressor* the vector:
  $
    phi(t) = [y(t-1), y(t-2), dots, u(t-d), dots, u(t-d-p+1)]
  $ namely, the vector of data points we're going to multiply by the chosen parameters to compute our our predictor.
]
Using that structure, we can rewrite the predictor in a compact vectorial form:
$ hat(y)(t|t-1) = theta^T phi(t) $

That reflects on the cost function:

$
  cal(J)_N (theta) &= 1 / N sum_(t=1)^N (y(t+1) - hat(y)(t|t-1, theta))^2 \
  &= 1 / N sum_(t=1)^N (y(t+1) - theta^T phi(t) )^2 \
  &= 1 / N sum_(t=1)^N (y(t+1) - phi(t)^T theta )^2
$

#remark[
  Since $hat(y)$ is linear in $theta$, $cal(J)_N$ is *quadratic* in $theta$. Hence, the global minimum $hat(theta)_N$ is unique.\
  For example, with $m_theta = 2$, $cal(J)_N$ is a paraboloid in $RR^3$.
]

// TODO: plot Jn in R^3 as a paraboloid and show the global optimal

To find the global optimal $hat(theta)_N$ in general, the following two conditions must hold:
- $hat(theta)_N$ is a *stationary point* of $cal(J)_N$, therefore the gradient must be null:

$ evaluated(pdv(cal(J)_N (theta), theta))_(theta = hat(theta)_N) = 0 $



- $hat(theta)_N$ is a *minimum point* of $cal(J)_N$, therefore the hessian matrix must be positive definite:

$ evaluated(pdv(cal(J)_N (theta), theta, 2))_(theta = hat(theta)_N) succ 0 $

Let us analyze the gradient:

$
  pdv(cal(J)_N (theta), theta) &= vec(pdv(cal(J)_N (theta), a_1), dots, pdv(cal(J)_N (theta), b_(p-1))) \
  &= dv(,theta)[cal(J)_N (theta)] \
  &= dv(,theta)[1 / N sum_(t=1)^N (y(t) - phi(t)^TT theta )^2]\
  // invert summation and derivative
  &= 1 / N sum_(t=1)^N dv(,theta)[(y(t) - phi(t)^TT theta )^2]\
  &= 1 / N sum_(t=1)^N 2(y(t) - phi(t)^TT theta ) underbrace((dv(,theta)[y(t) - phi(t)^TT theta ]), bold(-phi(t)) \ "Not transposed, to make" \ "the gradient a column vector")\
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
  hat(theta)_N &= [sum_(t=1)^N phi(t)phi(t)^TT]^(-1) sum_(t=1)^N phi(t)y(t)\
$

This result is called *ordinary least squares formulas*. It's the only case of explicit solution for identification problems adn depends on the data only.

Let's analyze the hessian matrix to check if the solution is a minimum:

$

// start from the gradient and then derive it (divide in two summations)
pdv(cal(J)_N (theta), theta, 2) &= dv(,theta)[pdv(cal(J)_N (theta), theta)] \
&= dv(,theta)[underbrace(- 2 / N sum_(t=1)^N phi(t)y(t), "Does not depends on" theta) + 2 / N sum_(t=1)^N phi(t) phi(t)^TT theta] \

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
- if $pdv(cal(J)_N (theta), theta, 2) succ 0 $ the hessian matrix is positive definite and $hat(theta)_N$ is a minimum;
- if $pdv(cal(J)_N (theta), theta, 2) = 0 $ (degenerate case) the hessian matrix is singular and we have an infinite number of minimum points.

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
  #2,  C(z), A(z),
  -A(z), E(z) = 1,
  C(z)- A(z)
)
$

Let's call $F(z)z^(-k) := C(z)- A(z)$ the rest of our long division.
Our predictor is therefore:
$ hat(y)(t|t-1) = (C(z)-A(z)) / C(z) y(t) + B(z) / C(z) u(t-d) $

We can compute the prediction error:
$
  epsilon(t|t-1, theta) = A(z) / C(z) y(t) - B(z) / C(z)u(t-d)
$

#note-box()[We cannot apply OLS since $hat(theta)$ appears in the denominator of the prediction appears in C(z), so in the denominator.
  We cannot write $hat(y)(t|t-1, theta)$ as a linear function $phi(t)^T theta$.
  We have to use non linear function of $theta$ so to optimize the prediction error we're going to use an *heuristic / iterative numerical appproach*]

#note-box()[
  Iterative algorithms are guaranteed to converge to local minimima not global minima.
]

To update our parameters, there are various strategies. We could use *Newtons method* and find the tangent paraboloid around the current iteration of the estimate of , our parameter $Theta^((i))$, then take the minimum of the computed paraboloid as the new estimate of $J_N (theta)$

//TODO maybe add figure of parabolic gradient descent

$
  nu (Theta) = J_N (Theta^((i))) + pdv(J_N (Theta), Theta)_(Theta = Theta^((i))) (Theta - Theta^((i))) + 1 / 2 (Theta - Theta^((i)))^T pdv(J_N (Theta), Theta, 2)_(Theta = Theta^((i))) (Theta - Theta^((i)))
$

//TODO add the computation notes if we care

=== Results
After computing the expressions of the Gradient and Hessian matrix, we can choose to update our parameters $theta^((i))$ always in the direction of the gradient, since it's going to give us the direction toward a stationary point of our loss function, and use the hessian matrix convex, hopefully or approximated definite positive to go towards a minima.

We can choose different strategies to change the magnitude of our jump:
- *Newton's rule*: to have the learning rate of the algorithm according to the curvature of the paraboloid.
$ theta^((i+1)) = theta^((i)) - ["hessian"]^(-1) dot "gradient" $
- *Quasi-Newton*: to have the learning rate of the algorithm according to the curvature of the paraboloid.
  $ theta^((i+1)) = theta^((i)) - ["positive terms of the hessian"]^(-1) dot "gradient" $
- *Gradient descent*: to have a scalar and fixed learning rate. Fastest and simplest method.
$ theta^((i+1)) = theta^((i)) - eta ["gradient"] $

