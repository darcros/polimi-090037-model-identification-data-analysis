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
  Consider the class of $"ARMAX"(m, n, k, p)$ models
  $
    mat(delim: #none, align: #left, row-gap: #1em, column-gap: #0.5em,
      y(t) = #h(1em)
      , , #hide[+] a_1 y(y-1), + a_2 y(t-2), dots, + a_m y(t-m);
      , +b_0 u(t-k), +b_1 u(t-k-1), + b_2 u(t-k-2), dots, + b_m u(t-k-p);
      , +c_0 e(t), +c_1 e(t-1), + c_2 e(t-2), dots, + c_m e(t-n);
    )
  $
  where $e(t) ~ "WN"(0, lambda^2)$.

  Then
  $ theta = mat(a_1, a_2, dots, b_0, b_1, b_2, dots, c_0, c_1, c_2, dots)^transposed in Theta subset.eq RR^(m + n + p) $

  #note-box[
    $theta$ is not sufficient to fully specify the model.
    We also need to know
    - $lambda^2$: the variance of the white noise
    - $k$: the pure delay of the system
    - m, n, p: the degrees of the $"AR"$, $"MA"$ and $"X"$ parts

    The first two are not critical in the sense that they can be obtained "for free" through the minimaztion of our identification criterion (see @model-ident:step-3).
    While the latter must be specified ahead a priori.
  ]

  #note-box[
    $Theta != RR^(m + n + p)$ because some values in $RR^(m + n + p)$ are not admissible.

    Consider $y(t) = a_1 y(t-1) + e(t) ~ "AR"(1) = "ARMAX"(1, 0, 0, 0)$

    In order for $y(t)$ to be static, we need $|a_1| < 1$ (see @thm:stationarity).

    So $Theta = {a : |a| < 1} subset RR^1$ and our model class would be
    $cal(M)_theta = {"AR"(1), theta = a, |a| < 1}$
    // FIXME: the notation for the model class is a bit weird, but it is what the prof. wrote so...
  ]
]

=== Choice of the identification criterion <model-ident:step-3>

We use the _predictive approach_. // TODO: explain a little bit better

Our ideal objective would be
$ cal(J)(theta) = EE[(y(t+1) - hat(y)(t+1|t, theta))^2] $
but that is not computable, because we only have *one* realization.

So we use a sample-based version instead
$ cal(J)(theta) = 1 / N sum_(t=1)^N (y(t+1) - hat(y)(t|t-1, theta))^2 $

// TODO: explain better
#remark[
  We "look" 1-step ahead (instead of $2, 3, dots, k$) because

  $ epsilon(t+k|t) = t(t+k) - hat(y)(t+k|t) $
  $ "var"[e(t)] <= "var"[epsilon(t+k|t)] <= "var"[y(t)] $

  if the model is exact and we computer the optial 1-step ahead error
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
- $cal(J)(theta)$ is quadratic, this is the case for $"AR"$ and $"ARX"$ models
- $cal(J)(theta)$ is non-quadratic, this is the case for $"ARMA"$ and $"ARMAX"$ models

=== Model validation

// TODO: explain better the assumptions, why they could be wrong, what happens if they are
During the process we made a few assumptions that must be verified in order to validate the quality of the model:
- the system $cal(S)$ lies within the model class $cal(M)_theta$
- $m$, $n$, $p$ are specified a priori
- data have been "correctly" selected (they are informative)
