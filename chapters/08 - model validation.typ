#import "/prelude.typ": *

= Model validation
So far we have introduced an algorrithm to find the best parameters given a model. We're minimizing the empirical variance of the prediction error.
We need a quality assesment of the model. We need to check if the model is good enough to be used for prediction

$D_N = { (u(1), y(1)), (u(2), y(2)), \ldots, (u(N), y(N)) }$ our dataset of dimension $N$

$M_theta = {M(Theta), theta in Theta subset RR^(n_theta)}$ the model we want to validate. 

$ hat(Theta)_N = argmin_theta J_N (theta)$ with $J_N (theta) = 1/N sum_{i=1}^N (y(i) - M_theta(u(i)))^2 $.

//todo add figure to show different realizations of J_N generate different optimal parameters estimation

#remark()[
  $hat(theta)$ is going to change with the dataset since it depends on the realization of our loss function $J_N$ and therefore our process $y(t)$ and noise $e(t)$.
]

#theorem()[
  Under current assumptionss, as $N -> infinity$:
  $ J_N (hat(theta), s) -->_(N -> infinity) dash(J) (theta) = EE[epsilon (t|t-1, theta, s)^2] $ 
  Moreover, by letting 
  $ Delta = {theta^*, J(theta^* <= dash(J)(theta^*) forall theta} $
  be the set of global minima points of $ dash(J)(theta^*)$ we have
  $hat(theta_N) (s) -->_(N -> infinity) Delta$ with $PP(dot) = 1$
]

#corollary()[
  if $Delta= {theta^*}$ is a singleton we have that  $ hat(theta_N) (s) -->_(N -> infinity) theta^* "with "PP(dot) = 1 $
]

//todo add computations and diagrams about systems being the model class and delta being or not a singleton.

== Model order selection
Let's find the best dimension of $cal(M)_theta$ such that our system $cal(S) in cal(M)_theta$

// todo make this section actual code. 
$n := 1$ 

$"while" n<= n_max$

#h(0.7cm) $M_theta^((n)) = {M(theta), theta in Theta subset RR^(n_theta)}$

#h(0.7cm) $hat(theta)_N^((n)) = argmin_theta J_N^((n))(theta) $

#h(0.7cm) $J_N^((n))(theta) = 1/N sum_(i=1)^N (y(i) - M_theta^((n))(u(i)))^2 $


#figure(
)[
 

  #cetz.canvas({
    import suiji: gen-rng
      import cetz.draw: *
      import cetz-plot: *
      import "../util.typ": random-series
 
      let n_max =20
      let actual_model_oder=10

      let points = range(0,n_max).map(
        x => (x+3, 20/(1+x) +1 )
      ) 

      plot.plot(
        size: (5,5),
        axis-style: "school-book",
        x-label: [$n$ model order],
        x-max: n_max,
        y-max: n_max,
        x-tick-step: 50,
        y-tick-step: 50,
        y-label: [$J_N (theta)$],
             x-ticks: ((actual_model_oder, [$n_theta$]), 0),
        {
          plot.add(
            points, line: "spline"
          )

          plot.add-vline(actual_model_oder, style: (stroke:(dash:"dashed")))
            plot.add-fill-between(
            range(0, actual_model_oder+1).map(x => (x, 0)),
            range(0, actual_model_oder+1).map(x => (x, 20)),
            style: (fill: rgb("#ff00006c"), stroke: none), 
            label: ["underfitting"] // Highlight underfitting in red
            )
            plot.add-fill-between(
            range(actual_model_oder, actual_model_oder+11).map(x => (x, 0)),
            range(actual_model_oder, actual_model_oder+11).map(x => (x, 20)),
                        style: (fill: rgb("#0000ff46"), stroke: none),
            label: ["overfitting"]
            )

        }

      )
      


  }
    
  )
]

#note-box[We see that our loss function $J_N  (theta)$ is inversely proportional to the number of parameters in the model. The more parameters we have, the better we can fit the data. But this is not always a good thing since we may be fitting the noise of the data too. To avoid *overfitting*, we need to find a balance between the number of parameters and the goodness of fit.]


Three criteria for model order selections are
+ Whiteness test on the residuals (Anderson's test)
+ Cross validation
+ Identification of the model order penalties.


==== Whiteness test on the residuals
For a large enough $N$, $epsilon(t)$ is a white noise process, therefore we compute the covariance function or the spectrum of the prediction error to see if it has the same shape of a white noise's.

==== Cross validation with k fold
We split the dataset into k folds. We train the model on k-1 folds and test it on the last fold. We repeat this process k times, each time using a different fold as the test set. We then compute the average error over all k folds.

This is not adequeate for time series data since the data is not independent and we would be splicing together different time windows and we would be inserting in the dataset fake temporal correlations.
 We need to use a different approach for time series data.

==== Cross validation with model order penalties
Instead of minimizing J_N(theta), we can minimize a penalized version of it:

#figure(
)[
 

  #cetz.canvas({
    import suiji: gen-rng
      import cetz.draw: *
      import cetz-plot: *
      import "../util.typ": random-series
 
      let n_max =20
      let actual_model_oder=10

      let points = range(0,n_max).map(
        x => (x+3, 20/(1+x) +1 )
      ) 

      let AIC_points = range(0,n_max).map(
        x => (x+3, calc.ln((50/(1+x) +10 ) + 2*x/(n_max))
      ))
      let FPE_points = range(0,n_max).map(
        x => (x+3, 20/(1+x) +1 + (n_max+x)/(n_max - x))
      )
      let MDL_points = range(0,n_max).map(
        x => (x+3, calc.ln(100/(1+x) +10) +x*calc.ln(x+1)/(2*n_max))
      )

      plot.plot(
        size: (5,5),
        axis-style: "school-book",
        x-label: [$n$ model order],
        x-max: n_max,
        y-max: n_max,
        x-tick-step: 50,
        y-tick-step: 50,
        y-label: [$J_N (theta)$],
             x-ticks: ((actual_model_oder, [$n_theta$]), 0),
        {
          plot.add(
            points, line: "spline", label: ["J_N(theta)"]
          )

          plot.add-vline(actual_model_oder, style: (stroke:(dash:"dashed")))
          plot.add(
            AIC_points, line: "spline",  label: ["AIC"]
          )
          plot.add(
            FPE_points, line: "spline",  label: ["FPE"]
          )
          plot.add(
            MDL_points, line: "spline",  label: ["MDL"]
          )
          

        }

      )
      
  }
    
  )
]

+ *FPE*: Final Prediction Error $"FPE" = (N+ n)/(N-n) J_N(theta)$
+ *AIC*: Akaike Information Criterion $"AC"(n) = ln(J_N(theta)) + 2n/(N)$
+ *MDL*: Minimum Description Length $"MDL"(n) = ln(J_N(theta)) + n ln(N)/(2N)$

#note-box[
  If $cal(S) in cal(M)_theta$, and $cal(M)_theta$ is in the set of *ARX* model, MDL is right.
  If $cal(S) in cal(M)_theta$, and $cal(M)_theta$ is in the set of *ARMAX* model we prefer slightly to overfit so we tend to use AIC.
]
