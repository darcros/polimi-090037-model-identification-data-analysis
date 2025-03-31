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

//todo insert graphing of J_N(theta) in function of model order n_theta, delimiting underfitting from oveerfitting

#note-box[We see that our loss function J_N(theta) is inversely proportional to the number of parameters in the model. The more parameters we have, the better we can fit the data. But this is not always a good thing since we may be fitting the noise of the data too. To avoid *overfitting*, we need to find a balance between the number of parameters and the goodness of fit.]


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


+ *FPE*: Final Prediction Error $"FPE" = (N+ n)/(N-n) J_N(theta)$
+ *AIC*: Akaike Information Criterion $"AC"(n) = ln(J_N(theta)) + 2n/(N)$
+ *MDL*: Minimum Description Length $"MDL"(n) = ln(J_N(theta)) + n ln(N)/(2N)$


