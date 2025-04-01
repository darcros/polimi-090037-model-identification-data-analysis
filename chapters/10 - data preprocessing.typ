#import "/prelude.typ": *


= Data preprocessing
Let's now tackle non-stationary processes. An approach to this problem is preprocess the data to make extract the stationary part of the process. This is called *data preprocessing*.

Two possible phenomena can be observed in the data that make it non-stationary are:
- *trend*: a long term increase or decrease in the data.
- *seasonality*: a periodic pattern in the data.

== Trend removal


#figure(
  cetz.canvas({
    import suiji: gen-rng
    import cetz.draw: *
    import cetz-plot: *
    import "../util.typ": random-series
    let rng = gen-rng(2)
    let (rng, realization) = random-series(rng, 50, from: 0, to: 10, min: -10, max: 10)
    let linear_trend = range(-2, 10).map(x => (x, 2 * x + 10))
    let trended_realization = range(0, 50).map(x => (x / 5, 2 * x / 5 + 10 + realization.at(x).at(1)))


    plot.plot(
      size: (6, 6),
      axis-style: "school-book",
      x-min: -2,
      x-max: +10,
      y-max: +30,
      x-tick-step: none,
      x-label: [$t$],
      y-tick-step: none,
      y-label: [$v(t, overline(s))$],
      {
        plot.add(label: [SSP realization $tilde(y)(t)$], realization, line: "spline")
        plot.add(label: [linear trend $k t+m$], linear_trend, line: "spline")
        plot.add(label: [trended realization $y(t)$], trended_realization, line: "spline")
      },
    )
  }),
)


We have to estimate the $(hat(m), hat(k))$ parameters of the linear trend $k t + m$ and remove it from the data. This is done by estimating the linear regression of the data with OLS. The linear regression is given by the following equation:
$ (hat(m), hat(k)) = argmin_{m,k} sum_{i=1}^N (y_i - (k t_i + m))^2 $.

The OLS solution is given by the following vectorial equations:
$
  (hat(m), hat(k))^T = mat(sum^N(t^2), sum^N(t); sum^N(t), N)^{-1} dot mat(sum_((t=1))^N t y(t); sum_((t=1))^N y(t); )
$

Then, the detrended data is given by:
$ tilde(y)(t) = y(t) - (hat(m) + hat(k) t) $.

The detrended data is a stationary stochastic process. We can obtain the estimator of the trended realization by using the following equation:
$ y(t +1 | t) = hat(m) + hat(k) t + hat(tilde(y))(t +1 | t) $.


== Seasonality removal
#figure(
  cetz.canvas({
    import suiji: gen-rng
    import cetz.draw: *
    import cetz-plot: *
    import "../util.typ": random-series
    let rng = gen-rng(2)
    let (rng, realization) = random-series(rng, 50, from: 0, to: 10, min: -5, max: 5)
    let trend(x) = (
      return 20 * calc.sin(x) + 10 * calc.sin(2 * x) + 5 * calc.sin(4 * x)
    )

    let seasonal_trend = range(-2, 10).map(x => (x, trend(x)))
    let trended_realization = range(0, 50).map(x => (x / 5, trend(x / 5) + realization.at(x).at(1)))

    plot.plot(
      size: (6, 6),
      axis-style: "school-book",
      x-min: -2,
      x-max: +10,
      y-max: +30,
      x-tick-step: none,
      x-label: [$t$],
      y-tick-step: none,
      y-label: [$v(t, overline(s))$],
      {
        plot.add(label: [SSP realization $tilde(y)(t)$], realization, line: "spline")
        plot.add(label: [seasonal trend $k sin(t)$], seasonal_trend, line: "spline")
        plot.add(label: [trended realization $y(t)$], trended_realization, line: "spline")
      },
    )
  }),
)

$ y(t) = tilde(y)(t) + s(t) "with" s(t) = s(t+(2pi) / T k ), k in ZZ $.

The way to remove seasonality, and therefore periodic signals from the data, is to use a the spectral density of the data to find the main frequencies of the dataset by identifying peaks in the spectrum. If there is more than one peak, the process could be multiseasonal.

#figure(
  cetz.canvas({
    import cetz-plot: *
    import cetz.draw: *
    import suiji: gen-rng
    import "../util.typ": random-series
    let rng = gen-rng(2)
    let (rng, realization) = random-series(rng, 50, from: 0, to: 10, min: -5, max: 5)

    let trend(x) = (
      return 20 * calc.sin(x) + 10 * calc.sin(2 * x) + 5 * calc.sin(4 * x)
    )

    let seasonal_trend = range(-2, 10).map(x => (x, trend(x)))
    let trended_realization = range(0, 50).map(x => (x / 5, trend(x / 5) + realization.at(x).at(1)))

    let dft(omega) = (
      (1 / trended_realization.len())
        * trended_realization
          .map(value => {
            let (t, y) = value


            let re = y * calc.cos(omega * t)
            let im = y * calc.sin(omega * t)

            let mod = (re * re + im * im)
            return mod
          })
          .sum()
    )


    let dft_results = range(0, 200).map(omega => (omega / 100, dft(omega / 100)))

    plot.plot(
      size: (6, 6),
      axis-style: "school-book",
      x-min: 0,
      x-tick-step: none,
      x-label: [Frequency],
      y-tick-step: none,
      y-label: [Magnitude],
      {
        plot.add(label: [DFT Magnitude], dft_results, line: "spline")
      },
    )
  }),
)
