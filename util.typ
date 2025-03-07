#import "@preview/suiji:0.3.0": uniform

#let rng-iter(rng, iterations, function) = {
  let values = ()
  let r = rng

  let i = 0
  while (i < iterations) {
    let (ret-r, value) = function(r, i)
    values.push(value)
    r = ret-r
    i += 1
  }

  return (r, values)
}

#let random-pts-in-circle(rng, n, radius: 1.0) = rng-iter(
  rng,
  n,
  (r, i) => {
    let (r, rho) = uniform(r, low: 0.0, high: radius)
    let (r, theta) = uniform(r, low: 0.0, high: 2 * calc.pi)

    let x = rho * calc.cos(theta)
    let y = rho * calc.sin(theta)
    let point = (x, y)

    return (r, point)
  },
)

#let random-series(rng, length, from: 0, to: 1, min: -1.0, max: +1.0) = rng-iter(
  rng,
  length,
  (r, i) => {
    let time = from + i * ((to - from) / length)
    let (r, value) = uniform(r, low: min, high: max)
    return (r, (time, value))
  },
)

// produces a new time series that is the average of the input series
// it assumes that all series have the same time points
#let avg-series(series-arr) = {
  return array
    .zip(..series-arr)
    .map(samples => {
      let values = samples.map(s => s.at(1))
      let avg = values.sum() / values.len()

      let (time, _v) = samples.at(0)
      return (time, avg)
    })
}
