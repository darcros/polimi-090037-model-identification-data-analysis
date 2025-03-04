#import "@preview/suiji:0.3.0": uniform

#let random-pts-in-circle(rng, n, radius: 1.0) = {
  let pts = ()

  let i = 0
  while (i < n) {
    let (r, rho) = uniform(rng, low: 0.0, high: radius)
    let (r, theta) = uniform(rng, low: 0.0, high: 2 * calc.pi)
    pts.push((rho, theta))

    rng = r
    i += 1
  }

  pts = pts.map(p => {
    let (rho, theta) = p
    return (
      rho * calc.cos(theta),
      rho * calc.sin(theta)
    )
  })

  return (rng, pts)
}

#let random-series(rng, length, from: 0, to: 1, min: -1.0, max: +1.0) = {
  let values = ()

  let i = 0
  while (i < length) {
    let t = from + i * ((to - from) / length)
    
    let (r, v) = uniform(rng, low: min, high: max)
    values.push((t, v))

    rng = r
    i += 1
  }

  return (rng, values)
}

#let avg-series(series-arr) = {
  return array.zip(..series-arr).map(samples => {
    let values = samples.map(s => s.at(1))
    let avg = values.sum() / values.len()
  
    let (time, _v) = samples.at(0)
    return (time, avg)
  })
}
