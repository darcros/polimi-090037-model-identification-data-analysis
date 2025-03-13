// random numbers
#import "@preview/suiji:0.3.0"

// drawings and plots
#import "@preview/cetz:0.3.3"
#import "@preview/cetz-plot:0.1.1"

// diagrams
#import "@preview/fletcher:0.5.6"

// theorem boxes
#import "@preview/theorion:0.3.2": *
#import cosmos.rainbow: *

#let (properties-counter, properties-box, properties, show-properties) = make-frame(
  "properties",
  (en: "Properties", it: "Proprietà"),
  inherited-levels: 2, // just to make it the same as the other boxes
  render: render-fn.with(fill: green.darken(40%)), // render-fn from cosoms.rainbow
)

// custom math shortcuts

#let argmin = $op("argmin", limits: #true)$
#let transposed = $sans(upright(T))$

#let WN = $op("WN")$
#let AR = $op("AR")$
#let MA = $op("MA")$
#let ARMA = $op("ARMA")$
#let ARMAX = $op("ARMAX")$

#let longdivision(divisor: none, quotient: none, steps: none) = {
  grid(
    columns: 2,
    inset: 5pt,
    grid.cell(
      rowspan: 2,
      steps,
    ),
    grid.vline(),
    divisor,
    grid.hline(),
    quotient
  )
}
