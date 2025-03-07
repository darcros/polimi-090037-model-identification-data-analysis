#import "/prelude.typ": *

= Modeling systems

#{
  import fletcher: diagram, node, edge
  figure(
    diagram(
      node-shape: "rect",
      node-stroke: 1pt,
      edge((0, 0), "rr", "-|>")[$u(t)$],
      node((2,0))[system],
      edge("r", "-|>")[$y(t)$],
      
      edge((0, 0), "r,d,r", "-|>"),
      node((2,1))[model],
      edge("r", "-|>")[$hat(y)(t)$],
    )
  )
}

We create models to predict the output $y(t)$ of a system.
The model takes the same input $u(t)$ as the system and produces an approximation of the output $hat(y)(t)$.

Models can be:
/ white box: derived from first principles and known physical laws 
/ black box: derived from historical data + statistical analysis
/ gray box: derived from a mix of historical data (+ analysis) and simplified physical laws

Some challenges in model identification are: 
- *Uncertainty*: #underline[structural] if the physical laws that govern the system are unknown or unavailable and #underline[parametric] if we are oblivious to its parameters.  
- *Complexity and model purpose*: the model must be complex enough to describe the model for the intended use but simple enough that computations and predictions are feasible  

== Modeling error

#{
  import fletcher: diagram, node, edge
  figure(
    diagram(
      node-shape: "rect",
      node-stroke: 1pt,
      edge((0, 0), "rr", "-|>")[$u(t)$],
      node((2,0))[system],
      edge("r,d", "-|>", label-pos: 75%)[$+$],
      
      edge((0, 0), "r,d,r", "-|>"),
      node((2,1))[model],
      edge("r", "-|>")[$-$],
      node((3,1), shape: circle, radius: 4pt, stroke: 1pt),
      edge("r", "-|>")[$epsilon(t)$],
    )
  )
}

== Static and dynamic systems

#definition(title: "Static system")[
  The output $y(t)$ of the system is determined _only_ by its input $u(t)$.

  The system does not change over time.
]

#definition(title: "Dynamic system")[
  The output $y(t)$ of the system is determined by
  - its input $u(t)$
  - the previous state of the system.

  The system changes over time. The change can be 
  - *Discrete*: described by difference equations: $y(t) = a dot y(t-1) + e(t)$
  - *Continuous*: described by differential equations $dot(y)(t) = a dot y(t) + e(t) $
]
