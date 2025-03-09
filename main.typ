#import "prelude.typ": *
#show: show-theorion
#show: show-properties

#import "template.typ": template

// TODO: find a glossary/abbreviations package
// abbreviations:
// - SP: stochastic process
// - SSP: stationary stochastic process
// TODO: find and list other abbreviations

#show: template.with(
  title: "Model Identification and Data Analysis",
  subtitle: "Course notes",
  academic-year: "2024-25",
  authors: (
    "Dario Crosa",
    "Francesco Genovese",
  ),
)

#show heading.where(level: 1): it => {
  pagebreak()
  it
}

#include "chapters/01 - modeling systems.typ"
#include "chapters/02 - stochastic processes.typ"
#include "chapters/03 - models.typ"
#include "chapters/04 - frequency domain.typ"
#include "chapters/05 - representations.typ"
#include "chapters/06 - prediction.typ"
#include "chapters/07 - model identification.typ"
