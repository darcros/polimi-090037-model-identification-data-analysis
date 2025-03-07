
#let is-upper(str) = {
  let A = "A".to-unicode()
  let Z = "Z".to-unicode()
  let char = str.to-unicode()

  return A <= char and char <= Z
}

#let multiline-title(title) = {
  let words = title.split()

  set text(size: 4em, weight: "semibold")
  set par(leading: 0.5em)

  for (i, word) in words.enumerate() {
    if is-upper(word.first()) {
      if (i != 0) {
        linebreak()
      }

      word
    } else {
      sym.space
      text(style: "italic", fill: luma(50%), weight: "regular", word)
    }
  }

  parbreak()
}

#let line-and-subtitle(subtitle) = context {
  let gap = 1.25em

  let t = text(size: 2em, top-edge: "x-height", bottom-edge: "baseline", subtitle)
  let (width, height) = measure(t)

  box(line(start: (0pt, -height / 2), length: 100% - width - gap, stroke: luma(50%)))
  h(1fr)
  t
}

#let template(
  title: none,
  subtitle: none,
  authors: (),
  academic-year: none,
  doc,
) = {
  if (title == none) {
    panic("Title is required")
  }

  // metdata
  set document(
    title: (title, subtitle).join(" - "),
    author: authors,
  )

  // title page
  multiline-title(title)

  if (subtitle != none) {
    line-and-subtitle(subtitle)
  }

  v(1fr)

  {
    set text(size: 1.25em)
    set list(marker: none, body-indent: 0pt)

    let h = text.with(weight: "bold")

    grid(
      columns: 2,
      row-gutter: 1.2em,
      column-gutter: 1.2em,
      ..if academic-year != none {
        (h[Academic Year], academic-year)
      },
      ..if authors.len() != 0 {
        let txt = if authors.len() == 1 [Author] else [Authors]
        (h(txt), list(..authors))
      }
    )
  }

  v(3cm)

  pagebreak(weak: true)

  // outline page
  show outline.entry.where(level: 1): set block(above: 1.2em)
  show outline.entry.where(level: 1): set text(weight: "bold")
  outline()

  pagebreak(weak: true)

  // content
  set page(numbering: "1", number-align: center)
  set text(font: "Libertinus Serif", lang: "en")

  set heading(numbering: "1.1")
  // set par(leading: 0.75em)
  set par(leading: 0.75em, spacing: 1.2em)

  doc
}
