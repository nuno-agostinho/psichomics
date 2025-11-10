# Create input to select a gene

Create input to select a gene

## Usage

``` r
selectizeGeneInput(
  id,
  label = "Gene",
  choices = NULL,
  multiple = FALSE,
  ...,
  placeholder = "Type to search for a gene..."
)
```

## Arguments

- id:

  Character: identifier

- label:

  Display label for the control, or `NULL` for no label.

- choices:

  List of values to select from. If elements of the list are named, then
  that name — rather than the value — is displayed to the user. It's
  also possible to group related inputs by providing a named list whose
  elements are (either named or unnamed) lists, vectors, or factors. In
  this case, the outermost names will be used as the group labels
  (leveraging the `<optgroup>` HTML tag) for the elements in the
  respective sublist. See the example section for a small demo of this
  feature.

- multiple:

  Is selection of multiple items allowed?

- ...:

  Arguments passed to the `options` list of
  [`selectizeInput()`](https://rdrr.io/pkg/shiny/man/selectInput.html)

- placeholder:

  Character: placeholder

## Value

HTML elements
