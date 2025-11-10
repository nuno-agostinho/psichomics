# Modified colour input with 100% width

Modified colour input with 100% width

## Usage

``` r
colourInputMod(...)
```

## Arguments

- ...:

  Arguments passed on to
  [`colourpicker::colourInput`](https://rdrr.io/pkg/colourpicker/man/colourInput.html)

  `inputId`

  :   The `input` slot that will be used to access the value.

  `label`

  :   Display label for the control, or \``NULL` for no label.

  `value`

  :   Initial value (can be a colour name or HEX code)

  `showColour`

  :   Whether to show the chosen colour as text inside the input, as the
      background colour of the input, or both (default).

  `palette`

  :   The type of colour palette to allow the user to select colours
      from. `square` (default) shows a square colour palette that allows
      the user to choose any colour, while `limited` only gives the user
      a predefined list of colours to choose from.

  `allowedCols`

  :   A list of colours that the user can choose from. Only applicable
      when `palette == "limited"`. The `limited` palette uses a default
      list of 40 colours if `allowedCols` is not defined. If the colour
      specified in `value` is not in the list, the default colour will
      revert to black.

  `allowTransparent`

  :   If `TRUE`, enables a slider to choose an alpha (transparency)
      value for the colour. When a colour with opacity is chosen, the
      return value is an 8-digit HEX code.

  `returnName`

  :   If `TRUE`, then return the name of an R colour instead of a HEX
      value when possible.

  `closeOnClick`

  :   If `TRUE`, then the colour selection panel will close immediately
      after selecting a colour.

  `width`

  :   The width of the input, e.g. `"400px"` or `"100%"`

## Value

HTML elements
