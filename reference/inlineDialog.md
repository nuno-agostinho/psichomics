# Alert in the style of a dialogue box with a button

Alert in the style of a dialogue box with a button

## Usage

``` r
inlineDialog(
  description,
  ...,
  buttonLabel = NULL,
  buttonIcon = NULL,
  buttonId = NULL,
  id = NULL,
  type = c("error", "warning"),
  bigger = FALSE
)

errorDialog(description, ...)

warningDialog(description, ...)
```

## Arguments

- description:

  Character: description

- ...:

  Extra parameters when creating the alert

- buttonLabel:

  Character: button label

- buttonIcon:

  Character: button icon

- buttonId:

  Character: button identifier

- id:

  Character: identifier

- type:

  Character: type of alert (error or warning)

- bigger:

  Boolean: wrap the `description` in a `h4` tag?

## Value

HTML elements
