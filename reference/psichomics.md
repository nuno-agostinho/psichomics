# Start graphical interface of psichomics

Start graphical interface of psichomics

## Usage

``` r
psichomics(
  ...,
  launch.browser = TRUE,
  shinyproxy = FALSE,
  testData = FALSE,
  cache = getAnnotationHubOption("CACHE")
)
```

## Arguments

- ...:

  Arguments passed on to
  [`shiny::runApp`](https://rdrr.io/pkg/shiny/man/runApp.html)

  `port`

  :   The TCP port that the application should listen on. If the `port`
      is not specified, and the `shiny.port` option is set (with
      `options(shiny.port = XX)`), then that port will be used.
      Otherwise, use a random port between 3000:8000, excluding ports
      that are blocked by Google Chrome for being considered unsafe:
      3659, 4045, 5060, 5061, 6000, 6566, 6665:6669 and 6697. Up to
      twenty random ports will be tried.

  `host`

  :   The IPv4 address that the application should listen on. Defaults
      to the `shiny.host` option, if set, or `"127.0.0.1"` if not. See
      Details.

  `workerId`

  :   Can generally be ignored. Exists to help some editions of Shiny
      Server Pro route requests to the correct process.

  `quiet`

  :   Should Shiny status messages be shown? Defaults to FALSE.

  `display.mode`

  :   The mode in which to display the application. If set to the value
      `"showcase"`, shows application code and metadata from a
      `DESCRIPTION` file in the application directory alongside the
      application. If set to `"normal"`, displays the application
      normally. Defaults to `"auto"`, which displays the application in
      the mode given in its `DESCRIPTION` file, if any.

  `test.mode`

  :   Should the application be launched in test mode? This is only used
      for recording or running automated tests. Defaults to the
      `shiny.testmode` option, or FALSE if the option is not set.

- launch.browser:

  If true, the system's default web browser will be launched
  automatically after the app is started. Defaults to true in
  interactive sessions only. The value of this parameter can also be a
  function to call with the application's URL.

- shinyproxy:

  Boolean: prepare visual interface to run in Shinyproxy?

- testData:

  Boolean: load with test data

- cache:

  Character: path to `AnnotationHub` cache (used to load alternative
  splicing event annotation)

## Value

`NULL` (function is only used to modify the Shiny session's state or
internal variables)

## Examples

``` r
if (FALSE) { # \dontrun{
psichomics()
} # }
```
