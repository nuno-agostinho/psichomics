# Change the value of a [`fileBrowserInput()`](https://nuno-agostinho.github.io/psichomics/reference/fileBrowserInput.md) on the client

Change the value of a
[`fileBrowserInput()`](https://nuno-agostinho.github.io/psichomics/reference/fileBrowserInput.md)
on the client

## Usage

``` r
updateFileBrowserInput(session, id, ..., value = NULL, ask = FALSE)
```

## Source

<https://github.com/wleepang/shiny-directory-input>

## Arguments

- session:

  Shiny session

- id:

  Character: identifier

- ...:

  Additional arguments passed to
  [`fileBrowser()`](https://nuno-agostinho.github.io/psichomics/reference/fileBrowser.md).
  Only used if `value = NULL`.

- value:

  Character: file or directory path

- ask:

  Boolean: ask user to pick a file using file browser?

## Value

`NULL` (function is only used to modify the Shiny session's state or
internal variables)

## Details

Sends a message to the client, telling it to change the value of the
input object. For
[`fileBrowserInput()`](https://nuno-agostinho.github.io/psichomics/reference/fileBrowserInput.md)
objects, this changes the value displayed in the text-field and triggers
a client-side change event. A directory selection dialogue is not
displayed.
