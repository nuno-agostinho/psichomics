# Create script for auto-completion of text input

Uses the JavaScript library `jquery.textcomplete`

## Usage

``` r
textSuggestions(id, words, novalue = "No matching value", char = " ")
```

## Arguments

- id:

  Character: input ID

- words:

  Character: words to suggest

- novalue:

  Character: string when there's no matching values

- char:

  Character to succeed accepted word

## Value

HTML string with the JavaScript script prepared to run

## Examples

``` r
words <- c("tumor_stage", "age", "gender")
psichomics:::textSuggestions("textareaid", words)
#> <script> textareaid_words = ["tumor_stage", "age", "gender"]; $("#textareaid").textcomplete([{
#>             match: /([a-zA-Z0-9_\.]{1,})$/,
#>             search: function(term, callback) {
#>                 var words = textareaid_words, sorted = [];
#>                 for (i = 0; i < words.length; i++) {
#>                     sorted[i] = fuzzy(words[i], term);
#>                 }
#>                 sorted.sort(fuzzy.matchComparator);
#>                 sorted = sorted.map(function(i) { return i.term; });
#>                 callback(sorted);
#>             },
#>             index: 1,
#>             cache: true,
#>             replace: function(word) {
#>             return word + " ";
#>         }}], { noResultsMessage: "No matching value"}); </script>
```
