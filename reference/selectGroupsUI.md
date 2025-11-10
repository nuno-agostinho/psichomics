# Group selection

Group selection interface and logic

## Usage

``` r
selectGroupsUI(
  id,
  label,
  type,
  placeholder = "Type to search groups",
  noGroupsLabel = NULL,
  groupsLabel = NULL,
  maxItems = NULL,
  returnAllDataLabel = NULL,
  returnAllDataValue = FALSE
)

selectGroupsServer(session, id, type, preference = NULL)

getSelectedGroups(input, id, type, filter = NULL)
```

## Arguments

- id:

  Character: identifier

- label:

  Character: `selectize` label

- type:

  Character: type of groups (either `Patients`, `Samples`, `ASevents` or
  `Genes`)

- placeholder:

  Character: `selectize` placeholder

- noGroupsLabel:

  Character: label to explicitly allow to select no groups (if `NULL`,
  this option is not displayed to the user)

- groupsLabel:

  Character: label to explicitly allow to select groups (only required
  if `noGroupsLabel` is not `NULL`)

- maxItems:

  Numeric: maximum number of groups to select

- returnAllDataLabel:

  Character: label to allow to return data outside selected groups as
  belonging to an outside group (if `NULL`, this option is not displayed
  to the user)

- returnAllDataValue:

  Boolean: default value to whether return all data or not (only
  required if `returnAllDataLabel` is not `NULL`)

- session:

  Shiny session

- preference:

  Character: name of groups to pre-select, when available (if `NULL`,
  all groups will be pre-selected)

- input:

  Shiny input

- filter:

  Character: get groups only if they are present in this argument (if
  TCGA-styled gene symbols, they will be "converted" to gene symbols
  alone)

## Value

`selectGroupsUI`: Interface for group selection

`selectGroupsServer`: Server logic for group selection

`getSelectedGroups`: List with selected groups (or `NULL` when no groups
are selected)

## Note

To allow the user to (explicitly) select no groups, pass the
`noGroupsLabel` and `groupsLabel` arguments.
