# Perform independent component analysis after processing missing values

Perform independent component analysis after processing missing values

## Usage

``` r
performICA(
  data,
  n.comp = min(5, ncol(data)),
  center = TRUE,
  scale. = FALSE,
  missingValues = round(0.05 * nrow(data)),
  alg.typ = c("parallel", "defaltion"),
  fun = c("logcosh", "exp"),
  alpha = 1,
  ...
)
```

## Arguments

- data:

  an optional data frame (or similar: see
  [`model.frame`](https://rdrr.io/r/stats/model.frame.html)) containing
  the variables in the formula `formula`. By default the variables are
  taken from `environment(formula)`.

- n.comp:

  number of components to be extracted

- center:

  a logical value indicating whether the variables should be shifted to
  be zero centered. Alternately, a vector of length equal the number of
  columns of `x` can be supplied. The value is passed to `scale`.

- scale.:

  a logical value indicating whether the variables should be scaled to
  have unit variance before the analysis takes place. The default is
  `FALSE` for consistency with S, but in general scaling is advisable.
  Alternatively, a vector of length equal the number of columns of `x`
  can be supplied. The value is passed to
  [`scale`](https://rdrr.io/r/base/scale.html).

- missingValues:

  Integer: number of tolerated missing values per column to be replaced
  with the mean of the values of that same column

- alg.typ:

  if `alg.typ == "parallel"` the components are extracted simultaneously
  (the default). if `alg.typ == "deflation"` the components are
  extracted one at a time.

- fun:

  the functional form of the \\G\\ function used in the approximation to
  neg-entropy (see ‘details’).

- alpha:

  constant in range \[1, 2\] used in approximation to neg-entropy when
  `fun == "logcosh"`

- ...:

  Arguments passed on to
  [`fastICA::fastICA`](https://rdrr.io/pkg/fastICA/man/fastICA.html)

## Value

ICA result in a `prcomp` object

## See also

Other functions to analyse independent components:
[`plotICA()`](https://nuno-agostinho.github.io/psichomics/reference/plotICA.md)

## Examples

``` r
performICA(USArrests)
#> $X
#>                Murder Assault UrbanPop    Rape
#> Alabama         5.412   65.24    -7.54  -0.032
#> Alaska          2.212   92.24   -17.54  23.268
#> Arizona         0.312  123.24    14.46   9.768
#> Arkansas        1.012   19.24   -15.54  -1.732
#> California      1.212  105.24    25.46  19.368
#> Colorado        0.112   33.24    12.46  17.468
#> Connecticut    -4.488  -60.76    11.46 -10.132
#> Delaware       -1.888   67.24     6.46  -5.432
#> Florida         7.612  164.24    14.46  10.668
#> Georgia         9.612   40.24    -5.54   4.568
#> Hawaii         -2.488 -124.76    17.46  -1.032
#> Idaho          -5.188  -50.76   -11.54  -7.032
#> Illinois        2.612   78.24    17.46   2.768
#> Indiana        -0.588  -57.76    -0.54  -0.232
#> Iowa           -5.588 -114.76    -8.54  -9.932
#> Kansas         -1.788  -55.76     0.46  -3.232
#> Kentucky        1.912  -61.76   -13.54  -4.932
#> Louisiana       7.612   78.24     0.46   0.968
#> Maine          -5.688  -87.76   -14.54 -13.432
#> Maryland        3.512  129.24     1.46   6.568
#> Massachusetts  -3.388  -21.76    19.46  -4.932
#> Michigan        4.312   84.24     8.46  13.868
#> Minnesota      -5.088  -98.76     0.46  -6.332
#> Mississippi     8.312   88.24   -21.54  -4.132
#> Missouri        1.212    7.24     4.46   6.968
#> Montana        -1.788  -61.76   -12.54  -4.832
#> Nebraska       -3.488  -68.76    -3.54  -4.732
#> Nevada          4.412   81.24    15.46  24.768
#> New Hampshire  -5.688 -113.76    -9.54 -11.732
#> New Jersey     -0.388  -11.76    23.46  -2.432
#> New Mexico      3.612  114.24     4.46  10.868
#> New York        3.312   83.24    20.46   4.868
#> North Carolina  5.212  166.24   -20.54  -5.132
#> North Dakota   -6.988 -125.76   -21.54 -13.932
#> Ohio           -0.488  -50.76     9.46   0.168
#> Oklahoma       -1.188  -19.76     2.46  -1.232
#> Oregon         -2.888  -11.76     1.46   8.068
#> Pennsylvania   -1.488  -64.76     6.46  -6.332
#> Rhode Island   -4.388    3.24    21.46 -12.932
#> South Carolina  6.612  108.24   -17.54   1.268
#> South Dakota   -3.988  -84.76   -20.54  -8.432
#> Tennessee       5.412   17.24    -6.54   5.668
#> Texas           4.912   30.24    14.46   4.268
#> Utah           -4.588  -50.76    14.46   1.668
#> Vermont        -5.588 -122.76   -33.54 -10.032
#> Virginia        0.712  -14.76    -2.54  -0.532
#> Washington     -3.788  -25.76     7.46   4.968
#> West Virginia  -2.088  -89.76   -26.54 -11.932
#> Wisconsin      -5.188 -117.76     0.46 -10.432
#> Wyoming        -0.988   -9.76    -5.54  -5.632
#> attr(,"scaled:center")
#>        Murder       Assault      UrbanPop          Rape 
#> -2.842171e-16  9.094947e-15 -6.252776e-15  8.526513e-16 
#> 
#> $K
#>               [,1]         [,2]        [,3]        [,4]
#> [1,] -0.0005031233  0.003185718 -0.01243588 -0.40479570
#> [2,] -0.0120064074  0.004176394  0.01051799  0.01584251
#> [3,] -0.0005589971 -0.069430563  0.03121730 -0.02366681
#> [4,] -0.0009066803 -0.014266122 -0.15162666  0.02942629
#> 
#> $W
#>            [,1]        [,2]        [,3]        [,4]
#> [1,] -0.3272940 -0.94429966 -0.02925497 -0.01791523
#> [2,]  0.2522182 -0.07988422 -0.66871761  0.69485341
#> [3,]  0.4959778 -0.20255298  0.70075675  0.47108204
#> [4,]  0.7637216 -0.24675688 -0.24678106 -0.54308332
#> 
#> $A
#>             [,1]      [,2]      [,3]      [,4]
#> [1,] -0.83164439 27.496932 -1.679598 -1.641239
#> [2,]  3.92134419 77.722168  4.499120  7.331951
#> [3,] -0.07676045  2.141101 10.241290 -2.358271
#> [4,]  1.58637674  2.204880 -8.796661 -4.895097
#> 
#> $S
#>                        IC1         IC2         IC3          IC4
#> Alabama        -0.09449110  0.83631141 -0.00732968  1.294389562
#> Alaska          0.40717676  1.17895294 -3.42844353 -1.472300260
#> Arizona         1.81467992  0.97323601  0.21042435 -1.247542710
#> Arkansas        0.51751879  0.06841125 -0.86269229  0.698388250
#> California      0.37335669  1.28203521  0.09011010 -2.204951917
#> Colorado       -0.52095560  0.69068271 -0.84785891 -1.950820618
#> Connecticut     0.39534489 -0.95597259  1.44144313 -0.189028385
#> Delaware        2.26855195  0.03793316  0.94412073 -0.048947617
#> Florida         0.38374628  1.94884283  0.80884355  0.221349760
#> Georgia        -2.26982407  1.27146999  0.01142288  1.726777023
#> Hawaii         -2.28989784 -0.80332272  1.05337559 -0.732120046
#> Idaho           1.12285326 -1.02787566 -0.63867754 -0.171810595
#> Illinois        0.42067854  0.82452406  1.31953235 -0.107225346
#> Indiana        -1.07939230 -0.35600590 -0.12717930 -0.062665585
#> Iowa           -0.07580912 -1.44366582 -0.21560434 -0.004088263
#> Kansas         -0.50939954 -0.53997831  0.15398870 -0.051928527
#> Kentucky       -1.43215202 -0.31158418 -0.38413986  1.206082283
#> Louisiana      -0.71157202  1.19679959  0.69121407  1.500412771
#> Maine           0.87333007 -1.44250298 -0.27814040  0.424552704
#> Maryland        1.38081687  1.17730955 -0.12991985  0.021268558
#> Massachusetts   0.46849075 -0.47242248  1.63193282 -0.643342293
#> Michigan       -0.42258871  1.26091697 -0.34422904 -0.636895460
#> Minnesota      -0.25719862 -1.17099477  0.13848102 -0.440875337
#> Mississippi     0.01538190  1.06549723 -0.33767297  2.597548784
#> Missouri       -0.69882097  0.36257553 -0.27822425 -0.512054042
#> Montana        -0.29457702 -0.68316335 -0.64991598  0.375728830
#> Nebraska       -0.10116024 -0.83951298 -0.15139963 -0.183897991
#> Nevada         -1.27110597  1.57279745 -0.95228726 -1.819042227
#> New Hampshire   0.09870188 -1.50035992 -0.11054592  0.169585797
#> New Jersey     -0.47133303 -0.02892030  1.94609988 -0.326020598
#> New Mexico      0.71420128  1.23930305 -0.36079679 -0.429574662
#> New York        0.13470029  0.99139572  1.39116891 -0.224911490
#> North Carolina  2.76529846  1.11471858 -0.35615731  1.962466327
#> North Dakota    0.60208537 -1.81642144 -0.89982085  0.357080335
#> Ohio           -1.17635193 -0.24490367  0.58890363 -0.290441964
#> Oklahoma       -0.04715184 -0.23761108  0.18455113 -0.177317602
#> Oregon          0.16446492 -0.14129310 -0.99959108 -1.433388290
#> Pennsylvania   -0.75022078 -0.59920918  0.93874765  0.195317552
#> Rhode Island    1.74872549 -0.64333488  2.54826451 -0.135742538
#> South Carolina  0.61350611  1.15009861 -0.73337283  1.611213973
#> South Dakota    0.24519607 -1.15945721 -1.08838647  0.428020446
#> Tennessee      -1.50413977  0.74866363 -0.57460381  0.744599918
#> Texas          -1.39696370  0.84367322  1.09704592  0.331636005
#> Utah           -0.07116232 -0.59954293  0.43862118 -1.426205230
#> Vermont         0.24941394 -1.62570154 -2.06741533  0.526774770
#> Virginia       -0.47254364 -0.02794304 -0.08470155  0.266068117
#> Washington      0.18937013 -0.34849636 -0.32680775 -1.442925300
#> West Virginia  -0.14793289 -1.11190127 -1.00540999  1.306085915
#> Wisconsin      -0.41971477 -1.37991255  0.53639669 -0.053434541
#> Wyoming         0.51886918 -0.32413848  0.07663573  0.454151753
#> 
```
