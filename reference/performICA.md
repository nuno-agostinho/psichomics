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
#>             [,1]       [,2]        [,3]       [,4]
#> [1,] -0.89631995 -0.1491553  0.18945838  0.3721139
#> [2,]  0.06307973  0.1622659  0.94818442 -0.2657763
#> [3,] -0.03126581 -0.8929384  0.02927882 -0.4481363
#> [4,] -0.43778299  0.3925346 -0.25336648 -0.7681589
#> 
#> $A
#>           [,1]      [,2]       [,3]      [,4]
#> [1,] 4.2248552  73.93836   2.597938  5.523499
#> [2,] 0.1163438  12.08863  -2.863837  6.128464
#> [3,] 0.5475492 -14.85697 -13.685532 -4.086213
#> [4,] 0.6544529 -31.18521   1.756066  1.100146
#> 
#> $S
#>                         IC1         IC2         IC3         IC4
#> Alabama         1.168794631 -0.48270770  0.88298662  0.07134714
#> Alaska          0.453125439  3.82579694  0.48598371 -0.63197393
#> Arizona         0.518312176  0.58217556 -1.31978847 -1.86854897
#> Arkansas        0.181431175  0.33807159  1.02913692 -0.54603796
#> California      0.592621774  1.32029054 -2.08381860 -0.46505340
#> Colorado        0.069833513  1.88641910 -1.23796256  0.42071141
#> Connecticut    -0.882018658 -1.31737189 -0.76616251 -0.28851450
#> Delaware       -0.003410582 -0.86815799 -0.57664453 -2.22604984
#> Florida         1.966340432 -0.39098749 -0.65827572 -0.44247074
#> Georgia         1.796604003 -0.51410971  1.13923255  2.22725865
#> Hawaii         -0.808034586 -0.53940360 -1.01366453  2.35863678
#> Idaho          -1.121343658  0.33939501  0.41835908 -1.09869504
#> Illinois        0.833188854 -0.88306854 -0.98500212 -0.40648541
#> Indiana        -0.323089221  0.11977853  0.09271818  1.08839298
#> Iowa           -1.392940644 -0.12148627  0.40277301  0.13838833
#> Kansas         -0.499279367 -0.18763230 -0.01971790  0.54092421
#> Kentucky        0.074595888 -0.29727126  1.25126962  1.44593623
#> Louisiana       1.645661171 -1.05772787  0.59016567  0.70170856
#> Maine          -1.313466666 -0.34238022  0.78145030 -0.80501172
#> Maryland        1.066037369  0.27218569 -0.14624897 -1.44157381
#> Massachusetts  -0.538016481 -1.14078328 -1.33465145 -0.38421075
#> Michigan        1.017886805  0.90448087 -0.57125200  0.33483464
#> Minnesota      -1.223711785 -0.12034115 -0.20039136  0.31435444
#> Mississippi     1.728873580 -0.84245720  2.07269978 -0.04451201
#> Missouri        0.210489326  0.60957838 -0.32882573  0.65984993
#> Montana        -0.580258314  0.21349671  0.80072708  0.30595149
#> Nebraska       -0.860879714  0.04801370  0.10234736  0.13364723
#> Nevada          0.973280042  2.13137131 -1.24684925  1.12272753
#> New Hampshire  -1.397699772 -0.31969895  0.49533730 -0.02589556
#> New Jersey      0.035802681 -1.42091063 -1.33958855  0.54938058
#> New Mexico      1.008588651  0.74547902 -0.39256452 -0.79596941
#> New York        0.976007480 -0.83067055 -1.15235219 -0.12816406
#> North Carolina  1.483349478 -0.63421523  1.55576728 -2.80082644
#> North Dakota   -1.722528143  0.14181133  1.14757466 -0.54308040
#> Ohio           -0.227168556 -0.32544848 -0.51143719  1.20658834
#> Oklahoma       -0.262391047 -0.10848791 -0.19863497  0.06409705
#> Oregon         -0.617371431  1.53077171 -0.57250334 -0.22051572
#> Pennsylvania   -0.421101715 -0.96430790 -0.24482479  0.82105017
#> Rhode Island   -0.545411878 -2.26601198 -1.40359932 -1.60674305
#> South Carolina  1.478644096 -0.01797509  1.47925588 -0.67679635
#> South Dakota   -1.073335848  0.41917616  1.18029914 -0.22668081
#> Tennessee       0.946888068  0.32792634  0.77506855  1.45005803
#> Texas           1.033801108 -0.82561285 -0.50752967  1.40314525
#> Utah           -0.942819398  0.25347766 -1.27597419  0.09846937
#> Vermont        -1.559859254  1.07837480  1.89721083 -0.24768643
#> Virginia        0.061801863 -0.04861596  0.26822174  0.47320082
#> Washington     -0.771900834  0.93425937 -0.91365404 -0.20666954
#> West Virginia  -0.757694079 -0.07356772  1.83388551  0.17963754
#> Wisconsin      -1.279250227 -0.68576641 -0.06757531  0.50948274
#> Wyoming        -0.196977744 -0.39515422  0.38702299 -0.49161362
#> 
```
