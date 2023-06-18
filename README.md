
<!-- README.md is generated from README.Rmd. Please edit that file -->

xtune.mJAM: External-informed tuning of regularized regression for
multi-population summary statistics

<!-- badges: start -->

<!-- badges: end -->

## 📗 Introduction

Regularized regression is a common approach to fit high-dimensional data
such as genetics data. In standard regularized regression, a single
penalty parameter \(\lambda\) applied equally to all regression
coefficients to control the amount of regularization in the model.
`xtune`, a previous study proposed a novel method that model the penalty
parameter as a log-linear function of the prior data \(**Z**\) to
introduce the feature-specific shrinkage parameters. To this end, we
extend the `xtune` from modeling with individual level data to summary
statistics by incorporating the `mJAM` method which utilizes the
Cholesky decomposition to perform the joint analysis of marginal summary
statistics for multi-population GWAS studies.

## 📙 Installation

`xtune.mJAM` can be installed from Github using the following command:

``` r
# install.packages("devtools")

library(devtools)
devtools::install_github("JingxuanH/xtune.mJAM")

library("xtune.mJAM")
```

## ✍ Citation

  - **`xtune`** package:

<!-- end list -->

``` r
citation("xtune.mJAM")
#> Warning in citation("xtune.mJAM"): no date field in DESCRIPTION file of package
#> 'xtune.mJAM'
#> Warning in citation("xtune.mJAM"): could not determine year for 'xtune.mJAM'
#> from package DESCRIPTION file
#> 
#> To cite package 'xtune.mJAM' in publications use:
#> 
#>   First Last (NA). xtune.mJAM: What the Package Does (One Line, Title
#>   Case). R package version 0.0.0.9000.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {xtune.mJAM: What the Package Does (One Line, Title Case)},
#>     author = {First Last},
#>     note = {R package version 0.0.0.9000},
#>   }
```

Feel free to contact `hejingxu@usc.edu` if you have any questions.

## 📘 Examples

``` r
library(xtune.mJAM)

## load the example data
data(example)
```

``` r
example$beta.gwas
#> [[1]]
#>   [1] -9.544658e-02 -2.437938e-02 -1.475468e-02 -6.582795e-03 -4.665708e-03
#>   [6] -2.810554e-03 -2.490440e-03 -2.471588e-04  2.212117e-04  3.073036e-03
#>  [11]  1.200192e-03  2.916514e-04 -7.499704e-05 -2.383013e-03 -4.597121e-04
#>  [16]  1.375390e-03  2.616503e-03 -2.301694e-04 -1.654624e-03 -4.666318e-04
#>  [21] -2.092457e-03 -1.762126e-03 -1.880288e-03 -1.334133e-03 -1.022047e-03
#>  [26] -2.257442e-03 -1.523862e-03 -1.270806e-03 -3.426662e-03  2.404001e-04
#>  [31] -9.339270e-04  1.111311e-03 -1.760466e-04 -1.629154e-03 -2.566436e-03
#>  [36] -5.207323e-03  1.033459e-03  3.339694e-05  4.167726e-04 -6.884198e-04
#>  [41] -2.143637e-03 -4.089908e-03 -5.443711e-04  1.377338e-03 -1.461230e-03
#>  [46] -3.903402e-03 -1.325282e-03 -2.742050e-03 -1.128299e-03 -3.877707e-03
#>  [51]  7.704498e-04 -1.495431e-03  8.669007e-04  1.183010e-03  3.548633e-03
#>  [56]  4.585852e-03  2.946819e-03  6.949832e-04 -1.798042e-03 -1.091522e-03
#>  [61] -1.759850e-04 -2.305814e-04  7.694513e-04 -1.329277e-03 -5.979818e-03
#>  [66] -2.329831e-03 -1.765296e-03 -1.971044e-03 -1.028607e-03  5.253111e-05
#>  [71] -7.143652e-04  2.688086e-04 -6.004806e-04 -5.602000e-04  1.646267e-03
#>  [76]  2.189403e-04  2.432982e-04 -1.927704e-03 -1.815775e-03 -2.269354e-03
#>  [81]  1.145841e-03  3.551188e-03 -6.516718e-04 -1.631185e-03  3.744704e-04
#>  [86]  4.074742e-04 -1.143057e-03 -1.977140e-03 -8.062651e-04 -1.692973e-03
#>  [91] -1.175792e-03  1.529926e-03  3.498734e-03 -6.601487e-04  1.083887e-03
#>  [96]  1.760026e-04  1.935054e-03  1.735156e-03  3.774327e-04 -1.408765e-03
#> 
#> [[2]]
#>   [1] -9.572282e-02 -4.431360e-02 -2.770634e-02 -1.484879e-02 -1.047627e-02
#>   [6] -6.523546e-03 -5.147496e-03 -4.176247e-03 -3.266516e-03 -1.406959e-03
#>  [11] -5.156374e-03 -3.809519e-03 -1.225351e-03  7.394067e-04 -7.587754e-04
#>  [16]  7.603743e-04  3.259096e-04  1.191372e-03 -1.781579e-03  1.413486e-04
#>  [21] -1.309462e-03  7.241712e-04 -1.254015e-03  1.020813e-03 -1.412290e-03
#>  [26]  6.811516e-04 -2.403577e-03 -2.508411e-03 -2.321398e-03 -3.356246e-03
#>  [31] -3.339043e-04 -2.177174e-03  2.384995e-03  2.847065e-04  1.241399e-03
#>  [36]  4.046761e-04 -1.071847e-03 -4.641617e-03 -2.505712e-03 -2.346008e-03
#>  [41] -1.623907e-04  3.184793e-05  6.367788e-04  3.033111e-03  9.367475e-04
#>  [46]  8.327162e-06  4.472803e-04  3.612153e-04  2.213042e-03  9.737805e-04
#>  [51] -8.518014e-04 -1.841310e-03 -5.331720e-04 -2.413834e-04 -3.821455e-04
#>  [56] -9.557130e-04 -5.738267e-04  9.277912e-04 -5.373631e-04  1.477894e-03
#>  [61] -2.260612e-04  3.439496e-04 -2.245204e-04  1.970623e-03  3.684888e-04
#>  [66] -9.424814e-04 -3.160068e-03 -1.998354e-03  1.076007e-03  5.638236e-04
#>  [71]  2.428916e-03  6.339957e-04  1.021333e-03  1.491315e-03 -3.388989e-04
#>  [76] -1.335919e-04 -7.745159e-04 -1.245713e-03 -1.968375e-03 -8.262028e-05
#>  [81] -3.348897e-04  8.713417e-04  7.228913e-04 -1.127191e-03  1.447199e-03
#>  [86]  1.159202e-03  1.498903e-03 -2.557623e-04  1.518759e-04  2.362655e-03
#>  [91] -1.089266e-03 -6.927802e-04 -1.718232e-03  2.191668e-03  5.634669e-04
#>  [96]  1.721026e-03  1.553591e-03 -3.766823e-05  1.833268e-03  3.442976e-04
#> 
#> [[3]]
#>   [1] -9.484370e-02 -1.284164e-02 -4.088803e-03 -3.672588e-03 -1.695742e-04
#>   [6] -1.708157e-03  1.975601e-03  1.229795e-03 -1.776184e-03  2.301582e-03
#>  [11]  2.148805e-04  1.304673e-03 -7.614513e-04 -1.125296e-03 -1.302509e-03
#>  [16] -5.120238e-04  1.884588e-04  5.976523e-04  1.539508e-03  6.308994e-04
#>  [21]  2.079912e-03 -1.888892e-04 -4.110778e-04  1.541191e-04 -1.276704e-03
#>  [26] -3.368428e-03 -2.657639e-04  2.547245e-05 -1.001576e-03  5.158864e-04
#>  [31] -8.598523e-04  3.355713e-03  9.919494e-04  1.208061e-03  1.409711e-03
#>  [36]  1.662640e-03  3.924411e-04  1.625773e-03  1.300101e-03 -1.420258e-03
#>  [41]  1.091149e-04  1.596124e-03 -6.434646e-05 -1.658803e-04 -2.458295e-04
#>  [46]  6.819346e-04  2.083511e-03 -5.206368e-04  1.538927e-03 -2.626486e-03
#>  [51]  3.976998e-03  1.284361e-03  9.493542e-04 -3.980565e-04 -6.979055e-04
#>  [56] -2.051027e-04  6.285126e-04  2.151048e-03 -1.426173e-03  8.747231e-04
#>  [61]  3.524664e-04 -1.190936e-03 -2.799227e-04 -8.936098e-04  1.434325e-03
#>  [66]  2.127943e-03  3.961711e-04  4.609820e-04 -8.259977e-05  1.371883e-03
#>  [71] -3.130622e-04 -1.715346e-04 -2.479200e-03 -3.377097e-03 -1.161349e-03
#>  [76] -1.225680e-03 -5.630998e-05  1.418771e-03  7.632455e-04  2.137508e-03
#>  [81] -6.094437e-04 -2.173850e-03 -1.628787e-03 -1.648740e-04 -1.211295e-03
#>  [86] -5.640841e-04  2.039230e-03  2.460614e-03 -2.858428e-04 -6.077265e-04
#>  [91]  5.749371e-05  1.757643e-03  6.942409e-04  1.751259e-03  1.747847e-03
#>  [96] -6.848935e-04 -9.715990e-04  1.499377e-03 -8.642548e-04 -2.132413e-03
example$N.Gx
#> [[1]]
#> [1] 5000
#> 
#> [[2]]
#> [1] 4000
#> 
#> [[3]]
#> [1] 6000
example$Geno[[1]][1:5,1:5]
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    0    0    0    1    1
#> [2,]    0    0    0    1    1
#> [3,]    1    0    0    0    0
#> [4,]    0    0    0    0    0
#> [5,]    0    0    0    0    0
example$Z
#>              z.1        z.2       z.3
#>   [1,] 0.1746261 0.14996633 0.1236590
#>   [2,] 0.2686071 0.17083853 0.2525722
#>   [3,] 0.2194457 0.19965372 0.2508865
#>   [4,] 0.2670951 0.22573726 0.1123414
#>   [5,] 0.2636908 0.17207661 0.1393634
#>   [6,] 0.2112726 0.25037016 0.2031992
#>   [7,] 0.2939146 0.27569587 0.1377040
#>   [8,] 0.1552494 0.09411183 0.2191566
#>   [9,] 0.1824191 0.28871891 0.1280374
#>  [10,] 0.1403017 0.20466488 0.1061048
#>  [11,] 0.1950397 0.11314435 0.2255481
#>  [12,] 0.2938074 0.17858719 0.2443251
#>  [13,] 0.2674978 0.29493786 0.2261661
#>  [14,] 0.2258577 0.16359818 0.2372788
#>  [15,] 0.2788787 0.19963010 0.1882092
#>  [16,] 0.1696446 0.25328636 0.2176982
#>  [17,] 0.1143128 0.29974026 0.1768178
#>  [18,] 0.1785005 0.24347494 0.1675308
#>  [19,] 0.1329873 0.26219060 0.2010904
#>  [20,] 0.1810242 0.30592834 0.2503564
#>  [21,] 0.2825436 0.19164803 0.1450465
#>  [22,] 0.2485199 0.10286553 0.1627199
#>  [23,] 0.2887422 0.10848126 0.1139737
#>  [24,] 0.2919417 0.29323306 0.1919809
#>  [25,] 0.3073872 0.24761576 0.2004315
#>  [26,] 0.1315076 0.29608316 0.1587184
#>  [27,] 0.2796999 0.22950807 0.2179854
#>  [28,] 0.2615012 0.19712612 0.1890784
#>  [29,] 0.2384061 0.22308834 0.1289153
#>  [30,] 0.2302499 0.10420692 0.1619324
#>  [31,] 0.3054734 0.29208572 0.2240549
#>  [32,] 0.1283207 0.21841958 0.1054338
#>  [33,] 0.3049942 0.29757360 0.2509593
#>  [34,] 0.2880733 0.23436061 0.1961782
#>  [35,] 0.2807284 0.22994810 0.2082215
#>  [36,] 0.1418487 0.23201893 0.1206113
#>  [37,] 0.1443431 0.18477946 0.1594201
#>  [38,] 0.2966172 0.11048401 0.1244825
#>  [39,] 0.2231527 0.22939431 0.1587599
#>  [40,] 0.1489505 0.16952963 0.2119822
#>  [41,] 0.1811742 0.25775668 0.1224224
#>  [42,] 0.1176727 0.28675663 0.1064061
#>  [43,] 0.2718719 0.24273141 0.2275978
#>  [44,] 0.2824980 0.27286841 0.1289717
#>  [45,] 0.2806353 0.26295859 0.1684632
#>  [46,] 0.2082315 0.25131368 0.2478401
#>  [47,] 0.3078014 0.30094544 0.1586812
#>  [48,] 0.2646237 0.28137602 0.1710912
#>  [49,] 0.2228656 0.14555305 0.1272185
#>  [50,] 0.1674981 0.29300196 0.1409033
#>  [51,] 0.1791428 0.24439352 0.1133406
#>  [52,] 0.2686299 0.11893851 0.2155191
#>  [53,] 0.2341893 0.30782138 0.2337109
#>  [54,] 0.1047598 0.29041172 0.2391565
#>  [55,] 0.2385781 0.27551286 0.2244945
#>  [56,] 0.1833427 0.29201752 0.2392620
#>  [57,] 0.2781494 0.29470350 0.1582250
#>  [58,] 0.2711335 0.12855799 0.2038887
#>  [59,] 0.2560756 0.14738342 0.1261880
#>  [60,] 0.2370992 0.24412596 0.2492055
#>  [61,] 0.1523484 0.14663348 0.1238225
#>  [62,] 0.2758508 0.27130894 0.1246483
#>  [63,] 0.2612414 0.24775053 0.1241275
#>  [64,] 0.1274422 0.12664750 0.1277786
#>  [65,] 0.2558731 0.24907435 0.2363469
#>  [66,] 0.2342989 0.20960831 0.1771360
#>  [67,] 0.1552964 0.16521960 0.1519556
#>  [68,] 0.2092823 0.30757543 0.1700828
#>  [69,] 0.2228660 0.24297274 0.1630662
#>  [70,] 0.1299753 0.27415938 0.1892589
#>  [71,] 0.2484171 0.15585731 0.2145981
#>  [72,] 0.3054255 0.15750547 0.1548911
#>  [73,] 0.2375979 0.25777257 0.1638035
#>  [74,] 0.2778770 0.23477048 0.1046263
#>  [75,] 0.2471966 0.27702985 0.1537813
#>  [76,] 0.1995326 0.30390575 0.1316240
#>  [77,] 0.2021037 0.22712910 0.1310631
#>  [78,] 0.1447975 0.25919561 0.1032559
#>  [79,] 0.2417915 0.15719861 0.2056977
#>  [80,] 0.2314565 0.30153152 0.1946312
#>  [81,] 0.1692178 0.17458490 0.1072159
#>  [82,] 0.1125210 0.30062881 0.1576848
#>  [83,] 0.2355005 0.17652502 0.1433036
#>  [84,] 0.2091958 0.19695254 0.1753948
#>  [85,] 0.2411232 0.25452771 0.1052252
#>  [86,] 0.2216430 0.29110011 0.2149357
#>  [87,] 0.2918014 0.16676764 0.2404025
#>  [88,] 0.3055394 0.23306770 0.2188631
#>  [89,] 0.2892530 0.15372204 0.2119888
#>  [90,] 0.1069607 0.12512025 0.2297680
#>  [91,] 0.2222811 0.21884164 0.2379773
#>  [92,] 0.2338383 0.22762561 0.1880415
#>  [93,] 0.2019732 0.16874147 0.1516012
#>  [94,] 0.3017644 0.21681968 0.1048655
#>  [95,] 0.1977890 0.12115918 0.1096869
#>  [96,] 0.2765695 0.21401896 0.1127103
#>  [97,] 0.2330259 0.25968293 0.2297008
#>  [98,] 0.1745934 0.22664838 0.1885552
#>  [99,] 0.1748943 0.28734024 0.1846492
#> [100,] 0.2689979 0.30021472 0.2393025
```

``` r
fit = xtune_mJAM(betas.Gx = example$beta.gwas, N.Gx =example$N.Gx, Geno = example$Geno, Z = example$Z, c = 0.5)
#> Warning in optim(alpha.old, likelihood.alpha.theta.EN,
#> likelihood.alpha.theta.gradient.EN, : bounds can only be used with method L-
#> BFGS-B (or Brent)
fit$penalty.vector
#>                [,1]
#>   [1,] 4.424177e-05
#>   [2,] 4.424177e-05
#>   [3,] 4.424177e-05
#>   [4,] 4.424177e-05
#>   [5,] 4.424177e-05
#>   [6,] 4.424177e-05
#>   [7,] 4.424177e-05
#>   [8,] 4.424177e-05
#>   [9,] 4.424177e-05
#>  [10,] 4.424177e-05
#>  [11,] 4.424177e-05
#>  [12,] 4.424177e-05
#>  [13,] 4.424177e-05
#>  [14,] 4.424177e-05
#>  [15,] 4.424177e-05
#>  [16,] 4.424177e-05
#>  [17,] 4.424177e-05
#>  [18,] 4.424177e-05
#>  [19,] 4.424177e-05
#>  [20,] 4.424177e-05
#>  [21,] 4.424177e-05
#>  [22,] 4.424177e-05
#>  [23,] 4.424177e-05
#>  [24,] 4.424177e-05
#>  [25,] 4.424177e-05
#>  [26,] 4.424177e-05
#>  [27,] 4.424177e-05
#>  [28,] 4.424177e-05
#>  [29,] 4.424177e-05
#>  [30,] 4.424177e-05
#>  [31,] 4.424177e-05
#>  [32,] 4.424177e-05
#>  [33,] 4.424177e-05
#>  [34,] 4.424177e-05
#>  [35,] 4.424177e-05
#>  [36,] 4.424177e-05
#>  [37,] 4.424177e-05
#>  [38,] 4.424177e-05
#>  [39,] 4.424177e-05
#>  [40,] 4.424177e-05
#>  [41,] 4.424177e-05
#>  [42,] 4.424177e-05
#>  [43,] 4.424177e-05
#>  [44,] 4.424177e-05
#>  [45,] 4.424177e-05
#>  [46,] 4.424177e-05
#>  [47,] 4.424177e-05
#>  [48,] 4.424177e-05
#>  [49,] 4.424177e-05
#>  [50,] 4.424177e-05
#>  [51,] 4.424177e-05
#>  [52,] 4.424177e-05
#>  [53,] 4.424177e-05
#>  [54,] 4.424177e-05
#>  [55,] 4.424177e-05
#>  [56,] 4.424177e-05
#>  [57,] 4.424177e-05
#>  [58,] 4.424177e-05
#>  [59,] 4.424177e-05
#>  [60,] 4.424177e-05
#>  [61,] 4.424177e-05
#>  [62,] 4.424177e-05
#>  [63,] 4.424177e-05
#>  [64,] 4.424177e-05
#>  [65,] 4.424177e-05
#>  [66,] 4.424177e-05
#>  [67,] 4.424177e-05
#>  [68,] 4.424177e-05
#>  [69,] 4.424177e-05
#>  [70,] 4.424177e-05
#>  [71,] 4.424177e-05
#>  [72,] 4.424177e-05
#>  [73,] 4.424177e-05
#>  [74,] 4.424177e-05
#>  [75,] 4.424177e-05
#>  [76,] 4.424177e-05
#>  [77,] 4.424177e-05
#>  [78,] 4.424177e-05
#>  [79,] 4.424177e-05
#>  [80,] 4.424177e-05
#>  [81,] 4.424177e-05
#>  [82,] 4.424177e-05
#>  [83,] 4.424177e-05
#>  [84,] 4.424177e-05
#>  [85,] 4.424177e-05
#>  [86,] 4.424177e-05
#>  [87,] 4.424177e-05
#>  [88,] 4.424177e-05
#>  [89,] 4.424177e-05
#>  [90,] 4.424177e-05
#>  [91,] 4.424177e-05
#>  [92,] 4.424177e-05
#>  [93,] 4.424177e-05
#>  [94,] 4.424177e-05
#>  [95,] 4.424177e-05
#>  [96,] 4.424177e-05
#>  [97,] 4.424177e-05
#>  [98,] 4.424177e-05
#>  [99,] 4.424177e-05
#> [100,] 4.424177e-05
fit$beta.est
#> 101 x 1 sparse Matrix of class "dgCMatrix"
#>                        s1
#> (Intercept) -7.232828e-02
#> V1          -9.555932e-02
#> V2          -2.026847e-04
#> V3           2.451485e-04
#> V4          -8.101502e-05
#> V5          -2.311370e-04
#> V6          -8.790772e-06
#> V7          -9.261223e-05
#> V8           5.453051e-04
#> V9          -5.887649e-04
#> V10          1.709910e-03
#> V11         -1.245186e-03
#> V12          3.861045e-04
#> V13         -4.953119e-04
#> V14         -5.168105e-04
#> V15         -3.959064e-04
#> V16          2.318450e-04
#> V17          6.141960e-04
#> V18          3.313143e-05
#> V19         -1.143891e-03
#> V20          3.543663e-04
#> V21          6.358815e-05
#> V22         -7.487644e-04
#> V23         -4.626212e-04
#> V24          4.744777e-04
#> V25         -9.950295e-05
#> V26         -8.550765e-06
#> V27         -8.688538e-04
#> V28          3.552889e-04
#> V29         -1.853643e-03
#> V30          7.634342e-04
#> V31         -4.725890e-04
#> V32          4.169912e-04
#> V33         -5.447777e-05
#> V34          2.603979e-04
#> V35         -4.527358e-04
#> V36         -1.015047e-03
#> V37          1.723925e-04
#> V38         -8.746706e-04
#> V39          3.007656e-04
#> V40         -6.711283e-04
#> V41         -9.487773e-05
#> V42          6.795280e-05
#> V43         -6.832356e-04
#> V44          2.053594e-03
#> V45         -1.103191e-04
#> V46          9.740768e-05
#> V47         -6.682217e-04
#> V48          2.989473e-04
#> V49          7.831457e-04
#> V50         -1.720797e-03
#> V51          1.397615e-03
#> V52         -6.214250e-04
#> V53         -4.111019e-04
#> V54         -8.170270e-04
#> V55          1.229605e-03
#> V56          6.686432e-04
#> V57          1.504261e-04
#> V58          1.181285e-03
#> V59         -1.383388e-03
#> V60          2.847549e-04
#> V61         -8.011839e-04
#> V62         -1.145220e-03
#> V63          1.044623e-03
#> V64          6.201910e-04
#> V65         -1.686661e-05
#> V66         -6.759514e-04
#> V67          1.152577e-03
#> V68         -2.008094e-03
#> V69         -8.805059e-05
#> V70          3.517997e-04
#> V71          3.454151e-04
#> V72         -4.514944e-04
#> V73         -1.176341e-03
#> V74          1.023275e-04
#> V75          5.017557e-04
#> V76         -9.184806e-04
#> V77         -4.537675e-04
#> V78         -5.668372e-04
#> V79          1.672296e-04
#> V80          4.191265e-04
#> V81         -1.023650e-03
#> V82          9.312757e-04
#> V83         -7.640672e-04
#> V84         -1.631839e-03
#> V85          7.917452e-04
#> V86         -1.321113e-04
#> V87          3.337623e-04
#> V88         -7.465717e-04
#> V89         -6.238098e-04
#> V90          8.715812e-04
#> V91         -1.281768e-03
#> V92          9.988738e-04
#> V93          9.134489e-04
#> V94         -1.104635e-04
#> V95         -1.439194e-04
#> V96         -5.042664e-04
#> V97          2.871347e-04
#> V98          2.016502e-04
#> V99          1.513089e-04
#> V100        -7.800611e-04
```