
The following is a list and description of the various data sets used in
the analyses.

\###Dataset 1

All individuals including the out group samples.Filtering is as
specified in [variant filtering](3_variant_filtering.md).

``` r
library(knitr)
```

    ## Warning: package 'knitr' was built under R version 4.4.3

``` r
data_frame <- data.frame(
  Location = c("Palmer Archipelago", "Low Island", "Bransfield Strait", "King George Island", "Elephant Island", "Outgroup (lineage 24)", "Total"),
  N = c(123, 14, 9, 9, 7, 2, 164)
)
kable(data_frame, caption = "Table 1: Individuals included in dataset 1.")
```

| Location              |   N |
|:----------------------|----:|
| Palmer Archipelago    | 123 |
| Low Island            |  14 |
| Bransfield Strait     |   9 |
| King George Island    |   9 |
| Elephant Island       |   7 |
| Outgroup (lineage 24) |   2 |
| Total                 | 164 |

Table 1: Individuals included in dataset 1.

\###Dataset 2

Dataset 2 was comprised only of individuals identified as being from the
same cluster (i.e. ancestral population), as there is suspicion other
clusters represent unidentified lineages, as cryptic speciation is
prolific within the *Doris “kerguelenensis”* species complex.The
individuals included are represented in table 2.

``` r
library(knitr)
data_frame <- data.frame(
  Location = c("Palmer Archipelago", "Low Island", "Bransfield Strait", "Elephant Island", "Total"),
  N = c(114, 14, 9, 7, 144)
)
kable(data_frame, caption = "Table 2: Individuals included in dataset 2.")
```

| Location           |   N |
|:-------------------|----:|
| Palmer Archipelago | 114 |
| Low Island         |  14 |
| Bransfield Strait  |   9 |
| Elephant Island    |   7 |
| Total              | 144 |

Table 2: Individuals included in dataset 2.

\###Dataset 3

Dataset 3 included all individuals except the two out group samples from
lineage 24. This dataset was created to replicate the results of dataset
1 and see if they would hold true.

``` r
library(knitr)
data_frame <- data.frame(
  Location = c("Palmer Archipelago", "Low Island", "Bransfield Strait", "King George Island", "Elephant Island", "Total"),
  N = c(123, 14, 9, 9, 7, 162)
)
kable(data_frame, caption = "Table 3: Individuals included in dataset 3.")
```

| Location           |   N |
|:-------------------|----:|
| Palmer Archipelago | 123 |
| Low Island         |  14 |
| Bransfield Strait  |   9 |
| King George Island |   9 |
| Elephant Island    |   7 |
| Total              | 162 |

Table 3: Individuals included in dataset 3.
