SUSY
----

SUSY computes synchrony as windowed cross-correlation based on two-dimensional time series, as described in Tschacher & Meier (2020).

[R package website](https://wtschacher.github.io/SUSY/)

----

Installation
----

```r
install.packages("SUSY", repos="https://wtschacher.github.io/SUSY/")
```

Usage
----

Note that the following example assumes that the source data are in a flat file and it has particular structure (column names in first row, whitespace as field separator, at least 5 columns). If you do not have such, then use the command in the comment below to mockup random data.

```r
library(SUSY)

## read in data from a flat file
data = read.csv(file.choose(), header=TRUE, sep=" ", na.strings=".")

## mockup random data if needed
#data = as.data.frame(replicate(5, runif(10000, 300, 330)))

## compute SUSY for column 2 and column 5
res = susy(data[, c(2, 5)], segment=30, Hz=15)
names(res)

## compute SUSY for columns 1-2 and 3-4
res = susy(data[, 1:4], segment=30, Hz=15)
names(res)

## print all SUSY computations
res

## subset (and print) susy object to single results
res[1]

## plot all SUSY computations, plot type 1
plot(res, type=1)

## plot only first SUSY computations, plot type 1 and 4
plot(res[1], type=c(1,4))

## plot only second SUSY computations, plot type 1, 2, 3, 4, 5
plot(res[2], type=1:5)

## compute SUSY for all permutations of columns
res = susy(data, segment=30, Hz=15, permutation=TRUE)
names(res)

## print legacy style
print(res, legacy=TRUE)

## export to flat file via data.frame and write.csv
df = as.data.frame(res)
df
```

[`susy` function manual](https://wtschacher.github.io/SUSY/reference/susy.html)
