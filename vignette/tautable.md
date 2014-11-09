---
title: "tables of critical values for Kendall’s tau"
author: "Cory Brunson"
output: html_document
---

Rank correlation takes place at the interface between statistics and combinatorics. Calculating rank correlation coefficients such as Kendall's \(\tau\) and Spearman's \(\rho\) is straightforward and, for small lists, fast; and the distributions for long lists rapidly approach their asymptotics. Kendall's \(\tau\) mimics Pearson's \(r\) by ranging from \(-1\) (perfect negative correlation) through \(0\) (perfect null correlation) to \(+1\) (perfect positive correlation) while remaining nonparametric--unlike \(\rho\), it makes no assumptions about the underlying distribution of values (if such values even exist) that the rankings reflect. Within the \(\tau\) framework, any of several equivalent statistics may be used to describe the (dis)similarity in two rankings of \(n\) objects, and each is easily recovered from the others: \(P\), the number of *concordant* pairs, or pairs of objects ranked the same relative to each other in both lists; \(Q\), the number of *discordant* pairs, so that \(P+Q={n\choose 2}\); \(S=P-Q\); and \(\tau={S}/{n\choose 2}\). (Ties are ignored in this discussion.)

There are at least two extant and two apparently orphaned implementations of \(\tau\) in R: [the `kendall` option in the `cor` function] [1] (`stats` package) and the eponymous function in [the `Kendall` package] [2]; and `cor.fk` in [the `pcaPP` projection pursuit package] [3] and [the `kendall` method for the `rpucor` function](http://www.r-tutor.com/gpu-computing/correlation/kendall-rank-coefficient) in [the `rpud` NVIDIA package] [4]. While it's easy to compute \(\tau\) values and confidence intervals with these packages, though, i haven't come across a computational resource for critical values, or quartiles, of \(\tau\). That is, given \(n\) objects and a threshold \(p\), what is the level of concordance (or, equivalently, discordance) between two uniformly random rankings that occurs with probability less than \(p\), and what value \(\tau_{n,p}\) corresponds to this concordance? [The normal asymptote \(\tau_{n,p}\approx\frac{\sqrt{2(2n+5)}}{3\sqrt{n(n-1)}}z_p\) (where \(Z\sim N(0,1)\)) offers an approximate solution for large values of \(n\)](http://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/kend_tau.htm), but there will always be a place for exact answers if they're available.

[1]: http://stat.ethz.ch/R-manual/R-patched/library/stats/html/cor.html
[2]: http://cran.r-project.org/web/packages/Kendall/Kendall.pdf
[3]: http://cran.r-project.org/web/packages/pcaPP/pcaPP.pdf
[4]: http://www.r-tutor.com/content/download



To find a specific critical value using any of these implementations of \(\tau\), since these functions do not receive values of \(\tau\), \(P\), etc. directly, one must gain control of the number of concordant pairs in the input rankings. This amounts to controlling the number of inversions \(k\) in a permutation \(\pi\) of \(n\), taking the identity permutation \((1,2,\ldots,n)\) as the second ranking:


```r
# Count inversions in a vector
vec.inv <- function(vec) {
    sum(unlist(sapply(1:(length(vec) - 1), function(i) {
        sapply((i + 1):length(vec), function(j) vec[i] > vec[j])
    })))
}
```

Whereas the identity has zero inversions, \(k\) here plays the role of \(Q\). (And, since \(\tau\) is symmetric, it's enough to consider \(k\) values from \(0\) to \(\left\lceil{n\choose 2}/{2}\right\rceil\).) To find such a permutation, recall (or [find out for yourself](), or just take my word for it) that the number of inversions of \(\pi\) equals the minimal length of a word for \(\pi\), and one minimal word for the longest permutation \(\pi_0=(n,n-1,\ldots,1)\) is
$$\prod_{j=1}^{n-1}\left(\prod_{i=j}^{1}s_i\right)=(s_1)(s_2s_1)\cdots(s_{n-1}s_{n-2}\cdots s_1)\text.$$
Therefore, a (highly non-arbitrary) permutation of length \(k\) obtains as the product of the first (rightmost) \(k\) generators in this word. Here is an implementation in R (with an option to leverage the symmetry for a bit of efficiency):


```r
# Generate a permutation of 1:n having a given number k of inversions
inv.vec <- function(n, k, bw.ok = TRUE) {
    if(k > choose(n, 2) | k < 0 | any(c(n, k) %% 1 != 0)) stop('invalid n or k')
    backwards <- (k > choose(n, 2) / 2) & bw.ok
    if(backwards) k <- choose(n, 2) - k
    vec <- 1:n; ninv <- 0
    while(ninv < k) {
        i <- vec[1]
        j <- if(k - ninv < n - i) (2 + k - ninv):n else {
            if(i == 1) c() else (2 + n - i):n
        }
        vec <- c(vec[-c(1, j)], i, vec[j])
        ninv <- ninv + {
            if(k - ninv < n - i) (k - ninv) else (n - i)
        }
        stopifnot(ninv == vec.inv(vec))
    }
    if(backwards) vec <- rev(vec)
    return(vec)
}
```

From there it is quick work to whittle down the options. For example, if i want to know the "value to beat" for a type-I error \(p<.01\) on a one-sided test for positive correlation between two rankings of \(n=30\) objects using the `Kendall` package, i can run the following:


```r
n <- 30; p <- .01
# range of revelant k (at most half of all possible discordant pairs)
k.ran <- c(0, floor(choose(n, 2) / 2))
# corresponding one-sided p-values
p.ran <- sapply(k.ran, function(k) Kendall(x = 1:n, y = inv.vec(n, k))$sl / 2)
# interpolate (linearly)
repeat {
    if(p <= p.ran[1]) {
        k <- NA; break
    }
    if(p > p.ran[2]) {
        k <- k.ran[2]; break
    }
    q <- (p - p.ran[1]) / diff(p.ran)
    k <- round(q * diff(k.ran) + k.ran[1])
    if(k == k.ran[1]) k <- k + 1 else if(k == k.ran[2]) k <- k - 1
    k.p <- Kendall(x = 1:n, y = inv.vec(n, k))$sl / 2
    if(k.p == p) break else if(k.p < p) {
        k.ran[1] <- k; p.ran[1] <- k.p
    } else if(k.p > p) {
        k.ran[2] <- k; p.ran[2] <- k.p
    }
    if(diff(k.ran) < 2) {
        k <- k.ran[1]; break
    }
}
tau.K  <- round((choose(n, 2) - 2 * k) / choose(n, 2), 4)
print(paste('need tau >', tau.K))
```

```
## [1] "need tau > 0.3057"
```



This result disagrees slightly with the value 0.301 in [Peter M Lee's critical values table](http://www.york.ac.uk/depts/maths/tables/kendall.pdf), possibly due to error in `Kendall`, which cites [this algorithm by Best and Gipps](http://www.jstor.org/stable/2347062). The same procedure with `Kendall(...)$sl` exchanged for `cor.test(..., alternative = 't', method = 'kendall')$p.value`, which is exact for small \(n\) and computed for the asymptotic Z-score otherwise, yields a critical value of 0.3011.

All this, however, makes for a very inefficient means to generate the titular table of critical values of \(\tau\). For that, rather than constructing a permutation of a desired length having a desired number of inversions, then iterating this construction over many lengths and inversion counts, it makes sense, provided it is possible, both

a. to produce the distribution of the inversion counts \((0,1,\ldots,{n\choose 2})\) for each length \(n\), and
b. to compute critical values along the (very probably) recursive process that yields these distributions, rather than producing each distribution from scratch.

Conveniently, Shashank at StackOverflow [has already made the key observation](http://stackoverflow.com/questions/19372991/number-of-n-element-permutations-with-exactly-k-inversions) that these distributions are the [Mahonian numbers](http://oeis.org/A008302) \(T(n,k)\), whose generating functions take the elegant form:
\[\prod_{j=0}^{n-1}\left(\sum_{i=0}^{j}x^i\right)=\sum_{k=0}^{{n\choose 2}}T_{n,k}x^k\text,\]
which implies that they satisfy a simple recursive formula:


```r
# Compute the Mahonian numbers
mahonians.recursive <- function(n) {
    if(n <= 0 | n%%1 != 0) stop('n must be a positive integer')
    if(n == 1) 1 else rowSums(sapply(0:(n - 1), function(i) {
        c(rep(0, i), mahonians.recursive(n - 1), rep(0, n - 1 - i))
    }))
}
print(paste('Example: the Mahonian numbers for n = 5 are',
            toString(mahonians.recursive(5))))
```

```
## [1] "Example: the Mahonian numbers for n = 5 are 1, 4, 9, 15, 20, 22, 20, 15, 9, 4, 1"
```

This formula can now serve as a basis for the construction of a critical values table for \(\tau\) with whatever maximum length \(n\) and suite of \(p\)-value thresholds one desires; though the external recursion slows it down dramatically, so for the implementation an internal loop is used. (For convenience, the code below produces a table of critical values of \(k\), from which those of \(\tau\) can be recovered by arithmetic.)


```r
# Rapid table of critical values (n a vector, alpha a vector)
tau.crit.table <- function (n, alpha, incl.len = TRUE, stat = "tau") {
    if (any(n <= 0 | n%%1 != 0)) 
        stop("n must be a positive integer")
    if (any(alpha > 0.5 | alpha <= 0)) 
        stop("alpha must be positive and below .5")
    vec <- 1
    m <- 1
    max.n <- max(n)
    mat <- matrix(NA, nr = if (min(n) == 1) 
        1
    else 0, nc = length(alpha))
    while (max.n > m) {
        vec <- rowSums(sapply(0:m, function(i) {
            c(rep(0, i), vec, rep(0, m - i))
        }))
        m <- m + 1
        if (m %in% n) {
            s <- cumsum(vec[1:floor(choose(m, 2)/2)])/sum(vec)
            i <- sapply(alpha, function(a) if (s[1] < a) 
                max(which(s < a)) - 1
            else NA)
            mat <- rbind(mat, i)
        }
    }
    stopifnot(length(vec) == choose(max.n, 2) + 1)
    colnames(mat) <- alpha
    rownames(mat) <- n
    if (stat == "tau") {
        mat <- (choose(n, 2) - 2 * mat)/choose(n, 2)
    }
    else if (stat == "P" | grepl("^concord", stat, ignore.case = TRUE)) {
        mat <- choose(n, 2) - mat
    }
    else if (stat == "S") {
        mat <- choose(n, 2) - 2 * mat
    }
    else if (!(stat %in% c("K", "k", "Q")) & !grepl("^discord|^inv", 
        stat, ignore.case = TRUE)) {
        stop("unknown statistic")
    }
    if (incl.len) 
        mat <- cbind(n = n, mat)
    return(mat)
}
```

As an example, here are the critical values (of \(k\) and of \(\tau\)) for lengths at intervals of 10 up to 100 and type-I error thresholds at negative powers of 10:


```r
# Ranges of n and of alpha
example.n <- seq(10, 100, 10); example.alpha <- 10 ^ (-1:-4)
# Critical value table for k
inv.crit.vals <- tau.crit.table(n = example.n,
                              alpha = example.alpha,
                              incl.len = TRUE, stat = 'inv')
print(inv.crit.vals)
```

```
##       n  0.1 0.01 0.001 1e-04
## 10   10   14    9     5     3
## 20   20   74   59    48    40
## 30   30  180  152   132   116
## 40   40  334  290   258   233
## 50   50  535  473   429   394
## 60   60  783  702   644   597
## 70   70 1080  978   904   845
## 80   80 1425 1300  1210  1137
## 90   90 1817 1669  1561  1474
## 100 100 2259 2084  1958  1856
```

```r
# Critical value table for tau
tau.crit.vals <- tau.crit.table(n = example.n,
                              alpha = example.alpha,
                              incl.len = TRUE)
print(tau.crit.vals)
```

```
##       n     0.1   0.01  0.001  1e-04
## 10   10 0.37778 0.6000 0.7778 0.8667
## 20   20 0.22105 0.3789 0.4947 0.5789
## 30   30 0.17241 0.3011 0.3931 0.4667
## 40   40 0.14359 0.2564 0.3385 0.4026
## 50   50 0.12653 0.2278 0.2996 0.3567
## 60   60 0.11525 0.2068 0.2723 0.3254
## 70   70 0.10559 0.1901 0.2513 0.3002
## 80   80 0.09810 0.1772 0.2342 0.2804
## 90   90 0.09263 0.1665 0.2205 0.2639
## 100 100 0.08727 0.1580 0.2089 0.2501
```

The agreement of several rows with their analogs in [Lee's table](http://www.york.ac.uk/depts/maths/tables/kendall.pdf) validates the procedure; other tables exist for others of the equivalent statistics i mentioned above, e.g. [the number of concordant pairs](http://faculty.washington.edu/heagerty/Books/Biostatistics/TABLES/Kendall.pdf) and [the numerator \(S\)](http://books.google.com/books?id=IMbVyKoZRh8C&pg=PA691#v=onepage&q&f=false), which the interested reader can check for themselves. Most of these functions, and a few more, are available [in package form](https://github.com/corybrunson/tautable), in accordance with [the Leek Group's advice](https://github.com/jtleek/rpackages#when-to-start-writing-an-r-package).

As a final test, here are the runtimes for tables up to maximum \(n\) values of 25 through 200 at intervals of 25, performed on a 3.2 GHz Intel Core i5, which hint at the slow combinatorial explosion that eventually undermines the usefulness of the method, though far beyond the point at which the normal approximaiton becomes adequate:


```r
N <- seq(25, 200, 25)
st.mat <- sapply(N, function(n) {
    system.time(tau.crit.table(n = n,
                               alpha = 10^(-1:-3),
                               incl.len = TRUE))
})
plot(x = N, y = st.mat[3, ], type = 'b', log = 'y',
     xlab = 'Maximum number of objects', ylab = 'Time elapsed')
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9.png) 
