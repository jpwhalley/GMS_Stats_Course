## Calibration / model checking for the O blood group finding

In class we analysed this table:

```
                      non-O blood group    O blood group
controls              3419.733             3233.261
severe malaria cases  3925.334             2737.664
```

The odds ratio computed from this table is 0.74 - suggesting O blood group is protective.  And a test of this table using `chisq.test()` or `fisher.test()` shows that this is highly statistically significant.  That means:

* *If* the data were truly sampled from a binomial distribution with given frequency in each row
* and *if* the true odds ratio was zero
* *then* a table with such a large odds ratio would almost never occur.

These are approximations which don't really hold.  Samples are not really independent (e.g. population structure is itself a reflection of distant relationships); frequencies might vary across populations; in a finite population it's not really credible that the true effect size is really zero (though it might be very close).

We would therefore like to verify that there isn't something going on we haven't thought of.  I.e. is our null model really sensible?

One of the advantages of genome-wide analysis is that we have a lot of data to answer this with.
Suppose we compare our finding with all other 'similar' variants in the genome.  They have all been
genotyped on the same set of samples.  So they represent the same sampling biases (if any) and the same demography.  However, due to recombination they also represent many independent draws from the genealogical history of the sample.

The file `recessive_counts_matching_O_frequency.csv` contains data for all the alleles in our study
that match the frequency of the O blood group mutation.

Load it:

```
X = read.csv(
    "recessive_counts_matching_O_frequency.csv.gz",
    header = T,      # Tell R our file has a row of column names
    as.is = T,       # Ask R to please not transmogrify any data
    check.names = F  # Also ask R not to mess around with column names 
)
```
Let's first sanity check I've computed frequencies right.  It's often best to first collect numerical data into a matrix - this is a container of numbers that is layed out nicely in memory for computation.  (Unlike a data frame which holds columns of data of different types).  The `as.matrix()` function does this:
```
data = as.matrix( X[, c( "controls:A", "controls:B", "cases:A", "cases:B" )])
head(data) # take a look
X$frequency = (data[,2] + data[,4]) / rowSums(data)
hist( X$frequency )
abline( v = X$frequency[ which( X$rsid == 'rs8176719' )], col = 'red' )
```

Compute the log-odds ratio for every variant:
```
X$logOR = log( (data[,1] * data[,4]) / ( data[,2] * data[,3]) )
```

Compute the standard error of the log OR (using the standard formula):
```
# for table
#   a b
#   c d
# the formula is sqrt( 1/a + 1/b + 1/c + 1/d )
X$SE = sqrt( rowSums( 1/data ))
```

Calibrate against our null model.  This is done by finding the mass of the normal distribution with the given standard error, outside the estimated log odds ratio.
```
X$pvalue = pnorm( abs(X$logOR), sd = X$SE, lower.tail = T ) * 2
```
*Note*: the `abs()` and multiplying by two in the above occur because there's both a left and right tail for the normal distribution.  We'd be interested in effects going in either direction.

Now make a quantile-quantile plot.  This plots ranked p-values against ranked expected quantiles:
```
pvalues = log10( sort( X$pvalue ) )
expected = log10( (1:length(pvalues)  / (length(pvalues)+1)))
plot(
    expected,
    xlab = "Expected log10 P-value",
    pvalues,
    ylab = "Observed log10 P-value",
    pch = 19 # nice round dots
)
abline( a = 0, b = 1, lty = 2 )
# Let's highlight O blood group:
w = which( pvalues == X$pvalue[ X$rsid == 'rs8176719' ])
points(
    expected[w],
    pvalues[w],
    col = 'red'
)
text(
    expected[w],
    pvalues[w],
    side = 4,
    label = "O bld grp"
)
```

## Summary

We have shown

* That O blood has an estimated protective effect in these data.
* That such a large estimated effect is unlikely *under the formal model assumptions* if the effect were really zero.
* That this observed odds ratio doesn't seem to be a spurious of effect of something we haven't thought of (that affects other variants genome-wide).
* And indeed that genome-wide variants seem largely compatible with the model assumptions.

In addition we know

* The O blood group mutation (rs8176719) is a functional mutation - a deletion of protein-coding sequence that alters red cell surface antigens.

(At this point I am starting to get convinced.)
