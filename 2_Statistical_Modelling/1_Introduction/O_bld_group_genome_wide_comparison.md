## Calibration / model checking for the O blood group finding

In class we analysed this table:

```
                      non-O blood group    O blood group
controls              3419.733             3233.261
severe malaria cases  3925.334             2737.664
```

The odds ratio computed from this table is 0.74 - suggesting O blood group is protective.  And a test of this table using `chisq.test()` or `fisher.test()` shows that this is highly statistically significant.  That means:

> *If* the data were truly sampled from a binomial distribution with given frequency in each row
> and *if* the true odds ratio was zero
> then the table would almost never occur.

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
    "recessive_counts_matching_O_frequency.csv",
    header = T, # Tell R our file has a row of column names
    as.is = T   # Ask R to please not transmogrify any data
)
```

Some technical points on this are as follows.

