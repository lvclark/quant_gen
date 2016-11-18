# Description

This repository contains a R functions that I created for estimating genetic correlations.  Comments in file headers describe what to provide for each argument.  Before using these functions, I suggest transforming the data if necessary to approximate a
normal distribution.

The function in genetic_correlation_traits.R estimates genetic correlations between traits measured
on the same field plots.  For example, if you had a collection of clonal genotypes or inbred lines, each was
planted multiple times across reps and/or locations, and on each plant you measured height, yield, and lodging
score, this function could estimate genetic correlation between each pair of traits.  The formula used is
equation 9 from Howe et al. (2000; doi:10.1007/s001220051525).

The function in genetic_correlation_sites.R estimates genetic correlations between environments within traits, i.e.
when the measurements were taken on different plants.  For example, if you had the same genotypes planted in replicated
field trials at two different locations and measured yield for each, the function would estimate genetic correlation
for yield between the two sites.  The formula used is equation 8 from Howe et al. (2000).