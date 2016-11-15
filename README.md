# Description

This repository contains an R function that I created for estimating genetic correlation between traits measured
on the same field plots.  For example, if you had a collection of clonal genotypes or inbred lines, each was
planted multiple times across reps and/or locations, and on each plant you measured height, yield, and lodging
score, this function could estimate genetic correlation between each pair of traits.  The formula used is
equation 9 from Howe et al. (2000; doi:10.1007/s001220051525).  Comments in the file header describe what to provide
for each argument.  Before using this function, I suggest transforming the data if necessary to approximate a
normal distribution.
