DistanceDistributionAnalysis
============================

R scripts dedicated to statistical analysis of distribution of distances using various graphical representations 
and non-parametric tests.

Scripts available:

    # distri-difference-ks-test_example-1.R: adapted to a particular case
    # distri-difference-generic.R: available soon


# STRATEGY: considering the set of observed distances between 2 categories of genomic regions 
# (i.e., CAT1 and CAT2), we want to evaluate if they are genetically linked or associated 
# (biased towards being found closer than expected by chance). In other words, if the positions 
# of CAT1 are random with respect to CAT2, the distribution of distances [CAT1-CAT2] should 
# follow a uniform distribution.
#
#
# STEPS:
#
#   1) considers the observed distribution of N distances between CAT1 and CAT2
#
#   2) generates a random list of N distances in the same space of potential distances. 
#     The expected distribution should be uniform. 
#
#   3) creates a number of representations to visualise the 2 distributions profiles, including:
#     frequency histograms
#     density plots
#     cumulative frequency plots
#     boxplots
#     qqplots...
#   
#   4) tests whether the distributions are significantly different or not. 
#   Non-parametric methods are used here, considering that a normal distribution cannot be assumed 
#   for the data, especially as the second object is already known to be uniform.
#   Also determine the direction of the shift if there is one and if possible e.g determine the 
#   negative skew of the observed distances.
