## distri-difference-ks-test_example-1.R
##
## Author: Nouri BEN ZAKOUR
## Affiliations: University of Queensland
##

## PURPOSE: Statistical evaluation of the distribution bias of distances between RECs regions and MGEs

##
## ATTENTION: The following script is tailored to suit a particular example. A more generic script is also available. 
## Data used here is specific to the generation of Panel A to D of Figure SX in Petty, Ben Zakour et al. 2014 (submitted).
##


# 1) input list of N real REC-MGE distances between REC and closest MGE we observed in the object "obs"
# N=137 <=> 115 non overlapping REC regions + 22 overlapping REC regions
# real REC-MGE distances were calculated separately

obs <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,170801,186128,211084,242215,253014,303737,337461,350992,
         492157,511598,532941,538831,533200,450664,443420,440807,393183,391627,367958,359850,338173,309913,296565,211592,52655,
         60464,64254,74168,87883,100516,40005,28024,7059,608,21563,45859,74610,78027,171177,212373,126899,113694,95449,85393,18452
         ,7446,3972,13563,27646,28321,34496,35337,37459,59899,60935,61193,61483,61936,70444,71957,154479,156624,163400,128304,7223,
         10852,228755,237573,242053,221465,211032,194760,178056,150238,91789,111287,132009,136389,141773,311311,297727,291897,
         200720,183559,141463,139466,72490,72281,68988,64811,15622,59191,73425,104602,107151,122394,155598,162262,218966,225967,
         253884,255604,263169,267865,282022,288138,300438,248507,247666,205767,201703,152115,142401,135613,85017)


# 2) generate a random set of N distances within the range (0,M) corresponding to all the possible distances observable between a REC and a MGE.
# In other words, we consider the same MGE positions, and we just randomely select REC-MGE distances fitting in the intervals.
# The maximal value M corresponds to the maximal distance a randomly selected REC can be from a real MGE 
# i.e. M = max_MGE-MGE distance between 2 MGEs minus the min REC_size observed, then divided by 2

# In this data, the maximal distance observed is between GI-leuX and PHI1, i.e. 1089739 bp, 
# and the min_REC_size observed is 9 bp

max_MGE_MGE <- 1089739                                  # maximal distance observed between 2 MGEs
min_REC_size <- 9                                       # minimum size of a REC region
M <- as.integer((max_MGE_MGE - min_REC_size)/2)         # maximal distance a randomly selected REC can be from a real MGE
N <- 137                                                # number of REC regions                     
rand <- as.integer(runif(N,0,M))                        # generates a uniformally distributed list of integers


# Selected representations of the input data
# 3.1) Generates panel A "Frequency histogram"

hist(obs, col=rgb(0,0,1,0.5), xlab="Distance (bp)")


# 3.2) Generates panel B "Empirical Cumulative Distribution Function Plots for observed vs theoritical uniform distributions"

plot(ecdf(obs), xlim = range(c(obs,rand)), col=rgb(0,0,1,0.5), xlab = "Distance (bp)", ylab = "Cumulative frequency")
abline(0,1/M, col=rgb(1,0,0,0.5), lwd=3, lty=2)         # theoritical CDF


# 4) Statistical evaluation of the difference (and possibly direction of difference) between the (obs) and (rand) data

# 4.1) Two-sample Kolmogorov–Smirnov test (K-S test)

    # PRINCIPLE: The two-sample Kolmogorov-Smirnov test compares the cumulative distribution of two data sets, 
    # and computes a D value and a P value that depend on the largest discrepancy between distributions. 
    # It is used to test whether two samples are drawn from an identical continuous distributions.
    # The possible values i) "two.sided", ii) "less" and iii) "greater" of alternatives specify the null hypothesis that 
    # the true distribution function of x is i) equal to, ii) not less than or iii) not greater than the hypothesized distribution function 
    # (one-sample case) or the distribution function of y (two-sample case), respectively. 
    # ADVANTAGE: It is a non-parametric test that doesn't assume a particular distribution for the data 
    # except that it should be continuous.
    # CAUTION: ties are not handled, unless data is resampled.

# K-S test without data resampling for preliminary hypothesis testing
    # test 1: null hypothesis is rejected
    # test 2: confirms the alternative = "greater", which includes distributions for which x is stochastically smaller than y 
    # (in other words, the Cumulative Distribution Function of x lies above that of y)

ks.test(obs,rand)
ks.test(obs,rand, alternative="greater")


# 4.2) One-sample K-S test, comparing the observed distribution to the cumulative uniform distribution function ("punif"):
   
# K-S test with data resampling
       
R=9999                                                  # the number of replicates - 1

D.values = numeric(R)                                   # to store the D value results
P.values = numeric(R)                                   # to store the P value results
for(i in 1:R) { 
  obsS = sample(obs, size=N, replace=T)                 # obs data resampling
  randS = as.integer(runif(N,0,M))                      # random data resamling
  kstest <- ks.test(obsS,randS,alternative="greater")   # K-S test
  D.values[i] = kstest$statistic                        # D value stored
  P.values[i] = kstest$p.value                          # P value stored
}


# Generates panel C: "Distribution of P-value for 10,000 replicates"       
hist(D.values, breaks=30, xlab="D value")

# Determine mean and standard deviation for D value       
mean(D.values)  #mean
sd(D.values)    #sd


# Generates panel D: "Distribution of P-value for 10,000 replicates"
       # x-axis changed to a logarithmic scale (log 10)
       
breaks = c(0, 1e-16, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10, 1e-09, 1e-08, 1e-07, 1e-06, 1e-05, 1e-04, 1e-03, 1e-02)
Pvalues_hist <- hist(P.values,
     breaks = c(0, 1e-16, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10, 1e-09, 1e-08, 1e-07, 1e-06, 1e-05, 1e-04, 1e-03, 1e-02))
plot(Pvalues_hist$count, type='h', lwd=50, lend=1, lty=1, col="black", xlab="P value", ylab="", xaxt="n")
axis(1,at=0.5:15.5,label=breaks)
       
# Determine mean and standard deviation for P value  
mean(P.values)  #mean
sd(P.values)    #sd
