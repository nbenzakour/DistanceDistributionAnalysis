## distri-difference-generic.R
##
## Author: Nouri BEN ZAKOUR
## Affiliations: University of Queensland
##


# STRATEGY: considering the set of observed distances between 2 categories of genomic regions (i.e., CAT1 and CAT2), 
# we want to evaluate if they are genetically linked or associated (biased towards being found closer than expected by chance). 
# In other words, if the positions of CAT1 are random with respect to CAT2, the distribution of distances [CAT1-CAT2] 
# should follow a uniform distribution.

# CUSTOMISATION: while application here is targeted towards comparing an observed distribution of distances vs a random one, 
# it is suitable to compare 2 distincts observed distributions of any kind (same number of items for both lists). 
# Depending on the type of distribution, it is the user responsibility to determine which tests available in this script 
# are suitable or not.


# STEPS:
  
  # 1) considers the distribution of n observed distances 

  # 2) generates a random list of n theoretical distances in the same space of potential distances. The expected distribution should be uniform.
  # This step is required to make visual exploratory comparisons between the observed distribution and a random uniform one.
  # In optimal conditions, subsequent statistical tests should resample the uniform distribution.

  # 3) creates a number of representations to visualise the 2 distributions profiles, including: 
      #     frequency histograms
      #     density plots
      #     cumulative frequency plots
      #     boxplots
  # Graphical representations can be performed using either the observed distribution only, or in comparison to a second distribution
  # (here the random distribution generated in step 2). This is particularly useful to visualise what trend the observed data follows.

  # 4) tests whether the distributions are significantly different or not. Non-parametric methods are used here, considering that 
  # a normal distribution cannot be assumed for the data, especially as the second object is already known to be uniform. 
  # Depending on the data, tests results might not be meaningful. Tests proposed in this script are:
      #     One- and two-sampled Kolmogorov-Smirnov (K-S) tests
      #     One-sampled K-S bootstrapped test
      #     Independent 2-group Mann-Whitney U Test
      #     Chi-square test


#### 1) considers the observed distribution of N distances of 2 categories of genomic regions (i.e., CAT1 and CAT2)

# obs <- read.csv("~/path/to/your/file/observed_distances.csv", head=TRUE) # customize to suit user data

obs <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,170801,186128,211084,242215,253014,303737,337461,350992,
         492157,511598,532941,538831,533200,450664,443420,440807,393183,391627,367958,359850,338173,309913,296565,211592,52655,
         60464,64254,74168,87883,100516,40005,28024,7059,608,21563,45859,74610,78027,171177,212373,126899,113694,95449,85393,18452
         ,7446,3972,13563,27646,28321,34496,35337,37459,59899,60935,61193,61483,61936,70444,71957,154479,156624,163400,128304,7223,
         10852,228755,237573,242053,221465,211032,194760,178056,150238,91789,111287,132009,136389,141773,311311,297727,291897,
         200720,183559,141463,139466,72490,72281,68988,64811,15622,59191,73425,104602,107151,122394,155598,162262,218966,225967,
         253884,255604,263169,267865,282022,288138,300438,248507,247666,205767,201703,152115,142401,135613,85017)


#### 2) generates a random list of n theoretical distances in the same space of potential [CAT1-CAT2] distances.
      # We consider only variation regarding the position of CAT1 elements with respect to CAT2 regions wich are fixed.
      # The maximal value M corresponds to the maximal distance max_CAT1_CAT2 a randomly selected CAT1 region can be from a real CAT2 region.


# In this data, the maximal distance observed is between GI-leuX and PHI1, i.e. 1089739 bp, 
# and the min_REC_size observed is 9 bp

max_CAT2_CAT2 <- 1089739                                  # maximal distance observed between 2 CAT2 regions
min_CAT1_size <- 9                                        # minimum size of a CAT1 region
M <- as.integer((max_CAT2_CAT2 - min_CAT1_size)/2)        # maximal distance a randomly selected CAT1 can be from fixed CAT2
N <- 137                                                  # number of REC regions                     
rand <- as.integer(runif(N,0,M))                          # generates a uniformally distributed list of integers


#### 3) selected representations of the input data
# 3.1) Frequency histograms 

hist(obs, col=rgb(0,0,1,0.5))                             # observed distri (blue)
hist(rand, add=T, col=rgb(1,0,0,0.5))                     # random distri (red)

# 3.2) Density plots

plot(density(obs))                                        # observed distrib (black) 
lines(density(rand), col=2)                               # random distrib (red)

# 3.3) Boxplots

obs_rand <-cbind(obs,rand)                                # binding observed distri and random distri in a 2-col table
boxplot(obs_rand, horizontal=TRUE)

# 3.4) qqplots

qqplot(obs,rand)

# 3.5) Empirical Cumulative Frequency Plots for: 
      # 3.5.1) observed distrib and random distrib

plot(ecdf(obs), xlim = range(c(obs,rand)), col=rgb(0,0,1,0.5), main = 'Empirical Cumulative Distribution obs vs rand', xlab = "Distance (bp)", 
     ylab = "Cumulative frequency") # in blue
plot(ecdf(rand), add=T, lty="dashed", col=rgb(1,0,0,0.5)) # in red
      # with this data CDF for (obs) is to the "left" of CDF for (rand), tested further for significance in 4)

      # 3.5.2) observed distrib and theoretical trend

plot(ecdf(obs), xlim = range(c(obs,rand)), col=rgb(0,0,1,0.5), main = 'Empirical Cumulative Distribution obs vs theoretical',
     xlab = "Distance (bp)", ylab = "Cumulative frequency")
abline(0,1/M, col=rgb(1,0,0,0.5), lwd=3, lty=2)           # theoritical CDF


#### 4) Statistical evaluation of the difference (and possibly direction of difference) between the (obs) and (rand) data
# 4.1) Two-sample Kolmogorov–Smirnov test (K-S test)

# PRINCIPLE: The Kolmogorov-Smirnov test compares the cumulative distribution of the two data sets, 
# and computes a D value and a P value that depend on the largest discrepancy between distributions. 
# It is used to test whether two samples are drawn from identical continuous distributions.
# The possible values i) "two.sided", ii) "less" and iii) "greater" of alternatives specify the null hypothesis that 
# the true distribution function of x is i) equal to, ii) not less than or iii) not greater than the hypothesized distribution function 
# (one-sample case) or the distribution function of y (two-sample case), respectively. 
# It is a non-parametric test that doesn't assume a particular distribution for the data except that it should be continuous.
# CAUTION: ties are not handled, unless data is resampled.

# K-S test without data resampling for preliminary hypothesis testing

ks.test(obs,rand)
      # results: by default alternative="two.sided". Here the p-value is low, null hypothesis is rejected 
      # (= distribution functions are different).

# The possible values i) "two.sided", ii) "less" and iii) "greater" of alternative specify the null hypothesis that the true distribution 
# function of x is i) equal to, ii) not less than or iii) not greater than the hypothesized distribution function (one-sample case) 
# or the distribution function of y (two-sample case), respectively. 

ks.test(obs,rand, alternative="greater") 
      # results: in the two-sample case alternative = "greater" includes distributions for which x is stochastically smaller than y 
      # (in other words, the CDF of x lies above that of y)

ks.test(obs,rand, alternative="less")
      # when alternative="less" => p-value > 0.05

### POTENTIAL CAVEAT: if there are ties in the data, depending on the random distribution values, the p-value might not be reliable
### SOLUTION: use the boostrapped version of the K-S to avoid the negative impact of ties (see 4.2 and 4.3)


# 4.2) One-sample K-S test, comparing the observed distribution to the cumulative uniform distribution function ("punif"):
      # IMPORTANT: Additional statistical support is provided here multiple iterations of the K-S test with data resampling

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


# Frequency distribution of P-value for 10,000 replicates       
hist(D.values, breaks=30, xlab="D value")

# Determine mean and standard deviation for D value       
mean(D.values)  #mean
sd(D.values)    #sd


# Frequency distribution of P-value for 10,000 replicates
# x-axis changed to a logarithmic scale (log 10)

breaks = c(0, 1e-16, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10, 1e-09, 1e-08, 1e-07, 1e-06, 1e-05, 1e-04, 1e-03, 1e-02)
Pvalues_hist <- hist(P.values,
                     breaks = c(0, 1e-16, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10, 1e-09, 1e-08, 1e-07, 1e-06, 1e-05, 1e-04, 1e-03, 1e-02))
plot(Pvalues_hist$count, type='h', lwd=50, lend=1, lty=1, col="black", xlab="P value", ylab="", xaxt="n")
axis(1,at=0.5:15.5,label=breaks)

# Determine mean and standard deviation for P value  
mean(P.values)  #mean
sd(P.values)    #sd


# 4.3) Bootstrap Kolmogorov-Smirnov
      # Requires the Matching package to run.
      # This function executes a bootstrap version of the univariate Kolmogorov-Smirnov test which provides correct coverage even when 
      # the distributions being compared are not entirely continuous. Ties are allowed with this test unlike the traditional K-S test.

install.packages("Matching", dependencies=TRUE)

ks.boot(obs,rand, nboots=1000)
ks <- ks.boot(obs,rand, nboots=1000)
summary(ks)
ksg <- ks.boot(obs,rand, nboots=1000, alternative="greater")
summary(ksg)

      # NOTE: ks.boot connot be used as one-sampled vs the actual CDF "punif" (contrary to ks.test),
      # as the second argument needs to be a numeric vector of data values 


# 4.4) Independent 2-group Mann-Whitney U Test
      # Like K-S, it also a non-parametric test.
      # PRINCIPLE: test first ranks all the values from low to high, and then computes a P value that depends on the discrepancy 
      # between the mean ranks of the two groups. wilcox.test(y,x) where y and x are numeric
      # In the Mann-Whitney-Wilcoxon Test, we can decide whether the population distributions are identical without assuming
      # them to follow the normal distribution. MW handles ties.
      # Here, the null hypothesis is that the 2 variables have identical data distribution
      # if p-value < 0.05 significance level, the null hypothesis is rejected meanning that the 2 distributions come from
      # non-identical populations.
      # CAUTION: Distributions should have similar shapes for the test to be interpretable. 
      # In the current case, the MW test might not be appropriate. 

wilcox.test(obs,rand)
wilcox.test(obs,rand,alternative="less")              # alternative hypothesis greater and less are inverted compared to K-S test 


# 4.6) Chi-test
# Chi-test are probably not the best tests to use here but is provided as an option

chisq.test(obs_rand)                                  # on 2-col table with obs and rand generated in 3.4

freq_obs_rand <- cbind(freqobs$counts,freqrand$counts)  # binding observed freq distri and random freq distri in a 2-col table
chisq.test(freq_obs_rand, sim=TRUE)                   # on 2-col table with obs and rand frequencies
