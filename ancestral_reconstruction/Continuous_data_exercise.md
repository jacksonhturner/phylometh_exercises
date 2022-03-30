Continuous Data Exercise
================
Jackson Turner
3/29/2022

``` r
library(ape)
library(geiger)
library(phytools)
```

    ## Loading required package: maps

``` r
library(OUwie)
```

    ## Loading required package: corpcor

    ## Loading required package: nloptr

    ## Loading required package: RColorBrewer

``` r
library(akima)
```

This exercise uses largely unpublished Co1 sequences from my PI,
Dr. Kevin Moulton. The tree I reference is a simple, one-gene phylogeny
of meniscus midges (Diptera:Dixidae) constructed with the Co1 genes of
several species. The Co1 sequences used to construct this tree were
aligned with MAFFT and this tree was built using RAxML. I personally
gathered the discrete data from these taxa and the continuous data was
made up for the explicit purpose of this exercise. While an attempt was
made to make the continuous data seem realistic and like what one would
expect if they were to actually measure this data, it is fake and has no
bearing outside this exercise.

Let’s load and clean our data and see how our tree changes with a
continuous character.

``` r
tree <- read.tree("C:/Users/jturn/OneDrive/Documents/Skool/EEB 587/phylometh_exercises/continuous_data/RAxML_bestTree.T2.tre")
continuous_data<-read.csv("C:/Users/jturn/OneDrive/Documents/Skool/EEB 587/phylometh_exercises/continuous_data/Co1_continuous.csv",stringsAsFactors=FALSE)
rownames(continuous_data)<-continuous_data[,1]
continuous_data<-continuous_data[,-1]
tree<-multi2di(tree)
tree<-chronos(tree) #to make ultrametric
```

    ## 
    ## Setting initial dates...
    ## Fitting in progress... get a first set of estimates
    ##          (Penalised) log-lik = -5.842152 
    ## Optimising rates... dates... -5.842152 
    ## Optimising rates... dates... -5.799642 
    ## 
    ## log-Lik = -5.763155 
    ## PHIIC = 123.54

``` r
cleaned_continuous<-treedata(tree,continuous_data,sort=TRUE)
data_vec <- cleaned_continuous$data[,1]
names(data_vec) <- rownames(cleaned_continuous$data)
contMap(tree,data_vec)
```

![](Continuous_data_exercise_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

We’ll use a Brownian Motion model to determine how our continuous trait
changes through time, and we’ll evaluate its AIC score compared to an OU
model.

``` r
BM1<-fitContinuous(tree, cleaned_continuous$data, model="BM")
cat("The rate of evolution is", BM1$tailpipe_length_width_th$opt$sigsq, "in units of","um^2 per length of time.","\n")
```

    ## The rate of evolution is 14.17226 in units of um^2 per length of time.

``` r
OU1 <- fitContinuous(tree, cleaned_continuous$data, model="OU")
par(mfrow=(c(1,2)))
plot(tree, show.tip.label=FALSE)
plot(tree)
```

![](Continuous_data_exercise_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
cat("The AIC score for BM1 is", BM1$tailpipe_length_width$opt$aic, "\n")
```

    ## The AIC score for BM1 is 98.93725

``` r
cat("The AIC score for OU1 is", OU1$tailpipe_length_width$opt$aic, "\n")
```

    ## The AIC score for OU1 is 100.9125

``` r
cat("The delta AIC for OU1 is", OU1$tailpipe_length_width$opt$aic-BM1$tailpipe_length_width$opt$aic, "\n")
```

    ## The delta AIC for OU1 is 1.975205

``` r
cat("The delta AIC for BM1 is", BM1$tailpipe_length_width$opt$aic-BM1$tailpipe_length_width$opt$aic, "\n")
```

    ## The delta AIC for BM1 is 0

For our data, the BM model works better than our OU model, as it has a
lower AIC score.

Let’s perform an ancestral reconstruction of our discrete traits,
rgroup, represnting a larval body type with dark scleritization and a
thick body. We’ll use the `ace` function to determine which ancestral
states to apply to the internal nodes of our tree.

Then, let’s create a tree with our new internal node labels.

``` r
discrete_data<-read.csv("C:/Users/jturn/OneDrive/Documents/Skool/EEB 587/phylometh_exercises/continuous_data/Co1_discrete.csv")
rownames(discrete_data)<-discrete_data[,1]
discrete_data<-discrete_data[,-1]
tree2<-tree
cleaned_discrete<-treedata(tree2,discrete_data,sort=TRUE)

barbs<-cleaned_discrete$data[,1]
meringo<-cleaned_discrete$data[,2]
rgroup<-cleaned_discrete$data[,3]
tailpipe_fused<-cleaned_discrete$data[,4]

one.discrete.char <- rgroup
reconstruction.info <- ace(one.discrete.char, tree, type="discrete", method="ML", CI=TRUE)
nodeBased.OUM.states <- colnames(reconstruction.info$lik.anc)[apply(reconstruction.info$lik.anc, 1, which.max)]

one_discrete_tree<-tree
one_discrete_tree$node.label<-vector()
for (i in 1:one_discrete_tree$Nnode){ 
  one_discrete_tree$node.label[i]<-nodeBased.OUM.states[i]
}
one_discrete_tree$node.label
```

    ##  [1] "0" "0" "1" "1" "1" "1" "1" "1" "1" "1" "1" "0" "0" "0" "0" "0" "0" "0" "0"

``` r
one_discrete_tree<-chronos(one_discrete_tree)
```

    ## 
    ## Setting initial dates...
    ## Fitting in progress... get a first set of estimates
    ##          (Penalised) log-lik = -17.33665 
    ## Optimising rates... dates... -17.33665 
    ## 
    ## log-Lik = -1e+100 
    ## PHIIC = 2e+100

``` r
nodecolr<-rep("green", length=one_discrete_tree$Nnode)
nodecolr[one_discrete_tree$node.label==1]<-"purple"
nodecolr[one_discrete_tree$node.label==0]<-"green"
plot(one_discrete_tree)
nodelabels(pch=16, col=nodecolr)
```

![](Continuous_data_exercise_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

Now, we’ll prepare to build our node-based BM and OU models, from which
we’ll extract corrected AIC scores to evaluate which model works best
for our data.

``` r
continuous_data_o<-read.csv("C:/Users/jturn/OneDrive/Documents/Skool/EEB 587/phylometh_exercises/continuous_data/Co1_continuous.csv",stringsAsFactors=FALSE)
check.identify(one_discrete_tree,continuous_data_o)
```

    ## The regime optima are unidentifiable.

    ## [1] 0

``` r
OU_data<-as.data.frame(continuous_data_o[,1])
OU_data$regime<-cleaned_discrete$data[,2]
OU_data$tailpipe_length_width<-cleaned_continuous$data[,1]
colnames(OU_data)<-c("Species","Reg","OU_data")

nodeBased.OUMV <- OUwie(one_discrete_tree, OU_data,model="OUMV", simmap.tree=FALSE,diagn=FALSE)
```

    ## Initializing... 
    ## Finished. Begin thorough search... 
    ## Finished. Summarizing results.

``` r
print(nodeBased.OUMV)
```

    ## 
    ## Fit
    ##        lnL      AIC     AICc      BIC model ntax
    ##  -51.94958 113.8992 118.1849 118.8778  OUMV   20
    ## 
    ## 
    ## Rates
    ##                  0         1
    ## alpha     2.167289  2.167289
    ## sigma.sq 47.443677 73.375890
    ## 
    ## Optima
    ##                  0         1
    ## estimate 11.921228 16.351558
    ## se        1.469015  2.887646
    ## 
    ## 
    ## Half life (another way of reporting alpha)
    ##         0         1 
    ## 0.3198223 0.3198223 
    ## 
    ## Arrived at a reliable solution

``` r
nodeBased.BM1 <- OUwie(one_discrete_tree, OU_data,model="BM1", simmap.tree=FALSE, diagn=FALSE)
```

    ## Initializing... 
    ## Finished. Begin thorough search... 
    ## Finished. Summarizing results.

``` r
nodeBased.BMS <- OUwie(one_discrete_tree, OU_data,model="BMS", simmap.tree=FALSE, diagn=FALSE)
```

    ## Initializing... 
    ## Finished. Begin thorough search... 
    ## Finished. Summarizing results.

``` r
nodeBased.OUM <- OUwie(one_discrete_tree, OU_data,model="OUM", simmap.tree=FALSE, diagn=FALSE)
```

    ## Initializing... 
    ## Finished. Begin thorough search... 
    ## Finished. Summarizing results.

``` r
nodeBased.OUMA <- OUwie(one_discrete_tree, OU_data,model="OUMA", simmap.tree=FALSE, diagn=FALSE)
```

    ## Initializing... 
    ## Finished. Begin thorough search... 
    ## Finished. Summarizing results.

``` r
nodeBased.OUMVA <- OUwie(one_discrete_tree, OU_data,model="OUMVA", simmap.tree=FALSE, diagn=FALSE)
```

    ## Initializing... 
    ## Finished. Begin thorough search... 
    ## Finished. Summarizing results.

Let’s determine lay out our corrected AIC scores so we can make a
decision on which model to use (Weirdly BM1 is the best for our data).
However, we’ll use our best OU model (OUM) since it, unlike our BM
models, has an alpha parameter to manipulate for this exercise.

``` r
modelAIC<-data.frame(cbind(nodeBased.BM1$AICc,nodeBased.BMS$AICc,nodeBased.OUM$AICc,nodeBased.OUMA$AICc,nodeBased.OUMV$AICc,nodeBased.OUMVA$AICc))
colnames(modelAIC)<-c("BM1","BMS","OUM","OUMA","OUMA","OUMVA")
rownames(modelAIC)<-"corrected AIC score:"
modelAIC #BM1 is nodeBased.OUM but using OUM since it's close and has an alpha parameter
```

    ##                           BM1      BMS      OUM     OUMA     OUMA    OUMVA
    ## corrected AIC score: 111.2331 113.9089 114.8043 116.0065 118.1849 120.1567

``` r
nodeBased.OUM
```

    ## 
    ## Fit
    ##       lnL      AIC     AICc      BIC model ntax
    ##  -52.0688 112.1376 114.8043 116.1205   OUM   20
    ## 
    ## 
    ## Rates
    ##                  0         1
    ## alpha     1.941445  1.941445
    ## sigma.sq 52.397632 52.397632
    ## 
    ## Optima
    ##                  0        1
    ## estimate 11.824077 16.62596
    ## se        1.659426  2.91117
    ## 
    ## 
    ## Half life (another way of reporting alpha)
    ##         0         1 
    ## 0.3570263 0.3570263 
    ## 
    ## Arrived at a reliable solution

``` r
OUwie.fixed(one_discrete_tree,OU_data,model="OUM",alpha=nodeBased.OUM$solution[1,],sigma.sq=nodeBased.OUM$solution[2,],theta=nodeBased.OUM$theta[,1])
```

    ## Calculating likelihood using fixed parameter values: 1.941445 1.941445 52.39763 52.39763 11.82408 16.62596

    ## 
    ## Fit
    ##       lnL      AIC     AICc      BIC model ntax
    ##  -52.0688 112.1376 114.8043 116.1205   OUM   20
    ## 
    ## 
    ## Rates
    ##                  0         1
    ## alpha     1.941445  1.941445
    ## sigma.sq 52.397632 52.397632
    ## 
    ## Optima
    ##                 0        1
    ## estimate 11.82408 16.62596
    ## se             NA       NA

Now we’ll make a plot showing a plot of how our likelihood values vary
at given alpha values.

``` r
alpha.values<-seq(from=0.5,to=3,length.out=10)
likelihood.values <- rep(NA, length(alpha.values))

for (i in sequence(length(alpha.values))) {
  likelihood.values[i] <- OUwie.fixed(one_discrete_tree, OU_data, model="OUM", alpha=rep(alpha.values[i],2), sigma.sq=nodeBased.OUM$solution[2,], theta=nodeBased.OUM$theta[,1])$loglik
}
```

``` r
plot(x=alpha.values, y=likelihood.values, xlab="Alpha", ylab="Likelihood", type="l", bty="n")
points(x=nodeBased.OUM$solution[1,1], y=nodeBased.OUM$loglik, pch=16, col="red")
text(x=nodeBased.OUM$solution[1,1], y=nodeBased.OUM$loglik, "unconstrained OUM model", pos=2, col="red")
abline(h=max(likelihood.values)-2, lty="dotted")
```

![](Continuous_data_exercise_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

Let’s create a plot that shows how our theta parameters of our OU model
change with each other.

``` r
require("akima")
nreps<-400
theta1.points<-c(nodeBased.OUM$theta[1,1], rnorm(nreps-1, nodeBased.OUM$theta[1,1], 5*nodeBased.OUM$theta[1,2])) #center on optimal value, have extra variance
theta2.points<-c(nodeBased.OUM$theta[2,1], rnorm(nreps-1, nodeBased.OUM$theta[2,1], 5*nodeBased.OUM$theta[2,2])) #center on optimal value, have extra variance
likelihood.values<-rep(NA,nreps)

for (i in sequence(nreps)) {
  likelihood.values[i] <- OUwie.fixed(one_discrete_tree, OU_data, model="OUM", alpha=nodeBased.OUM$solution[1,], sigma.sq=nodeBased.OUM$solution[2,], theta=c(theta1.points[i], theta2.points[i]))$loglik
}

likelihood.differences<-(-(likelihood.values-max(likelihood.values)))

interpolated.points<-interp(x=theta1.points, y=theta2.points, z= likelihood.differences, linear=FALSE, extrap=TRUE, xo=seq(min(theta1.points), max(theta1.points), length = 400), yo=seq(min(theta2.points), max(theta2.points), length = 400))
```

``` r
contour(interpolated.points, xlim=range(c(theta1.points, theta2.points)),ylim=range(c(theta1.points, theta2.points)), xlab="Theta 1", ylab="Theta 2", levels=c(2,5,10),add=FALSE,lwd=1, bty="n", asp=1)

points(x=nodeBased.OUM$theta[1,1], y=nodeBased.OUM$theta[2,1], col="red", pch=16)

points(x=OU_data$OU_data[which(OU_data$Reg==1)],y=rep(min(c(theta1.points, theta2.points)), length(which(OU_data$Reg==1))), pch=18, col=rgb(0,0,0,.3)) #the tip values in regime 1, plotted along x axis
points(y=OU_data$OU_data[which(OU_data$Reg==2)],x=rep(min(c(theta1.points, theta2.points)), length(which(OU_data$Reg==2))), pch=18, col=rgb(0,0,0,.3)) #the tip values in regime 2, plotted along y axis
```

![](Continuous_data_exercise_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

Unfortunately, the below regime-painted phylogeny using a simmap fails
to produce a legible graph, but we’re able to use the output to infer
that an OU built from our simmap has similar parameter values to our
non-simmaped phylogeny.

``` r
OU_data.ordered<-data.frame(OU_data[,2], OU_data[,2],row.names=OU_data[,1])
OU_data.ordered<- OU_data.ordered[one_discrete_tree$tip.label,]
z<-OU_data.ordered[,1]
names(z)<-rownames(OU_data.ordered)
tree.mapped<-make.simmap(tree,z,model="ER",nsim=1)
```

    ## make.simmap is sampling character histories conditioned on
    ## the transition matrix
    ## 
    ## Q =
    ##            0          1
    ## 0 -0.2841324  0.2841324
    ## 1  0.2841324 -0.2841324
    ## (estimated using likelihood);
    ## and (mean) root node prior probabilities
    ## pi =
    ##   0   1 
    ## 0.5 0.5

    ## Done.

``` r
leg<-c("black","red")
names(leg)<-c(1,2)
plotSimmap(tree.mapped,leg,pts=FALSE,ftype="off", lwd=1)
```

![](Continuous_data_exercise_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
simmapBased<-OUwie(tree.mapped,OU_data,model="OUM", simmap.tree=TRUE, diagn=FALSE)
```

    ## Initializing... 
    ## Finished. Begin thorough search... 
    ## Finished. Summarizing results.

``` r
print(simmapBased)
```

    ## 
    ## Fit
    ##        lnL      AIC     AICc      BIC model ntax
    ##  -50.38181 108.7636 111.4303 112.7466   OUM   20
    ## 
    ## 
    ## Rates
    ##                  0         1
    ## alpha     2.562857  2.562857
    ## sigma.sq 51.765394 51.765394
    ## 
    ## Optima
    ##                  0         1
    ## estimate 12.395488 23.493026
    ## se        1.028955  3.780731
    ## 
    ## 
    ## Half life (another way of reporting alpha)
    ##         0         1 
    ## 0.2704587 0.2704587 
    ## 
    ## Arrived at a reliable solution

``` r
print(nodeBased.OUM)
```

    ## 
    ## Fit
    ##       lnL      AIC     AICc      BIC model ntax
    ##  -52.0688 112.1376 114.8043 116.1205   OUM   20
    ## 
    ## 
    ## Rates
    ##                  0         1
    ## alpha     1.941445  1.941445
    ## sigma.sq 52.397632 52.397632
    ## 
    ## Optima
    ##                  0        1
    ## estimate 11.824077 16.62596
    ## se        1.659426  2.91117
    ## 
    ## 
    ## Half life (another way of reporting alpha)
    ##         0         1 
    ## 0.3570263 0.3570263 
    ## 
    ## Arrived at a reliable solution
