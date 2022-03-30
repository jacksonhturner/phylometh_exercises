EEB 587 Week 7 Exercise
================
Jackson Turner
3/11/2022

More verbose version coming soon.

Load our packages.

``` r
library(ape)
library(geiger)
library(phangorn)
library(corHMM)
```

    ## Loading required package: nloptr

    ## Loading required package: GenSA

Load and format our data.

``` r
tree <- read.tree("RAxML_bestTree.T2.tre")
plot.phylo(tree) #plot the tree if you want
```

![](Ancestral_Reconstruction_Exercise_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
discrete.data <- (read.csv(file="Co1_states.csv", stringsAsFactors=FALSE)) #read the csv file

#all states are 0 = absent & 1 = present
#barbs = if tailpipes have barbs 
#meringo = if dorsal coronae protrude to sides, paddles long and skinny
#rgroup = if heads and postspiracular process are black and larva is bulky
#tailpipe_fused = if tailpipe is connected to the postspiracular sclerite
#unknown_to_me = if i've never seen or worked with the larva in person

rownames(discrete.data) <- discrete.data[,1] #formatting the data
discrete.data<-discrete.data[,-1]
```

Visualize our data to validate that it’s not weird and format it some
more for the ancestral reconstruction models.

``` r
cleaned.discrete<-treedata(tree,data=discrete.data,sort=TRUE)# using the treedata() function to quality check it
(visualize_data<-cleaned.discrete) #look at it; it looks fine
```

    ## $phy
    ## 
    ## Phylogenetic tree with 20 tips and 18 internal nodes.
    ## 
    ## Tip labels:
    ##   pycta, johannseni, borkenti, distincta, gelhausi, rudis, ...
    ## 
    ## Rooted; includes branch lengths.
    ## 
    ## $data
    ##                  barbs meringo rgroup tailpipe_fused unknown_to_me
    ## pycta                0       0      1              1             0
    ## johannseni           0       0      0              0             1
    ## borkenti             0       0      1              1             0
    ## distincta            0       0      1              1             0
    ## gelhausi             0       0      1              1             0
    ## rudis                0       0      1              1             0
    ## turneri              0       0      1              1             0
    ## sierranevadensis     0       0      0              0             1
    ## torrenticola         0       0      1              1             0
    ## auricula             0       0      0              1             0
    ## merjalisco           0       0      0              0             1
    ## lobata               0       1      0              1             0
    ## xavia                0       1      0              1             0
    ## calciphala           1       0      0              0             0
    ## inextricata          1       0      0              0             0
    ## appalachiensis       1       0      0              0             0
    ## rhathyme             0       0      0              1             0
    ## arge                 0       1      0              1             0
    ## midahoensis          0       0      0              1             0
    ## fluvica              0       0      0              0             0

``` r
discrete.phyDat <- phangorn::phyDat((cleaned.discrete$data), type="USER",levels=c(0,1)) #formatting for functions later
```

Run our ancestral reconstruction model using parsimony.

``` r
anc.p <- phangorn::ancestral.pars(tree, discrete.phyDat[,1]) #using parsimony for ancestral reconstruction
plotAnc(tree, anc.p, 1) #behold
```

![](Ancestral_Reconstruction_Exercise_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

There’s no uncertainty, which is the dissimilarity between produced
trees. This is because parsimony produces reconstructions that require
the same amount of changes no matter how many trees you make (are
equally parsimonious) (Losos, 1999).

``` r
anc.ml <- ancestral.pml(pml(tree, discrete.phyDat), type="ml") #using maximum likelihood for ancestral reconstruction
plotAnc(tree, anc.ml, 1) #pretty much the same but with more uncertainty
```

![](Ancestral_Reconstruction_Exercise_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

There is uncertainty in the maximum likelihood model since doesn’t
minimize the number of changes, but instead uses the data to determine
the most likely tree, leading to some constituent trees being different
from others.

Markov models are used to estimate transition rates between two discrete
character states. Here I estimate the transition rates between the
states for the “no barbs” trait using equation 4 in Tuffley and Steel
(1994).

``` r
k<-2 #number of states
a1<-0.1 #instantaneous rate of change
t1<-1 #amount in units time represented per branch length

rate_of_0_to_0<-(1/k)+(((k-1)/(k))*exp(-k*a1*t1)) #rate of no barbs to no barbs
cat("The transition rate from no barbs to no barbs is",rate_of_0_to_0)
```

    ## The transition rate from no barbs to no barbs is 0.9093654

``` r
rate_of_0_to_1<-(1/k)-(((k-1)/(k))*exp(-k*a1*t1)) #rate of no barbs to barbs
cat("The transition rate from no barbs to barbs is",rate_of_0_to_1)
```

    ## The transition rate from no barbs to barbs is 0.09063462

Now for the “barbs” trait.

``` r
a2<-0.2 #instantaneous rate of change
t2<-1 #amount in units time represented per branch length

rate_of_1_to_1<-(1/k)+(((k-1)/(k))*exp(-k*a2*t2)) #rate of barbs to barbs
cat("The transition rate from no barbs to barbs is",rate_of_1_to_1)
```

    ## The transition rate from no barbs to barbs is 0.83516

``` r
rate_of_1_to_0<-(1/k)-(((k-1)/(k))*exp(-k*a2*t2)) #rate of barbs to no barbs
cat("The transition rate from barbs to no barbs is",rate_of_1_to_0)
```

    ## The transition rate from barbs to no barbs is 0.16484

If the parameters estimating the instantaneous rate of change and the
amount in unit time represented per branch length were the equal, then
so would the transition rate probabilities.

My traits are all variable as defined in Lewis (2001), so the MKV model
would make sense to use given my data. I attempt to use it below:

``` r
DE_given_TA<-0.8 #the probability of no barbs (assumed to be variable) given the time represented per branch length and our rate of change
E<-0.9

probability_no_barbs_given_para<-DE_given_TA/E
cat("The probability of no barbs (variable) given the MkV model is",probability_no_barbs_given_para)
```

    ## The probability of no barbs (variable) given the MkV model is 0.8888889

``` r
DE_given_TA2<-0.05#the probability of barbs (assumed to be variable) given the time represented per branch length and our rate of change
E2<-0.9

probability_barbs_given_para<-DE_given_TA2/E2
cat("The probability of barbs (variable) given the MkV model is",probability_barbs_given_para)
```

    ## The probability of barbs (variable) given the MkV model is 0.05555556

The divergence order test (DOT), which is based on the average age of
the nodes on a tree, weighted by the absolute magnitude of the contrast
at each node for a particular trait, could be used to test the evolution
of our states (Ackerly et al., (2006)).

Last updated 3/12/22.
