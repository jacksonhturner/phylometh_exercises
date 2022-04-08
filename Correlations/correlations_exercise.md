Correlations Exercise
================
Jackson Turner
4/7/2022

## Initial

``` r
library(geiger)
```

    ## Loading required package: ape

``` r
library(corHMM)
```

    ## Loading required package: nloptr

    ## Loading required package: GenSA

``` r
library(ape)
tree.primates <- read.tree(text="((((Homo:0.21,Pongo:0.21):0.28,Macaca:0.49):0.13,Ateles:0.62):0.38,Galago:1.00);") #using examples from ape ?pic
X <- c(4.09434, 3.61092, 2.37024, 2.02815, -1.46968)
Y <- c(4.74493, 3.33220, 3.36730, 2.89037, 2.30259)
names(X) <- names(Y) <- c("Homo", "Pongo", "Macaca", "Ateles", "Galago")
pic.X <- pic(X, tree.primates)
pic.Y <- pic(Y, tree.primates)

require("corHMM")
data(primates)
ls()
```

    ## [1] "pic.X"         "pic.Y"         "primates"      "tree.primates"
    ## [5] "X"             "Y"

``` r
print(primates)
```

    ## $tree
    ## 
    ## Phylogenetic tree with 60 tips and 58 internal nodes.
    ## 
    ## Tip labels:
    ##   Homo_sapiens, Pan_paniscus, Pan_troglodytes, Gorilla_gorilla, Pongo_pygmaeus, Pongo_pygmaeus_abelii, ...
    ## 
    ## Rooted; includes branch lengths.
    ## 
    ## $trait
    ##                    Genus_sp T1 T2
    ## 1      Cercocebus_torquatus  1  1
    ## 2    Cercopithecus_aethiops  0  1
    ## 3        Cercopithecus_mona  0  0
    ## 4   Cercopithecus_nictitans  0  0
    ## 5        Colobus_angolensis  0  1
    ## 6           Colobus_guereza  0  0
    ## 7         Colobus_polykomos  0  1
    ## 8           Gorilla_gorilla  0  0
    ## 9              Homo_sapiens  0  0
    ## 10         Hylobates_agilis  0  0
    ## 11       Hylobates_concolor  0  0
    ## 12     Hylobates_gabriellae  0  0
    ## 13        Hylobates_hoolock  0  0
    ## 14        Hylobates_klossii  0  0
    ## 15            Hylobates_lar  0  0
    ## 16     Hylobates_leucogenys  0  0
    ## 17         Hylobates_moloch  0  0
    ## 18       Hylobates_muelleri  0  0
    ## 19       Hylobates_pileatus  0  0
    ## 20    Hylobates_syndactylus  0  0
    ## 21       Macaca_brunnescens  1  1
    ## 22          Macaca_cyclopis  0  1
    ## 23      Macaca_fascicularis  1  1
    ## 24             Macaca_hecki  1  1
    ## 25             Macaca_maura  1  1
    ## 26           Macaca_mulatta  1  1
    ## 27        Macaca_nemestrina  1  1
    ## 28             Macaca_nigra  1  1
    ## 29        Macaca_nigriscens  1  1
    ## 30          Macaca_ochreata  1  1
    ## 31           Macaca_silenus  1  1
    ## 32          Macaca_sylvanus  1  1
    ## 33          Macaca_tonkeana  1  1
    ## 34   Mandrillus_leucophaeus  1  1
    ## 35        Mandrillus_sphinx  1  1
    ## 36         Nasalis_larvatus  0  0
    ## 37             Pan_paniscus  1  1
    ## 38          Pan_troglodytes  1  1
    ## 39             Papio_anubis  1  1
    ## 40       Papio_cynocephalus  1  1
    ## 41          Papio_hamadryas  1  1
    ## 42           Pongo_pygmaeus  0  0
    ## 43    Pongo_pygmaeus_abelii  0  0
    ## 44      Presbytis_francoisi  0  0
    ## 45        Presbytis_phayrei  0  0
    ## 46          Presbytis_senex  0  0
    ## 47        Procolobus_badius  1  1
    ## 48        Pygathrix_nemaeus  0  1
    ## 49      Pygathrix_roxellana  0  1
    ## 50  Rhinopithecus_avunculus  0  1
    ## 51      Rhinopithecus_bieti  0  1
    ## 52   Semnopithecus_entellus  0  1
    ## 53 Trachypithecus_cristatus  0  0
    ## 54 Trachypithecus_francoisi  0  0
    ## 55      Trachypithecus_geei  0  0
    ## 56    Trachypithecus_johnii  0  0
    ## 57   Trachypithecus_phayrei  0  0
    ## 58  Trachypithecus_pileatus  0  0
    ## 59   Trachypithecus_vetulus  0  0
    ## 60         Macaca_arctoides  0  1

``` r
require(phytools)
```

    ## Loading required package: phytools

    ## Loading required package: maps

``` r
primates$trait[which(grepl("Hylobates",primates$trait[,1])),2]<-1

trait1<-primates$trait[,2]
names(trait1)<-primates$trait[,1]
primates$tree <- ape::multi2di(primates$tree)
plotSimmap(make.simmap(primates$tree, trait1), pts=FALSE, fsize=0.8)
```

    ## make.simmap is sampling character histories conditioned on
    ## the transition matrix
    ## 
    ## Q =
    ##             0           1
    ## 0 -0.01076402  0.01076402
    ## 1  0.01076402 -0.01076402
    ## (estimated using likelihood);
    ## and (mean) root node prior probabilities
    ## pi =
    ##   0   1 
    ## 0.5 0.5

    ## Done.

    ## no colors provided. using the following legend:
    ##         0         1 
    ##   "black" "#DF536B"

![](correlations_exercise_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
rate.mat.er<-corHMM:::rate.mat.maker(rate.cat=1, hrm=FALSE, ntraits=1, nstates=2, model="ER")
print(rate.mat.er)
```

    ##    1  2
    ## 1 NA  1
    ## 2  1 NA

``` r
# Same rates

pp.er<-corHMM(primates$tree,primates$trait[,c(1,2)],rate.cat=1,rate.mat=rate.mat.er,node.states="marginal")
```

    ## State distribution in data:
    ## States:  1   2   
    ## Counts:  28  32  
    ## Beginning thorough optimization search -- performing 0 random restarts 
    ## Finished. Inferring ancestral states using marginal reconstruction.

``` r
print(pp.er)
```

    ## 
    ## Fit
    ##       -lnL     AIC     AICc Rate.cat ntax
    ##  -23.41535 48.8307 48.89967        1   60
    ## 
    ## Legend
    ##   1   2 
    ## "0" "1" 
    ## 
    ## Rates
    ##            (1,R1)     (2,R1)
    ## (1,R1)         NA 0.01076323
    ## (2,R1) 0.01076323         NA
    ## 
    ## Arrived at a reliable solution

``` r
# Forward and backward rates are the same; 1 and 2 will be about as common after lots of time.

rate.mat.ard<-corHMM:::rate.mat.maker(rate.cat=1, hrm=FALSE, ntraits=1, nstates=2, model="ARD")
print(rate.mat.ard)
```

    ##    1  2
    ## 1 NA  2
    ## 2  1 NA

``` r
# Rates are different.

pp.ard<-corHMM(primates$tree,primates$trait[,c(1,2)],rate.cat=1,rate.mat=rate.mat.ard,node.states="marginal")
```

    ## State distribution in data:
    ## States:  1   2   
    ## Counts:  28  32  
    ## Beginning thorough optimization search -- performing 0 random restarts 
    ## Finished. Inferring ancestral states using marginal reconstruction.

``` r
print(pp.ard)
```

    ## 
    ## Fit
    ##      -lnL     AIC     AICc Rate.cat ntax
    ##  -23.4031 50.8062 51.01673        1   60
    ## 
    ## Legend
    ##   1   2 
    ## "0" "1" 
    ## 
    ## Rates
    ##            (1,R1)     (2,R1)
    ## (1,R1)         NA 0.01012564
    ## (2,R1) 0.01144043         NA
    ## 
    ## Arrived at a reliable solution

``` r
# Backward (2 to 1) has higher rate than forward; more 2 will appear after lots of time.

# Which of these three is better?
# pp.ard is probably better than pp.er because the rates are probably at least slightly different from each other in nature.

rate.mat.er.4state<-corHMM:::rate.mat.maker(rate.cat=1, hrm=FALSE, ntraits=1, nstates=4, model="ER")
print(rate.mat.er.4state)
```

    ##    1  2  3  4
    ## 1 NA  1  1  1
    ## 2  1 NA  1  1
    ## 3  1  1 NA  1
    ## 4  1  1  1 NA

``` r
rate.mat.ard.4state<-corHMM:::rate.mat.maker(rate.cat=1, hrm=FALSE, ntraits=1, nstates=4, model="ARD")
print(rate.mat.ard.4state)
```

    ##    1  2  3  4
    ## 1 NA  4  7 10
    ## 2  1 NA  8 11
    ## 3  2  5 NA 12
    ## 4  3  6  9 NA

``` r
rate.mat.gtr.4state<-rate.mat.ard.4state
rate.mat.gtr.4state<-corHMM:::rate.par.eq(rate.mat.gtr.4state, c(1,4))
rate.mat.gtr.4state<-corHMM:::rate.par.eq(rate.mat.gtr.4state, c(2,6))
rate.mat.gtr.4state<-corHMM:::rate.par.eq(rate.mat.gtr.4state, c(3,8))
rate.mat.gtr.4state<-corHMM:::rate.par.eq(rate.mat.gtr.4state, c(4,6))
rate.mat.gtr.4state<-corHMM:::rate.par.eq(rate.mat.gtr.4state, c(5,7))
rate.mat.gtr.4state<-corHMM:::rate.par.eq(rate.mat.gtr.4state, c(6,7))
print(rate.mat.gtr.4state)
```

    ##    1  2  3  4
    ## 1 NA  1  2  3
    ## 2  1 NA  4  5
    ## 3  2  4 NA  6
    ## 4  3  5  6 NA

``` r
fourstate.trait<-rep(NA,Ntip(primates$tree))
for(i in sequence(Ntip(primates$tree))) {
  if(primates$trait[i,2]==0 && primates$trait[i,3]==0) {
    fourstate.trait[i]<-0
  }
  if(primates$trait[i,2]==0 && primates$trait[i,3]==1) {
    fourstate.trait[i]<-1
  }
  if(primates$trait[i,2]==1 && primates$trait[i,3]==0) {
    fourstate.trait[i]<-2
  }
  if(primates$trait[i,2]==1 && primates$trait[i,3]==1) {
    fourstate.trait[i]<-3
  }
}
fourstate.data<-data.frame(Genus_sp=primates$trait[,1], T1=fourstate.trait)

print(rayDISC(primates$tree, fourstate.data, ntraits=1, rate.mat= rate.mat.gtr.4state, node.states="marginal", model="ARD"))
```

    ## State distribution in data:
    ## States:  0   1   2   3   
    ## Counts:  18  10  11  21  
    ## Initializing... 
    ## Finished. Beginning thorough search... 
    ## Finished. Inferring ancestral states using marginal reconstruction. 
    ## 
    ## Fit
    ##       -lnL      AIC    AICc ntax
    ##  -44.78056 101.5611 103.146   60
    ## 
    ## Rates
    ##             0            1            2            3
    ## 0          NA 1.000000e+02 2.906960e-03 1.312130e-02
    ## 1 1.00000e+02           NA 4.286561e-06 4.112289e-03
    ## 2 2.90696e-03 4.286561e-06           NA 7.582560e-10
    ## 3 1.31213e-02 4.112289e-03 7.582560e-10           NA
    ## 
    ## Arrived at a reliable solution

``` r
print(rayDISC(primates$tree, fourstate.data, ntraits=1, model="ER", node.states="marginal"))
```

    ## State distribution in data:
    ## States:  0   1   2   3   
    ## Counts:  18  10  11  21  
    ## Initializing... 
    ## Finished. Beginning thorough search... 
    ## Finished. Inferring ancestral states using marginal reconstruction. 
    ## 
    ## Fit
    ##       -lnL      AIC     AICc ntax
    ##  -52.96386 107.9277 107.9967   60
    ## 
    ## Rates
    ##             0           1           2           3
    ## 0          NA 0.006599861 0.006599861 0.006599861
    ## 1 0.006599861          NA 0.006599861 0.006599861
    ## 2 0.006599861 0.006599861          NA 0.006599861
    ## 3 0.006599861 0.006599861 0.006599861          NA
    ## 
    ## Arrived at a reliable solution

``` r
print(rayDISC(primates$tree, fourstate.data, ntraits=1, rate.mat=rate.mat.er.4state, node.states="marginal", model="ARD"))
```

    ## State distribution in data:
    ## States:  0   1   2   3   
    ## Counts:  18  10  11  21  
    ## Initializing... 
    ## Finished. Beginning thorough search... 
    ## Finished. Inferring ancestral states using marginal reconstruction. 
    ## 
    ## Fit
    ##       -lnL      AIC     AICc ntax
    ##  -52.96386 107.9277 107.9967   60
    ## 
    ## Rates
    ##             0           1           2           3
    ## 0          NA 0.006598762 0.006598762 0.006598762
    ## 1 0.006598762          NA 0.006598762 0.006598762
    ## 2 0.006598762 0.006598762          NA 0.006598762
    ## 3 0.006598762 0.006598762 0.006598762          NA
    ## 
    ## Arrived at a reliable solution

``` r
rate.mat.ard.4state<-corHMM:::rate.mat.maker(rate.cat=1, hrm=FALSE, ntraits=1, nstates=4, model="ARD")
print(rate.mat.ard.4state)
```

    ##    1  2  3  4
    ## 1 NA  4  7 10
    ## 2  1 NA  8 11
    ## 3  2  5 NA 12
    ## 4  3  6  9 NA

``` r
print(corHMM:::rate.mat.maker(rate.cat=1, hrm=FALSE, ntraits=2, nstates=2, model="ARD"))
```

    ##       (0,0) (0,1) (1,0) (1,1)
    ## (0,0)    NA     3     5    NA
    ## (0,1)     1    NA    NA     7
    ## (1,0)     2    NA    NA     8
    ## (1,1)    NA     4     6    NA

``` r
rate.mat.pag94<-corHMM:::rate.par.drop(rate.mat.ard.4state, drop.par=c(3,5,8,10))
print(rate.mat.pag94)
```

    ##    1  2  3  4
    ## 1 NA  3  5 NA
    ## 2  1 NA NA  7
    ## 3  2 NA NA  8
    ## 4 NA  4  6 NA

## Route 1

``` r
pagelsim<-function(v00,v10,v01,v11,silence_output=FALSE){
 
  r00t10<-6.9
  r00t01<-0.00121
  r10t11<-0.0172
  r10t00<-6.9
  r01t11<-0.00000000076
  r01t00<-0.0121
  r11t10<-0.017
  r11t01<-0.00000000076
  r_sum<-sum(c(r00t10,r00t01,r10t11,r10t00,r01t11,r01t00,r11t10,r11t01))
  sr00t10<-r00t10/r_sum
  sr00t01<-r00t01/r_sum
  sr10t11<-r10t11/r_sum
  sr10t00<-r10t00/r_sum
  sr01t11<-r01t11/r_sum
  sr01t00<-r01t00/r_sum
  sr11t10<-r11t10/r_sum
  sr11t01<-r11t01/r_sum
  sr_sum<-c(sr00t10,sr00t01,sr10t11,sr10t00,sr01t11,sr01t00,sr11t10,sr11t01)
   
steps_val<-0
ding<-FALSE

while((ding==FALSE)&(v00>0)&(v10>0)&(v01>0)&(v11>0)){
 if(silence_output==TRUE){
  r1<-runif(1,0,1)
 
  if((r1>0)&(r1<sr_sum[1])){
    v00<-v00-1
    v10<-v10+1
  }
  if((r1>sum(sr_sum[1]))&(r1<sum(sr_sum[1:2]))){
    v00<-v00-1
    v01<-v01+1
  }
  if((r1>sum(sr_sum[1:2]))&(r1<sum(sr_sum[1:3]))){
    v10<-v10-1
    v11<-v11+1
  }
  if((r1>sum(sr_sum[1:3]))&(r1<sum(sr_sum[1:4]))){
    v10<-v10-1
    v00<-v00+1
  }
  if((r1>sum(sr_sum[1:4]))&(r1<sum(sr_sum[1:5]))){
    v01<-v01-1
    v11<-v11+1
  }
  if((r1>sum(sr_sum[1:5]))&(r1<sum(sr_sum[1:6]))){
    v01<-v01-1
    v11<-v11+1
    ding<-TRUE
  }
  if((r1>sum(sr_sum[1:6]))&(r1<sum(sr_sum[1:7]))){
    v11<-v11-1
    v10<-v10+1
  }
  if((r1>sum(sr_sum[1:7]))&(r1<sum(sr_sum[1:8]))){
    v11<-v11-1
    v01<-v01+1
  }
  steps_val<-steps_val+1
  states_end<<-c(v00,v10,v01,v11)
 }else{
  r1<-runif(1,0,1)
 
  if((r1>0)&(r1<sr_sum[1])){
    v00<-v00-1
    v10<-v10+1
  }
  if((r1>sum(sr_sum[1]))&(r1<sum(sr_sum[1:2]))){
    v00<-v00-1
    v01<-v01+1
  }
  if((r1>sum(sr_sum[1:2]))&(r1<sum(sr_sum[1:3]))){
    v10<-v10-1
    v11<-v11+1
  }
  if((r1>sum(sr_sum[1:3]))&(r1<sum(sr_sum[1:4]))){
    v10<-v10-1
    v00<-v00+1
  }
  if((r1>sum(sr_sum[1:4]))&(r1<sum(sr_sum[1:5]))){
    v01<-v01-1
    v11<-v11+1
  }
  if((r1>sum(sr_sum[1:5]))&(r1<sum(sr_sum[1:6]))){
    v01<-v01-1
    v11<-v11+1
    ding<-TRUE
  }
  if((r1>sum(sr_sum[1:6]))&(r1<sum(sr_sum[1:7]))){
    v11<-v11-1
    v10<-v10+1
  }
  if((r1>sum(sr_sum[1:7]))&(r1<sum(sr_sum[1:8]))){
    v11<-v11-1
    v01<-v01+1
  }
  if((ding==TRUE)){
    cat("Allele 01 has transferred to allele 11 after",steps_val,"time steps. The numbers of each allele are as follows: 00 =",v00,", 10 =",v10,", 01 =",v01,", 11 =",v11)
  }
  if((v00==0)|(v10==0)|(v01==0)|(v11==0)){
    cat("An allele has become extinct from this population after",steps_val,"time steps. The numbers of each allele are as follows: 00 =",v00,", 10 =",v10,", 01 =",v01,", 11 =",v11)
  }
  steps_val<-steps_val+1
  states_end<<-c(v00,v10,v01,v11)
 }}
}

pagelsim(100,100,100,100) # try it out! change the number of starting alleles and the rates within the function. 
```

    ## Allele 01 has transferred to allele 11 after 1148 time steps. The numbers of each allele are as follows: 00 = 136 , 10 = 66 , 01 = 99 , 11 = 99

``` r
# It tells you how many time steps it takes for 01 to transfer to 11 unless one of the alleles becomes extinct first.

# Is 00 ever lost? This should tell us hopefully:

n_attempts<-0
states_end<-c(100,100,100,100)
n_final<-as.vector(NULL)

for(i in 1:100){
  n_attempts<-0
  states_end<-c(100,100,100,100)
while(states_end[1]>0){
  pagelsim(100,100,100,100,TRUE)
  n_attempts<-n_attempts+1
  if(states_end[1]==0){
    cat("Run",i,": Allele 00 has become extinct after",n_attempts,"simulation runs.")
  }
}
  n_final[i]<-n_attempts
  cat("\n")
}
```

    ## Run 1 : Allele 00 has become extinct after 5 simulation runs.
    ## Run 2 : Allele 00 has become extinct after 247 simulation runs.
    ## Run 3 : Allele 00 has become extinct after 14 simulation runs.
    ## Run 4 : Allele 00 has become extinct after 6 simulation runs.
    ## Run 5 : Allele 00 has become extinct after 22 simulation runs.
    ## Run 6 : Allele 00 has become extinct after 12 simulation runs.
    ## Run 7 : Allele 00 has become extinct after 15 simulation runs.
    ## Run 8 : Allele 00 has become extinct after 68 simulation runs.
    ## Run 9 : Allele 00 has become extinct after 6 simulation runs.
    ## Run 10 : Allele 00 has become extinct after 182 simulation runs.
    ## Run 11 : Allele 00 has become extinct after 82 simulation runs.
    ## Run 12 : Allele 00 has become extinct after 11 simulation runs.
    ## Run 13 : Allele 00 has become extinct after 140 simulation runs.
    ## Run 14 : Allele 00 has become extinct after 32 simulation runs.
    ## Run 15 : Allele 00 has become extinct after 13 simulation runs.
    ## Run 16 : Allele 00 has become extinct after 46 simulation runs.
    ## Run 17 : Allele 00 has become extinct after 114 simulation runs.
    ## Run 18 : Allele 00 has become extinct after 7 simulation runs.
    ## Run 19 : Allele 00 has become extinct after 5 simulation runs.
    ## Run 20 : Allele 00 has become extinct after 20 simulation runs.
    ## Run 21 : Allele 00 has become extinct after 225 simulation runs.
    ## Run 22 : Allele 00 has become extinct after 3 simulation runs.
    ## Run 23 : Allele 00 has become extinct after 15 simulation runs.
    ## Run 24 : Allele 00 has become extinct after 39 simulation runs.
    ## Run 25 : Allele 00 has become extinct after 2 simulation runs.
    ## Run 26 : Allele 00 has become extinct after 10 simulation runs.
    ## Run 27 : Allele 00 has become extinct after 46 simulation runs.
    ## Run 28 : Allele 00 has become extinct after 240 simulation runs.
    ## Run 29 : Allele 00 has become extinct after 15 simulation runs.
    ## Run 30 : Allele 00 has become extinct after 23 simulation runs.
    ## Run 31 : Allele 00 has become extinct after 144 simulation runs.
    ## Run 32 : Allele 00 has become extinct after 64 simulation runs.
    ## Run 33 : Allele 00 has become extinct after 17 simulation runs.
    ## Run 34 : Allele 00 has become extinct after 7 simulation runs.
    ## Run 35 : Allele 00 has become extinct after 60 simulation runs.
    ## Run 36 : Allele 00 has become extinct after 114 simulation runs.
    ## Run 37 : Allele 00 has become extinct after 30 simulation runs.
    ## Run 38 : Allele 00 has become extinct after 106 simulation runs.
    ## Run 39 : Allele 00 has become extinct after 17 simulation runs.
    ## Run 40 : Allele 00 has become extinct after 10 simulation runs.
    ## Run 41 : Allele 00 has become extinct after 41 simulation runs.
    ## Run 42 : Allele 00 has become extinct after 52 simulation runs.
    ## Run 43 : Allele 00 has become extinct after 62 simulation runs.
    ## Run 44 : Allele 00 has become extinct after 67 simulation runs.
    ## Run 45 : Allele 00 has become extinct after 76 simulation runs.
    ## Run 46 : Allele 00 has become extinct after 71 simulation runs.
    ## Run 47 : Allele 00 has become extinct after 14 simulation runs.
    ## Run 48 : Allele 00 has become extinct after 93 simulation runs.
    ## Run 49 : Allele 00 has become extinct after 219 simulation runs.
    ## Run 50 : Allele 00 has become extinct after 19 simulation runs.
    ## Run 51 : Allele 00 has become extinct after 40 simulation runs.
    ## Run 52 : Allele 00 has become extinct after 78 simulation runs.
    ## Run 53 : Allele 00 has become extinct after 94 simulation runs.
    ## Run 54 : Allele 00 has become extinct after 16 simulation runs.
    ## Run 55 : Allele 00 has become extinct after 42 simulation runs.
    ## Run 56 : Allele 00 has become extinct after 29 simulation runs.
    ## Run 57 : Allele 00 has become extinct after 1 simulation runs.
    ## Run 58 : Allele 00 has become extinct after 45 simulation runs.
    ## Run 59 : Allele 00 has become extinct after 79 simulation runs.
    ## Run 60 : Allele 00 has become extinct after 42 simulation runs.
    ## Run 61 : Allele 00 has become extinct after 14 simulation runs.
    ## Run 62 : Allele 00 has become extinct after 47 simulation runs.
    ## Run 63 : Allele 00 has become extinct after 21 simulation runs.
    ## Run 64 : Allele 00 has become extinct after 40 simulation runs.
    ## Run 65 : Allele 00 has become extinct after 37 simulation runs.
    ## Run 66 : Allele 00 has become extinct after 14 simulation runs.
    ## Run 67 : Allele 00 has become extinct after 7 simulation runs.
    ## Run 68 : Allele 00 has become extinct after 65 simulation runs.
    ## Run 69 : Allele 00 has become extinct after 63 simulation runs.
    ## Run 70 : Allele 00 has become extinct after 13 simulation runs.
    ## Run 71 : Allele 00 has become extinct after 59 simulation runs.
    ## Run 72 : Allele 00 has become extinct after 133 simulation runs.
    ## Run 73 : Allele 00 has become extinct after 59 simulation runs.
    ## Run 74 : Allele 00 has become extinct after 42 simulation runs.
    ## Run 75 : Allele 00 has become extinct after 38 simulation runs.
    ## Run 76 : Allele 00 has become extinct after 12 simulation runs.
    ## Run 77 : Allele 00 has become extinct after 9 simulation runs.
    ## Run 78 : Allele 00 has become extinct after 21 simulation runs.
    ## Run 79 : Allele 00 has become extinct after 142 simulation runs.
    ## Run 80 : Allele 00 has become extinct after 74 simulation runs.
    ## Run 81 : Allele 00 has become extinct after 126 simulation runs.
    ## Run 82 : Allele 00 has become extinct after 73 simulation runs.
    ## Run 83 : Allele 00 has become extinct after 71 simulation runs.
    ## Run 84 : Allele 00 has become extinct after 27 simulation runs.
    ## Run 85 : Allele 00 has become extinct after 71 simulation runs.
    ## Run 86 : Allele 00 has become extinct after 36 simulation runs.
    ## Run 87 : Allele 00 has become extinct after 41 simulation runs.
    ## Run 88 : Allele 00 has become extinct after 35 simulation runs.
    ## Run 89 : Allele 00 has become extinct after 77 simulation runs.
    ## Run 90 : Allele 00 has become extinct after 76 simulation runs.
    ## Run 91 : Allele 00 has become extinct after 238 simulation runs.
    ## Run 92 : Allele 00 has become extinct after 40 simulation runs.
    ## Run 93 : Allele 00 has become extinct after 29 simulation runs.
    ## Run 94 : Allele 00 has become extinct after 24 simulation runs.
    ## Run 95 : Allele 00 has become extinct after 38 simulation runs.
    ## Run 96 : Allele 00 has become extinct after 55 simulation runs.
    ## Run 97 : Allele 00 has become extinct after 8 simulation runs.
    ## Run 98 : Allele 00 has become extinct after 2 simulation runs.
    ## Run 99 : Allele 00 has become extinct after 51 simulation runs.
    ## Run 100 : Allele 00 has become extinct after 30 simulation runs.

``` r
hist(n_final,xlab="Number of pagelsim runs before 00 extinction",ylab="Frequency",main="Frequency of the number of pagelsim runs before 00 extinction") # Yes, and not infrequently
```

![](correlations_exercise_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->
