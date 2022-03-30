library(ape)
library(geiger)
library(phytools)
library(OUwie)
library(akima)

#geo<-get(data(geospiza))

tree <- read.tree("RAxML_nodeBased.OUMTree.T2.tre")
continuous_data<-read.csv("Co1_continuous.csv",stringsAsFactors=FALSE)
rownames(continuous_data)<-continuous_data[,1]
continuous_data<-continuous_data[,-1]
tree<-multi2di(tree)
tree<-chronos(tree) #to make ultrametric

#tree<-geo$phy
#continuous_data<-geo$dat

cleaned_continuous<-treedata(tree,continuous_data,sort=TRUE)
data_vec <- cleaned_continuous$data[,1]
names(data_vec) <- rownames(cleaned_continuous$data)
contMap(tree,data_vec)

BM1<-fitContinuous(tree, cleaned_continuous$data, model="BM")
cat("The rate of evolution is", BM1$tailpipe_length_width_th$opt$sigsq, "in units of","um^2 per length of time.","\n")
      
OU1 <- fitContinuous(tree, cleaned_continuous$data, model="OU")
par(mfrow=(c(1,2)))
plot(tree, show.tip.label=FALSE)
ou.tree <- rescale(tree, model="OU",OU1$tailpipe_length_width$opt$alpha)
plot(ou.tree)

cat("The AIC score for BM1 is", BM1$tailpipe_length_width$opt$aic, "\n")
cat("The AIC score for OU1 is", OU1$tailpipe_length_width$opt$aic, "\n")
cat("The delta AIC for BM1 is", BM1$tailpipe_length_width$opt$aic-OU1$tailpipe_length_width$opt$aic, "\n")
cat("The delta AIC for OU1 is", OU1$tailpipe_length_width$opt$aic-OU1$tailpipe_length_width$opt$aic, "\n")

#wow_factor<-(c("nice","nice","rad","rad","rad","rad","nice","nice","rad","nice","rad","nice","nice","rad"))
#row.names(wow_factor)<-row.names(continuous_data)
#colnames(wow_factor)<-"how_wow"
#cleaned_discrete<-treedata(tree,wow_factor,sort=TRUE)

discrete_data<-read.csv("Co1_discrete.csv")
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

#tree$tip.label<-tailpipe_fused
#making internal node labels for OUwie
#tree<-makeNodeLabel(tree,"n",nodeBased.OUM.states)

#OUframe<-as.data.frame(rownames(cleaned_continuous$data))
#OUframe$regime<-one.discrete.char
#OUframe$cont_OU_data<-rep("tailpipe_length_width",20)

one_discrete_tree<-tree
one_discrete_tree$node.label<-vector()
for (i in 1:one_discrete_tree$Nnode){ 
  one_discrete_tree$node.label[i]<-nodeBased.OUM.states[i]
}
one_discrete_tree$node.label
one_discrete_tree<-chronos(one_discrete_tree)

nodecolr<-rep("green", length=one_discrete_tree$Nnode)
nodecolr[one_discrete_tree$node.label==1]<-"purple"
nodecolr[one_discrete_tree$node.label==0]<-"green"
plot(one_discrete_tree)
nodelabels(pch=16, col=nodecolr)

continuous_data_o<-read.csv("Co1_continuous.csv",stringsAsFactors=FALSE)
check.identify(one_discrete_tree,continuous_data_o)

OU_data<-as.data.frame(continuous_data_o[,1])
OU_data$regime<-cleaned_discrete$data[,2]
OU_data$tailpipe_length_width<-cleaned_continuous$data[,1]
colnames(OU_data)<-c("Species","Reg","OU_data")

nodeBased.OUMV <- OUwie(one_discrete_tree, OU_data,model="OUMV", simmap.tree=FALSE,diagn=FALSE)
print(nodeBased.OUMV)

nodeBased.BM1 <- OUwie(one_discrete_tree, OU_data,model="BM1", simmap.tree=FALSE, diagn=FALSE)
nodeBased.BMS <- OUwie(one_discrete_tree, OU_data,model="BMS", simmap.tree=FALSE, diagn=FALSE)
nodeBased.OUM <- OUwie(one_discrete_tree, OU_data,model="OUM", simmap.tree=FALSE, diagn=FALSE)
nodeBased.OUMA <- OUwie(one_discrete_tree, OU_data,model="OUMA", simmap.tree=FALSE, diagn=FALSE)
nodeBased.OUMVA <- OUwie(one_discrete_tree, OU_data,model="OUMVA", simmap.tree=FALSE, diagn=FALSE)

modelAIC<-data.frame(cbind(nodeBased.BM1$AICc,nodeBased.BMS$AICc,nodeBased.OUM$AICc,nodeBased.OUMA$AICc,nodeBased.OUMV$AICc,nodeBased.OUMVA$AICc))
colnames(modelAIC)<-c("BM1","BMS","OUM","OUMA","OUMA","OUMVA")
rownames(modelAIC)<-"corrected AIC score:"
modelAIC #BM1 is nodeBased.OUM but using OUM since it's close and has an alpha parameter

nodeBased.OUM
OUwie.fixed(one_discrete_tree,OU_data,model="OUM",alpha=nodeBased.OUM$solution[1,],sigma.sq=nodeBased.OUM$solution[2,],theta=nodeBased.OUM$theta[,1])

alpha.values<-seq(from=0.5,to=3,length.out=10)
likelihood.values <- rep(NA, length(alpha.values))

for (i in sequence(length(alpha.values))) {
  likelihood.values[i] <- OUwie.fixed(one_discrete_tree, OU_data, model="OUM", alpha=rep(alpha.values[i],2), sigma.sq=nodeBased.OUM$solution[2,], theta=nodeBased.OUM$theta[,1])$loglik
}

plot(x=alpha.values, y=likelihood.values, xlab="Alpha", ylab="Likelihood", type="l", bty="n")
points(x=nodeBased.OUM$solution[1,1], y=nodeBased.OUM$loglik, pch=16, col="red")
text(x=nodeBased.OUM$solution[1,1], y=nodeBased.OUM$loglik, "unconstrained OUM model", pos=4, col="red")
abline(h=(which(likelihood.values)==max)-2, lty="dotted") #FIX ME

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

contour(interpolated.points, xlim=range(c(theta1.points, theta2.points)),ylim=range(c(theta1.points, theta2.points)), xlab="Theta 1", ylab="Theta 2", levels=c(2,5,10),add=FALSE,lwd=1, bty="n", asp=1)

points(x=nodeBased.OUM$theta[1,1], y=nodeBased.OUM$theta[2,1], col="red", pch=16)

points(x=OU_data$OU_data[which(OU_data$Reg==1)],y=rep(min(c(theta1.points, theta2.points)), length(which(OU_data$Reg==1))), pch=18, col=rgb(0,0,0,.3)) #the tip values in regime 1, plotted along x axis
points(y=OU_data$OU_data[which(OU_data$Reg==2)],x=rep(min(c(theta1.points, theta2.points)), length(which(OU_data$Reg==2))), pch=18, col=rgb(0,0,0,.3)) #the tip values in regime 2, plotted along y axis

OU_data.ordered<-data.frame(OU_data[,2], OU_data[,2],row.names=OU_data[,1])
OU_data.ordered<- OU_data.ordered[one_discrete_tree$tip.label,]
z<-OU_data.ordered[,1]
names(z)<-rownames(OU_data.ordered)
tree.mapped<-make.simmap(tree,z,model="ER",nsim=1)
leg<-c("black","red")
names(leg)<-c(1,2)
plotSimmap(tree.mapped,leg,pts=FALSE,ftype="off", lwd=1)

simmapBased<-OUwie(tree.mapped,OU_data,model="OUMV", simmap.tree=TRUE, diagn=FALSE)
print(simmapBased)
print(nodeBased.OUM)
