# Hamadryas 

library("geomorph") 
library("phytools") 
library("geiger") 
library("MASS")
library("plyr")

setwd("~/Desktop/Projects/Hamadryas/Hama_data")

(SessID <- sessionInfo())
# save(list = ls(), file = "Hama_data_1.R")
load("Hama_data_1.R")
# packageDescription("ape")$Version


## Read in the tree and make ultrametric
Hama <- read.nexus("../Tree_data/Hamydryas_ml.tre")
plot(Hama)
is.ultrametric(Hama)

# make tree ultrametric, lamdba is the smoothing parameter. 
lam <- 10^(-1:6)
cv <- sapply(lam, function(x) sum(attr(chronopl(Hama, lambda = x, CV = TRUE), "D2")))
plot(lam, cv, pch = 19, ylab = "cross-validation score", xlab = expression(paste(lambda)), las = 1, cex = 1.5) # lowest CV 

Hama2 <- chronopl(phy = Hama, lambda = 0.1, CV = TRUE, eval.max = 1e3, iter.max = 1e4)
is.ultrametric(Hama2)

plot(Hama2)

plot(attr(Hama2, "D2"), type = "l", xlab = "", ylab = "") # plot of the D2i value
mtext("Taxon", 1, line = 2)
mtext(expression(paste("D"[i]^2)), 2, line = 2)


##### Morphometrics
#### Fore-wing
### Read in data
fore <- readland.tps("dorsal_data/Hama_dorsal_3_aligned.TPS", specID = "ID")
str(fore)
dim(fore) # 50 landmarks with X & Y coordinates, for 71 samples

### conduct generalized Procrustes superimposition
fore_gpa <- gpagen(A = fore, PrinAxes = FALSE, Proj = TRUE, ProcD = TRUE, max.iter = 1e4) # used Procrustes distance for sliding in TPSrelW
plotOutliers(fore_gpa$coords)
fw <- two.d.array(fore_gpa$coords)
dim(fw)
table(rownames(fw))
match(Hama2$tip.label, row.names(fw)) # confirm names in tree and data match

# examine correlations of characters within wings as well
plotAllSpecimens(fore_gpa$coords, mean = TRUE)
str(fore_gpa)

# pdf(file = "FW-pca.pdf", bg = "white")
plotTangentSpace(fore_gpa$coords, label = TRUE) # PC1 40.8%, PC2 25.7%
# dev.off()


### Now calculate mean shape because we can only use one indidivual per species.
groupF <- as.factor(unique(rownames(fw)))

pF <- dim(fore_gpa$coords)[1]
kF <- dim(fore_gpa$coords)[2]
YF <- array(NA, dim = c(pF, kF, length(levels(groupF))))
dimnames(YF)[[3]] <- levels(groupF)
dim(YF)

for(i in 1:length(levels(groupF))){
	grpF <- fw[which(groupF == levels(groupF)[i]), ]
	fooF <- arrayspecs(grpF, pF, kF)
	YF[, , i] <- mshape(fooF)
}
YF
dim(YF)
plot(YF[, 1, ], YF[, 2, ], pch = 19)

YF2d <- two.d.array(YF)
head(YF2d)
dim(YF2d)

# I think that my function worked correctly, but I want to be super duper extra safe, so I'm going to double check. 
velutina <- fw[c(17, 68:70), ]
v1F <- arrayspecs(velutina, p = 50, k = 2)
v2F <- mshape(v1F)

plot(YF[, , "velutina"], pch = 19)
points(v2F[, 1], v2F[, 2], pch = 15) #good.

# pdf(file = "FW-species.pdf", bg = "white")
plotTangentSpace(A = YF, label = TRUE) # PC1 = 37%, PC2 = 28%
# dev.off()

plotGMPhyloMorphoSpace(phy = Hama2, A = YF) 

# Now confirm the set of landmarks that correspond to leading, outer, and trailing edge of the wing. 
plot(YF[, 1, ], YF[, 2, ], pch = 19)
points(YF[1:6, 1, ], YF[1:6, 2, ], col = "red", pch = 19) #edge of wing
points(YF[7:33, 1, ], YF[7:33, 2, ], pch = 19, col = "Dodgerblue") # leading edge
points(YF[34:50, 1, ], YF[34:50, 2, ], pch = 19, col = "dark green")# trailing edge points

## Import forewing centroid size data
fore_cs <- read.delim("dorsal_data/Hama_dorsal_3_centroid.txt", sep = "\t", header = TRUE)
head(fore_cs)
str(fore_cs)
fore_cs
table(fore_cs$species)

## calculate mean centrois size by species
fw_cs <- with(fore_cs, tapply(csize, species, mean))
fw_cs <- fw_cs[-4] #remove arete, which is not in the tree
fw_cs <- fw_cs[Hama2$tip.label] #make sure the labels are in the same order as the tree
fw_cs 


#####
##### Import hind wing data
#####

hind <- readland.tps("ventral_data/Hama_vent_3.tps", specID = "ID")
str(hind)
dim(hind)

hind_gpa <- gpagen(A = hind, PrinAxes = FALSE, Proj = TRUE, ProcD = TRUE, max.iter = 1e4) # used Procrustes distance for sliding in TPSrelW
hw <- two.d.array(hind_gpa$coords)
dim(hw)
table(rownames(hw))
match(Hama2$tip.label, row.names(hw)) # confirm names in tree and data match

groupH <- as.factor(unique(rownames(hw)))

pH <- dim(hind_gpa$coords)[1]
kH <- dim(hind_gpa$coords)[2]
YH <- array(NA, dim = c(pH, kH, length(levels(groupH))))
dimnames(YH)[[3]] <- levels(groupH)
dim(YH)

for(i in 1:length(levels(groupH))){
	grpH <- hw[which(groupH == levels(groupH)[i]), ]
	fooH <- arrayspecs(grpH, pH, kH)
	YH[, , i] <- mshape(fooH)
}
YH
dim(YH)
plot(YH[, 1, ], YH[, 2, ], pch = 19)

YH2d <- two.d.array(YH)
head(YH2d)
dim(YH2d)

plotTangentSpace(A = YH, label = TRUE, warpgrids = TRUE, verbose = FALSE) # PC1 = 60%, PC2 = 14%

# hind wing centroid size
hind_gpa$Csize


hind_cs <- as.data.frame(names(hind_gpa$Csize))

hind_cs$Csize <- (hind_gpa$Csize)
hind_cs <- rename(hind_cs, c("names(hind_gpa$Csize)" = "species"))
hind_cs
table(hind_cs$species)

hw_cs <- with(hind_cs, tapply(Csize, species, mean))
hw_cs

hw_cs <- hw_cs[Hama2$tip.label] #make sure the labels are in the same order as the tree
hw_cs 


#####
##### enter covariate data
#####
## import range data
ranges <- read.csv("Hamadryas_range.csv", header = TRUE, row.names = 1)
ranges$Range <- as.factor(ranges$Range)
str(ranges)
# Small <= 50000; 50001 > Medium < 100000; 100001 > Large 

H.range <- matrix(ranges$Range, dimnames = list(row.names(ranges)))
H.range <- as.factor(H.range)
names(H.range) <- row.names(ranges)
H.range <- H.range[Hama2$tip.label]
H.range

## Import habitat data
cu <- c(NA, 0, 0.5, 0.5, 0, 0.5, 0, 0, 0, 0.5, 0.5, 0, 0, 0 , NA, 1, 1) #0 = understory, 0.5 = mixed, 1 = canopy
cu <- matrix(cu, dimnames = list(rownames(YF2d)))
cu <-cu[Hama2$tip.label, ]
cu1 <- cu[-c(3, 9)]
cu1

Hama.pruned <- treedata(phy = Hama2, data = cu1, sort = TRUE, warnings = TRUE) # Safety check to ensure that the covariate data and tree tips match up
plot(Hama.pruned$phy)



#####
##### The comparative methods analysis
#####

## Forewing
physignal(phy = Hama2, A = YF, iter = 1e4, seed = 876234872) # K = 0.716, P << 0.001

# phylogenetic signal for centroid size
phylosig(tree = Hama2, x = fw_cs, nsim = 1e4, method = "K", test = TRUE) # K = 0.30, P = 0.856

## Hind wing
physignal(phy = Hama2, A = YH, iter = 1e4, seed = 826342) # K = 0.8, P = 0.0018

phylosig(tree = Hama2, x = hw_cs, nsim = 1e4, method = "K", test = TRUE) # K = 0.314, P = 0.81

## ancestral state reconstruction
anc.ML(tree = Hama.pruned$phy, x = cu1, maxit = 1e4, model = "BM")
fastAnc(tree = Hama.pruned$phy, x = cu1, vars = TRUE, CI = TRUE)
# pdf(file = "Hama_FastAnc.pdf", bg = "white")
contMap(Hama.pruned$phy, x = cu1, res = 1000)
# dev.off()

# pdf(file = "Hama_fan.pdf", bg = "white")
contMap(Hama.pruned$phy, x = cu1, res = 1000, type = "fan", legend = FALSE)
# dev.off()





## phlyogenetic morphological covartiation
# leading and trailing edge of FW
LE_TE_pls2 <- phylo.integration(A = YF[7:33, , ], A2 = YF[34:50, , ], phy = Hama2, iter = 1e4, seed = 98234) # r-PLS = 0.867, P = 0.001

# internal and terminal (vein) of HW
HW_pls <- phylo.integration(A = YH[1:6, , ], A2 = YH[7:12, , ],	phy = Hama2, iter = 1e4, seed = 876234) # r-pls = 0.94, P = <<0.001

#accounting for phylogeny, are FW and HW "integrated"?
FW_HW <- phylo.integration(A = YF, A2 = YH, phy = Hama2, iter = 1e4, seed = 987324) # r-pls = 0.85, P = 0.0028


### compare rate of morphological change
## Forewing
# change rate based on range size
fw_range_rate <- compare.evol.rates(phy = Hama2, A = YF, gp = H.range, iter = 1e4) # no sig difference in FW rates based on habitat size 

# change rate based on habitat
fw_habitat_rate <- compare.evol.rates(phy = Hama2, A = YF, gp = cu1, iter = 1e4) # no association with habitat


## hind wing
# rate of HW change based on range
Hw_range_rate <- compare.evol.rates(phy = Hama2, A = YH, gp = H.range, iter = 1e4) # HW rate is sig fast

# rate of HW change based on habitat
Hw_habitat_rate <- compare.evol.rates(phy = Hama2, A = YH, gp = cu1, iter = 1e4) # Yup. HW has different rates based on habitat



## making plots
nspecies <- length(Hama2$tip.label)
cu_color <- cu
cu_color[cu == 0.0] <- "brown"
cu_color[cu == 0.5] <- "orange"
cu_color[cu == 1.0] <- "dark green"


# pdf(file = "Hama_hab_tree.pdf", bg = "white")
plot.phylo(Hama2, font = 3, label.offset = 0.07, edge.width = 3)
points(x = rep(1.035, nspecies), y = 1:nspecies, pch = 19, col = cu_color, cex = 1.5)
# dev.off()

anc.ML(tree = Hama.pruned$phy, x = cu1, maxit = 1e4, model = "BM")
fastAnc(tree = Hama.pruned$phy, x = cu1, vars = TRUE, CI = TRUE)

# pdf(file = "Hama_FastAnc.pdf", bg = "white")
contMap(Hama.pruned$phy, x = cu1, res = 1000)
# dev.off()

# pdf(file = "Hama_fan.pdf", bg = "white")
contMap(Hama.pruned$phy, x = cu1, res = 1000, type = "fan", legend = FALSE)
# dev.off()



# To get rate for whole group (not by subgroup), use a dummy variable.
dummy <- c(rep(0, 17))
dummy <- matrix(dummy, dimnames = list(names(hw_cs)))
dummy <- dummy[Hama2$tip.label, ]
dummy <- as.factor(dummy)
dummy

fw_dummy <- compare.evol.rates(phy = Hama2, A =YF, gp = dummy, iter = 1e4)

hw_dummy <- compare.evol.rates(phy = Hama2, A = YH, gp = dummy, iter = 1e4)



range.color <- as.integer(H.range)
range.color[range.color == "3"] <- "red"
range.color[range.color == "2"] <- "green"
range.color[range.color == "1"] <- "dodgerblue"

# pdf(file = "Hamadryas_phylo2.pdf")
plot.phylo(Hama2, font = 3, label.offset = 0.03, edge.width = 3)
points(rep(0.935, length(Hama2$tip.label)), 1 : length(Hama2$tip.label), pch = 19, col = range.color, cex = 1.5)
points(rep(0.98, length(Hama2$tip.label)), 1 : length(Hama2$tip.label), pch = 19, col = cu_color, cex = 1.5)
legend("topleft", legend = c("Range", "Small", "Medium", "Large", "Habitat", "Ground", "Canopy", "Mixed"), ncol = 2, pch = 19, col = c(0, "red", "green", "dodgerblue", 0, "darkred", "darkgreen", "orange"), bty = "n", pt.cex = 1.5)
segments(x0 = 0.02, x1 = 0.135, y0 = 16.6, y1 = 16.6, lty = 1, lwd = 2)
segments(x0 = 0.239, x1 = 0.36, y0 = 16.6, y1 = 16.6, lty = 1, lwd = 2)
# dev.off()
