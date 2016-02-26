# Hamadryas 

library("geomorph") 
library("phytools") 
library("geiger") 
library("dplyr")
library("plyr")

setwd("~/Desktop/Projects/Hamadryas/Hama_data")

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

## conduct generalized Procrustes superimposition
fore_gpa <- gpagen(A = fore, Proj = TRUE, ProcD = TRUE, ShowPlot = TRUE) # used Procrustes distance for sliding in TPSrelW
fw <- two.d.array(fore_gpa$coords)
dim(fw)
table(rownames(fw))
match(Hama2$tip.label, row.names(fw)) # confirm names in tree and data match

plotTangentSpace(fore_gpa$coords, label = TRUE)

# plot lda by species
fw.pc <- prcomp(fw)
summary(fw.pc)

fw.pc.shape <- cbind(rownames(fw), fw.pc$x[, 1:67]) # remove final 4 dimensions because they are empty.
head(fw.pc.shape)
dim(fw.pc.shape)

fw.lda <- lda(as.matrix(fw.pc.shape[, 2:68]), as.factor(fw.pc.shape[, 1]), method = "mle")



## Now calculate mean shapet because we can only use one indidivual per species.

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

physignal(phy = Hama2, A = YF, iter = 1e4)

plotTangentSpace(A = YF, label = TRUE, warpgrids = TRUE, verbose = FALSE) # PC1 = 37%, PC2 = 28%

plotGMPhyloMorphoSpace(phy = Hama2, A = YF) 

# Now confirm the set of landmarks that correspond to leading, outer, and trailing edge of the wing. 
plot(YF[, 1, ], YF[, 2, ], pch = 19)
points(YF[1:6, 1, ], YF[1:6, 2, ], col = "red", pch = 19) #edge of wing
points(YF[7:33, 1, ], YF[7:33, 2, ], pch = 19, col = "Dodgerblue") # leading edge
points(YF[34:50, 1, ], YF[34:50, 2, ], pch = 19, col = "dark green")# trailing edge points

## Import forewing centroid size data
fore_cs <- read.delim("dorsal_data/Hama_dorsal_3_centroid.txt", sep = "\t", header = TRUE)
head(fore_cs)
fore_cs
table(fore_cs$species)

t1 <- fore_cs %>% group_by(species) %>% summarize(mean_csize = mean(csize))
t1
fw_cs <- matrix(t1$mean_csize, dimnames = list(t1$species))
fw_cs <- fw_cs[-4, ]

fw_cs <- fw_cs[Hama2$tip.label]
fw_cs 

# phylogenetic signal for centroid size
physignal(phy = Hama2, A = fw_cs, iter = 1e4, method = "Kmult")













### Import hind wing data
hind <- readland.tps("ventral_data/Hama_vent_3.tps", specID = "ID")
str(hind)
dim(hind)

hind_gpa <- gpagen(A = hind, Proj = TRUE, ProcD = TRUE, ShowPlot = TRUE) # used Procrustes distance for sliding in TPSrelW
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

# I think that my function worked correctly, but I want to be super duper extra safe, so I'm going to double check. 
velutinaH <- hw[c(17, 68:70), ]
vH1 <- arrayspecs(velutinaH, p = 50, k = 2)
vH2 <- mshape(vH1)

plot(YH[, , "velutina"], pch = 19)
points(vH2[, 1], vH2[, 2], pch = 15) #good.


physignal(phy = Hama2, A = YH, iter = 1e4, method = "Kmult")

plotTangentSpace(A = YH, label = TRUE, warpgrids = TRUE, verbose = FALSE) # PC1 = 37%, PC2 = 28%

# hind wing centroid size
hind_gpa$Csize


hind_cs <- as.data.frame(names(hind_gpa$Csize))

hind_cs$Csize <- (hind_gpa$Csize)
hind_cs <- rename(hind_cs, c("names(hind_gpa$Csize)" = "species"))
hind_cs
table(hind_cs$species)

t2 <- hind_cs %>% group_by(species) %>% summarize(mean_csize = mean(Csize))
t2



hw_cs <- matrix(t1$mean_csize, dimnames = list(t2$species))
hw_cs <- fw_cs[-4, ]

fw_cs <- fw_cs[Hama2$tip.label]
fw_cs 








### enter covariate data
# import range data
ranges <- read.csv("Hamadryas_range.csv", header = TRUE, row.names = 1)
ranges$Range <- as.factor(ranges$Range)
str(ranges)

H.range <- matrix(ranges$Range, dimnames = list(row.names(ranges)))
H.range <- as.factor(H.range)
names(H.range) <- row.names(ranges)
H.range
H.range <- H.range[Hama2$tip.label]

# Import habitat data
cu <- c(NA, 0, 0.5, 0.5, 0, 0.5, 0, 0, 0, 0.5, 0.5, 0, 0, 0 , NA, 1, 1) #0 = understory, 0.5 = mixed, 1 = canopy

cu <- matrix(cu, dimnames = list(rownames(YF2d)))

cu <-cu[Hama2$tip.label, ]
cu1 <- cu[-c(3, 9)]
cu1

Hama.pruned <- treedata(phy = Hama2, data = cu1, sort = TRUE, warnings = TRUE) Safety check to ensure that the covariate data and tree tips match up
plot(Hama.pruned$phy)

anc.ML(tree = Hama.pruned$phy, x = cu1, maxit = 1e4, model = "BM")
fastAnc(tree = Hama.pruned$phy, x = cu1, vars = TRUE, CI = TRUE)
# pdf(file = "Hama_FastAnc.pdf", bg = "white")
contMap(Hama.pruned$phy, x = cu1, res = 1000)
# dev.off()

# pdf(file = "Hama_fan.pdf", bg = "white")
contMap(Hama.pruned$phy, x = cu1, res = 1000, type = "fan", legend = FALSE)
# dev.off()


cu_color <- cu
cu_color[cu == 0.0] <- "brown"
cu_color[cu == 0.5] <- "orange"
cu_color[cu == 1.0] <- "dark green"

nspecies <- length(Hama2$tip.label)

cuf <- matrix(cu, dimnames = list(row.names(ranges)))
cuf <- as.factor(cuf)
names(cuf) <- row.names(ranges)
cuf
cuf <- H.range[Hama2$tip.label]



# compare canopy and understory rates

(LE_TE_pls2 <- phylo.pls(A = YF[7:33, , ], A2 = YF[34:50, , ], phy = Hama2, iter = 1e4)) # r-PLS 0.867

(fw_habitat_rate2 <- compare.evol.rates(phy = Hama2, A = Y, gp = cuf, iter = 1e4)) 

(fw_range_rate2 <- compare.evol.rates(phy = Hama2, A = Y, gp = H.range, iter = 1e4))




## test for phylogenetic signal
# forewing
(fw_ps <- physignal(phy = Hama2, A = YF2d, iter = 1e4)) # fore-wing shape is not randomly distributed through the tree

#centroid size
(fw_cs_ps <- physignal(phy = Hama2, A = fw_cs, iter = 1e4)) # fore-wing centroid size is not randomly distributed through the tree

# PCA on forewing shape data
plotTangentSpace(A = fore_gpa$coords, label = TRUE, warpgrids = TRUE, verbose = FALSE) # PC1 46%, PC2 24%

# now to ask if there is evidence for convergent evolution. If there is, taxa on different branches will converge towards the same point in tangent space (which is a PCA on shape data).

plotGMPhyloMorphoSpace(phy = Hama2, A = fore_gpa$coords) 
# nb, for the figures in the paper I used the above code as a starting point and modified it to make it a bit easier to read.


# compare rates of canopy vs. understory?
# allometry
(fw_allo <- procD(A = fore_gpa$coords, sz = fw_cs, label = Hama2$tip.label, method = "RegScore", mesh = TRUE, iter = 1e4)) # the difference in size is not associated with a difference in shape. But this does not account for phylogeny
















##### Now the real comparative stuff

# pdf(file = "Hama_hab_tree.pdf", bg = "white")
plot.phylo(Hama2, font = 3, label.offset = 0.07, edge.width = 3)
points(x = rep(1.035, nspecies), y = 1:nspecies, pch = 19, col = cu_color, cex = 1.5)
# dev.off()


Hama.pruned <- treedata(phy = Hama2, data = cu1, sort = TRUE, 
	warnings = TRUE) # removes alicia and julitta from data set (we have no ecological information for these)
plot(Hama.pruned$phy)

anc.ML(tree = Hama.pruned$phy, x = cu1, maxit = 1e4, model = "BM")
fastAnc(tree = Hama.pruned$phy, x = cu1, vars = TRUE, CI = TRUE)

# pdf(file = "Hama_FastAnc.pdf", bg = "white")
contMap(Hama.pruned$phy, x = cu1, res = 1000)
# dev.off()

# pdf(file = "Hama_fan.pdf", bg = "white")
contMap(Hama.pruned$phy, x = cu1, res = 1000, type = "fan", legend = FALSE)
# dev.off()


##########
# Accounting for phlyogeny, are fore- and hind-wing "integrated"?
(F_H_pls <- phylo.pls(A1 = YF, A2 = YH, phy = Hama2, warpgrids = TRUE, iter = 1e4)) # even when correcting for phylogeny there is no association

# examine correlations of characters within wings as well
plotAllSpecimens(fore_gpa$coords, mean = TRUE)
str(fore_gpa)

plot(fore_gpa$coords[, 1, ], fore_gpa$coords[, 2, ], pch = 19, xlim = c(-0.21, 0.21), ylim = c(-0.1, 0.16))
# vein points
points(fore_gpa$coords[1:6, 1, ], fore_gpa$coords[1:6, 2, ], pch = 19, col = "red")
# leading edge points
points(fore_gpa$coords[7:33, 1, ], fore_gpa$coords[7:33, 2, ], pch = 19, col = "Dodgerblue")
# trailing edge points
points(fore_gpa$coords[34:50, 1, ], fore_gpa$coords[34:50, 2, ], pch = 19, col = "dark green")

# comparing hw internal and wing vein landmarks
plot(hind_gpa$coords[, 1, ], hind_gpa$coords[, 2, ], pch = 19, xlim = c(-0.4, 0.4), ylim = c(-0.3, 0.3))
points(hind_gpa$coords[1:6, 1, ], hind_gpa$coords[1:6, 2, ], col = "dodgerblue", pch = 19)
points(hind_gpa$coords[7:12, 1, ], hind_gpa$coords[7:12, 2, ], col = "red", pch = 19)

# we only want to compare "like" landmarks, such as type I only comapred to type I, and semi-landmarks compared only to semi-landmarks

# compare leading edge and trailing edge of the wing
(LE_TE_pls <- phylo.pls(A1 = fore_gpa$coords[7:33, , ], 
	A2 = fore_gpa$coords[34:50, , ], phy = Hama2, warpgrids = TRUE, 
	iter = 1e4))

(HW_pls <- phylo.pls(A1 = hind_gpa$coords[1:6, , ], A2 = 
	hind_gpa$coords[7:12, , ],	phy = Hama2, warpgrids = TRUE, 
	iter = 1e4))








# forewing rates
(fw_habitat_rate <- compare.evol.rates(phy = Hama2, 
	A = fore_gpa$coords, gp = cuf, iter = 1e4))

# hindwing rates
(hw_habitat_rate <- compare.evol.rates(phy = Hama2, 
	A = hind_gpa$coords, gp = cuf, iter = 1e4))

# To get rate for whole group (not by subgroup), use a dummy variable.
dummy <- c(NA, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, NA, 0, 0)
dummy <- matrix(dummy, dimnames = list(names(hind_gpa$Csize)))
dummy <- dummy[Hama2$tip.label, ]
dummy <- dummy[-c(3, 9)]
dummy <- as.factor(dummy)

(fw_dummy <- compare.evol.rates(phy = Hama2, A = fore_gpa$coords, 
	gp = dummy, iter = 1e4))

(hw_dummy <- compare.evol.rates(phy = Hama2, A = hind_gpa$coords, 
	gp = dummy, iter = 1e4))


# What about the size of species ranges? Hypothesis that morphological evolution is accelerated in species with reduced ranges.  
# Small <= 50000; 50001 > Medium < 100000; 100001 > Large 
ranges <- read.csv("Hamadryas_range.csv", header = TRUE, row.names = 1)
ranges$Range <- as.factor(ranges$Range)
str(ranges)

H.range <- matrix(ranges$Range, dimnames = list(row.names(ranges)))
H.range <- as.factor(H.range)
names(H.range) <- row.names(ranges)
H.range
H.range <- H.range[Hama2$tip.label]


(fw_range_rate <- compare.evol.rates(phy = Hama2, A = fore_gpa$coords, 
	gp = H.range, iter = 1e4)) #runs in support to the notion that reduced range species have accelerated rates of Evolution

(hw_range_rate <- compare.evol.rates(phy = Hama2, A = hind_gpa$coords, 
	gp = H.range, iter = 1e4)) # for hindwing the small range has the lowest rate. 

# Figure 2
Hama.range <- c("Large", "Small", "Small", "Medium", "Large", "Large", "Medium", "Medium", "Small", "Large", "Large", "Large", "Large", "Large", "Large", "Small", "Large")
Hama.range <- matrix(Hama.range, dimnames = list(Hama2$tip.label))


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


#####
