# Hamadryas 

library("geomorph") # v2.1.1
library("phytools") # v0.4-31, loads ape v3.2 
library("geiger") # v2.0.3
setwd("~/Desktop/Projects/Hamadryas/Hama_data")

# save(list = ls(), file = "Hama_data_1.R")
load("Hama_data_1.R")
# packageDescription("ape")$Version

Hama <- read.nexus("../Tree_data/Hamydryas_ml.tre")
plot(Hama)
is.ultrametric(Hama)

# make tree ultrametric, lamdba is the smoothing parameter. 
lam <- 10^(-1:6)
cv <- sapply(lam, function(x) sum(attr(chronopl(Hama, lambda = x, CV = TRUE), "D2")))
plot(lam, cv, pch = 19, ylab = "cross-validation score", xlab = expression(paste(lambda)), las = 1, cex = 1.5) # lowest CV 

Hama2 <- chronopl(phy = Hama, lambda = 0.1, CV = TRUE, eval.max = 1e3, iter.max = 1e4)
is.ultrametric(Hama2)

plot(attr(Hama2, "D2"), type = "l", xlab = "", ylab = "") #plot of the D2i value
mtext("Taxon", 1, line = 2)
mtext(expression(paste("D"[i]^2)), 2, line = 2)

# pdf(file = "Hamadryas_phylo.pdf")
plot.phylo(Hama2, font = 3, label.offset = 0.03, edge.width = 3)
# dev.off()




##### Morphometrics

#### Fore-wing
# read in data
fore <- readland.tps("dorsal_data/Hama_dorsal_aligned.tps", specID = "ID")
str(fore)
dim(fore) # 50 landmarks with X & Y coordinates, for 17 samples

fore_gpa <- gpagen(A = fore, Proj = TRUE, ProcD = TRUE, ShowPlot = TRUE) # used Procrustes distance for sliding in TPSrelW
fw <- two.d.array(fore_gpa$coords)
match(Hama2$tip.label, row.names(fw)) # confirm names in tree and data match

# Because I constructed the sliders in a different program the centroid sizes calculated by gpagen are not correct. Here I import them manually (lame, I know). 
fore_cs <- c(3.33294570842831e+3, 3.16568469860016e+3, 3.39261539540194e+3, 3.78430193381669e+3, 3.48021469444961e+3, 3.33777956879712e+3, 2.68538836886851e+3,3.17658593061965e+3, 
3.39194656501948e+3, 3.61375985829088e+3, 3.41355030541168e+3, 2.85453923363221e+3, 4.05270026183310e+3, 3.10518493006033e+3, 3.01038815277082e+3, 3.56187305827559e+3, 3.64253669986741e+3) 
fw_cs <- matrix(fore_cs, dimnames = list(names(fore_gpa$Csize)))

## test for phylogenetic signal
# fore-wing
(fw_ps <- physignal(phy = Hama2, A = fw, iter = 5e3, method = "Kmult")) # fore-wing shape is not randomly distributed through the tree

#centroid size
(fw_cs_ps <- physignal(phy = Hama2, A = fw_cs, iter = 5e3, method = "Kmult")) # fore-wing centroid size is not randomly distributed through the tree

# PCA on fore-wing shape data
plotTangentSpace(A = fore_gpa$coords, label = TRUE, warpgrids = TRUE, verbose = FALSE) # PC1 46%, PC2 24%

# now to ask if there is evidence for convergent evolution. If there is, taxa on different branches will converge towards the same point in tangent space (which is a PCA on shape data).
plotGMPhyloMorphoSpace(phy = Hama2, A = fore_gpa$coords, labels = TRUE) #nb, for the figures in the paper I used the above code as a starting point and modified it to make it a bit easier to read.





# compare rates of canopy vs. understory?
# Hypothetical example
c_u <- factor(c(0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))
names(c_u) <- Hama2$tip.label
compare.evol.rates(phy = Hama2, A = fore_gpa$coords, gp = c_u, iter = 1e3)

# compare morphological disparity for groups, canopy vs. understory. 

# allometry
plotAllometry(A = fore_gpa$coords, sz = fw_cs, label = Hama2$tip.label, method = "RegScore", mesh = TRUE, iter = 5e3) # the difference in size is not associated with a difference in shape





#### Hind-wing
# read in data
hind <- readland.tps("ventral_data/Hama_hind_struc.tps", specID = "ID")
str(hind)
dim(hind) # 12 landmarks with X & Y for 17 individuals

hind_gpa <- gpagen(A = hind, Proj = TRUE, ProcD = TRUE, ShowPlot = TRUE)
hw <- two.d.array(hind_gpa$coords)
match(Hama2$tip.label, row.names(hw))
match(row.names(fw), row.names(hw))

hw_cs <- matrix(hind_gpa$Csize, dimnames = list(names(hind_gpa$Csize)))

# phylogenetic signal
physignal(phy = Hama2, A = hw, iter = 5e3, method = "Kmult") # hw shape is not randomly distributed through the tree

physignal(phy = Hama2, A = hw_cs, iter = 5e3, method = "Kmult") # hindwing centroid size (as a proxy for size) is randomly distributed through the phylogeny, though shape is not.

plotTangentSpace(A = hind_gpa$coords, label = TRUE) # PC1 (54%), PC2 (16%)

plotGMPhyloMorphoSpace(phy = Hama2, A = hind_gpa$coords, labels = TRUE)






##### Now the real comparative stuff
cu <- c(NA, 0, 0.5, 0.5, 0, 0.5, 0, 0, 0, 0.5, 0.5, 0, 0, 0 , NA, 1, 1)

cu <- matrix(cu, dimnames = list(names(hind_gpa$Csize)))
cu <- cu[Hama2$tip.label, ]
cu1 <- cu[-c(3, 9)]
cu1

cu_color <- cu
cu_color[cu == 0] <- "brown"
cu_color[cu == 0.5] <- "orange"
cu_color[cu == 1] <- "dark green"

nspecies <- length(Hama2$tip.label)

# pdf(file = "Hama_hab_tree.pdf", bg = "white")
plot.phylo(Hama2, font = 3, label.offset = 0.07, edge.width = 3)
points(x = rep(1.035, nspecies), y = 1:nspecies, pch = 19, col = cu_color, cex = 1.5)
# dev.off()


Hama.pruned <- treedata(phy = Hama2, data = cu1, sort = TRUE, warnings = TRUE) # removes alicia and julitta from data set (we have no ecological information for these)
plot(Hama.pruned$phy)


anc.ML(tree = Hama.pruned$phy, x = cu1, maxit = 1e4, model = "BM")
fastAnc(tree = Hama.pruned$phy, x = cu1, vars = TRUE, CI = TRUE)


# pdf(file = "Hama_FastAnc.pdf", bg = "white")
contMap(Hama.pruned$phy, x = cu1, res = 1000)
# dev.off()



# Accounting for phlyogeny, are fore- and hind-wing "integrated"?
phylo.pls(A1 = fore_gpa$coords, A2 = hind_gpa$coords, phy = Hama2, warpgrids = TRUE, iter = 5e3) # even when correcting for phylogeny there is no association


# compare canopy and understory rates
summary(cu1)
cuf <- as.factor(cu)



compare.evol.rates(phy = Hama2, A = fore_gpa$coords, gp = cuf, iter = 5e3)

compare.evol.rates(phy = Hama2, A = hind_gpa$coords, gp = cuf, iter = 5e3)

# To get rate for whole group (not by subgroup), use a dummy variable.
dummy <- c(NA, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 , NA, 0, 0)

dummy <- matrix(dummy, dimnames = list(names(hind_gpa$Csize)))
dummy <- dummy[Hama2$tip.label, ]
dummy <- dummy[-c(3, 9)]
dummy <- as.factor(dummy)

compare.evol.rates(phy = Hama2, A = fore_gpa$coords, gp = dummy, iter = 5e3)

compare.evol.rates(phy = Hama2, A = hind_gpa$coords, gp = dummy, iter = 5e3)


# What about the size of species ranges? Hypothesis that morphological evolution is accelerated in species with reduced ranges.  
ranges <- read.csv("Hamadryas_range.csv", header = TRUE, row.names = 1)
ranges$Range <- as.factor(ranges$Range)
str(ranges)

H.range <- matrix(ranges$Range, dimnames = list(row.names(ranges)))
H.range <- as.factor(H.range)
names(H.range) <- row.names(ranges)
H.range

compare.evol.rates(phy = Hama2, A = fore_gpa$coords, gp = H.range, iter = 5e3)

compare.evol.rates(phy = Hama2, A = hind_gpa$coords, gp = H.range, iter = 5e3) # for hindwing the small range has the lowest rate. 
# Try this with more precise factors (better quantification)


# sound producing vs. mute rates


which(fore_gpa == "alicia") # alica is #1, julitta is 15
fore_gpa_1 <- fore_gpa$coords[, , -1]
fore_gpa_1 <- fore_gpa_1[, , -14]

fore_gpa.2d <- two.d.array(fore_gpa_1)

fact <- as.factor(sort(cu1))
str(fact)
# morphological disparity for groups
morphol.disparity(A = fore_gpa_1, groups = fact, iter = 5e3) 


procD.lm(fore_gpa.2d ~ fact, iter = 5e3, RRPP = TRUE)
# keep in mind we will have limited power because of the small size of our data set. 

# phylogenetic regression of shape on ecological variable? Not for GM data













# options(error = recover) # c to escape, error = normal to terminate
# define.modules(fw, nmodules = 2)
# phylogenetic regression for shape based on behavior / ecology? 




# compare integration of leading edge of forewing to trailing edge of forewing
morphol.integr(fore_gpa$coords[7:33, , ], fore_gpa$coords[34:50, , ], method = "RV", iter = 5e3) # the leading edge of forewing seems to go with a trailing edge. This is innane and just means that forwing shape differs among species. Will delete is later versions of code. 




# THIS DOES NOT ACCOUNT FOR PHYLOGENY. Integaration between two modules (2D vs. 3D array should give same results) using PLS and RV
morphol.integr(A1 = fore_gpa$coords, A2 = hind_gpa$coords, method = "PLS", warpgrids = TRUE, iter = 5e3, label = NULL) # fore-wing and hind-wing shape do not appear to be significantly correlated 

morphol.integr(A1 = fore_gpa$coords, A2 = hind_gpa$coords, method = "RV", warpgrids = TRUE, iter = 5e3, label = NULL)

# however, once we correct for phylogeny we find that this means shit p = 0.5

two.b.pls(A1 = fw, A2 = hw, iter = 1e4, warpgrids = TRUE, verbose = FALSE)