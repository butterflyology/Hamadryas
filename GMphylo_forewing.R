x <- two.d.array(fore_gpa$coords)

N <- length(Hama2$tip.label)

x <- x[Hama2$tip.label, ]
names <- row.names(x)
anc.states <- NULL
    
for (i in 1:ncol(x)) {
    options(warn = -1)
    tmp <- as.vector(fastAnc(Hama2, x[, i]))
    anc.states <- cbind(anc.states, tmp)
}
colnames(anc.states) <- NULL
row.names(anc.states) <- 1:length(tmp)
all.data <- rbind(x, anc.states)
pcdata <- prcomp(all.data)$x

cu <- c(2, 0, 0.5, 0.5, 0, 0.5, 0, 0, 0, 0.5, 0.5, 0, 0, 0 , 2, 1, 1)

cu <- matrix(cu, dimnames = list(names(hind_gpa$Csize)))
cu <- cu[Hama2$tip.label, ]


cu_color <- cu
cu_color[cu == 0] <- "brown"
cu_color[cu == 0.5] <- "orange"
cu_color[cu == 1] <- "dark green"
cu_color[cu == 2] <- "black"

# might want to rotate 
# pdf(file = "FW_GMphylo_2a_2.pdf", bg = "white")
plot(pcdata, type = "n", xlim = c(-0.051, 0.051), ylim = c(-0.035, 0.06), las = 1, xlab = "PC1 (37%)", ylab = "PC2 (28%)")

# plot(pcdata, type = "n", xlim = c(-0.08, 0.08), ylim = c(-0.05, 0.05), las = 1, asp = 1) force the 1:1 aspect ratio

for (i in 1:nrow(Hama2$edge)) {
	lines(pcdata[(Hama2$edge[i, ]), 1], pcdata[(Hama2$edge[i, ]), 2], 	type = "l", pch = 19, col = "black", lwd = 3)
}
    
points(pcdata[1:N, ], pch = 19, col = cu_color, cex = 2)
points(pcdata[(N + 1):nrow(pcdata), ], pch = 21, bg = "white", cex = 1.25)

legend("left", legend = c("Ground", "Canopy", "Mixed", "Unknown"), pch = 19, col = c("brown", "dark green", "orange", "black"), bty = "n", pt.cex = 1.5)

# text(pcdata[1:17, 1], pcdata[1:17, 2], rownames(pcdata[1:17, ]), cex = 1, adj = c(-0.05, 1.8), font = 3)
# dev.off()
