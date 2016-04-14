# hind wing 2d array = YH2d
N <- length(rownames(YH2d))

x <- YH2d[Hama2$tip.label, ]
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

cu.h <- c(0, 0, 2, 0, 0.5, 0, 0, 0, 2, 0, 0, 0.5, 1, 1 , 0.5, 0.5, 0.5)

cu.h <- matrix(cu.h, dimnames = list(rownames(x)))
cu.h <- cu.h[Hama2$tip.label, ]


cu_color <- cu.h
cu_color[cu.h == 0] <- "brown"
cu_color[cu.h == 0.5] <- "orange"
cu_color[cu.h == 1] <- "dark green"
cu_color[cu.h == 2] <- "black"


# pdf(file = "../Hama_figures/HW_GMphylo_2b.pdf", bg = "white")
plot(pcdata, type = "n", xlim = c(-0.13, 0.1), ylim = c(-0.04, 0.05), las = 1, xlab = "PC1 (60%)", ylab = "PC2 (14%)")

# plot(pcdata, type = "n", xlim = c(-0.08, 0.08), ylim = c(-0.05, 0.05), las = 1, asp = 1) force the 1:1 aspect ratio

for (i in 1:nrow(Hama2$edge)) {
	lines(pcdata[(Hama2$edge[i, ]), 1], pcdata[(Hama2$edge[i, ]), 2], 	type = "l", pch = 19, col = "black", lwd = 3)
}
    
points(pcdata[1:N, ], pch = 19, col = cu_color, cex = 2)
points(pcdata[(N + 1):nrow(pcdata), ], pch = 21, bg = "white", cex = 1.25)

text(pcdata[1:17, 1], pcdata[1:17, 2], rownames(pcdata[1:17, ]), cex = 1, adj = c(-0.05, 1.8), font = 3)
# dev.off()
