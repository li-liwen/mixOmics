## Customizing individual node colors with color.node.individual
# Useful for highlighting specific nodes (e.g., based on log fold change)
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene
nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)

# Calculate log fold change for each variable (binary genotype: WT vs KO)
logfc.lipid <- rowMeans(nutrimouse$lipid[nutrimouse$genotype == "KO", ]) - 
               rowMeans(nutrimouse$lipid[nutrimouse$genotype == "WT", ])
logfc.gene <- rowMeans(nutrimouse$gene[nutrimouse$genotype == "KO", ]) - 
              rowMeans(nutrimouse$gene[nutrimouse$genotype == "WT", ])

# Create color gradient: blue (negative) to red (positive)
color.range <- range(c(logfc.lipid, logfc.gene))
cols.lipid <- colorRampPalette(c("blue", "white", "red"))(100)[
    findInterval(logfc.lipid, seq(color.range[1], color.range[2], length.out = 100))]
cols.gene <- colorRampPalette(c("blue", "white", "red"))(100)[
    findInterval(logfc.gene, seq(color.range[1], color.range[2], length.out = 100))]

# Create named color vector for individual nodes
color.node.individual <- c(setNames(cols.lipid, rownames(nutrimouse$lipid)),
                           setNames(cols.gene, rownames(nutrimouse$gene)))

# Plot network with custom individual colors
jpeg('example-custom-colors-network.jpeg')
network(nutri.res, comp = 1:3, cutoff = 0.6, color.node.individual = color.node.individual)
dev.off()

## network representation for objects of class 'rcc'
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene
nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)

\dontrun{
# may not work on the Linux version, use Windows instead
# sometimes with Rstudio might not work because of margin issues,
# in that case save it as an image
jpeg('example1-network.jpeg', res = 600, width = 4000, height = 4000)
network(nutri.res, comp = 1:3, cutoff = 0.6)
dev.off()

## Changing the attributes of the network

# sometimes with Rstudio might not work because of margin issues,
# in that case save it as an image
jpeg('example2-network.jpeg')
network(nutri.res, comp = 1:3, cutoff = 0.45,
color.node = c("mistyrose", "lightcyan"),
shape.node = c("circle", "rectangle"),
color.edge = color.jet(100),
lty.edge = "solid", lwd.edge = 2,
show.edge.labels = FALSE)
dev.off()


## interactive 'cutoff' - select the 'cutoff' and "see" the new network
## only run this during an interactive session
if (interactive()) {
    network(nutri.res, comp = 1:3, cutoff = 0.55, interactive = TRUE)
}
dev.off()

## network representation for objects of class 'spls'
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50),
keepY = c(10, 10, 10))

# sometimes with Rstudio might not work because of margin issues,
# in that case save it as an image
jpeg('example3-network.jpeg')
network(toxicity.spls, comp = 1:3, cutoff = 0.8,
color.node = c("mistyrose", "lightcyan"),
shape.node = c("rectangle", "circle"),
color.edge = color.spectral(100),
lty.edge = "solid", lwd.edge =  1,
show.edge.labels = FALSE, interactive = FALSE)
dev.off()
}