
# GFLASSO modelling using the rice population characterized in 
# Genetic analysis of the metabolome exemplified using a rice population
# Gong et al., 2013 (PNAS)

library(readxl)
library(corrplot)
library(pheatmap)
library(RColorBrewer)
library(gflasso)

dir.create("data")
download.file("http://www.pnas.org/highwire/filestream/614522/field_highwire_adjunct_files/4/sd04.xls",
              destfile = "data/data.xls")
download.file("http://www.pnas.org/highwire/filestream/614522/field_highwire_adjunct_files/5/sd05.xls",
              destfile = "data/metabAnnotation.xls")

# Read files
genotype <- data.frame(read_xls("data/data.xls", sheet = 1, skip = 2),
                       row.names = 1)
metab <- data.frame(read_xls("data/data.xls", sheet = 2, skip = 2),
                    row.names = 1)
# The annotation has multiple records per peak (only metabs with mQTL)
annotation <- read_xls("data/metabAnnotation.xls", sheet = 1, skip = 1)

# Use mQTL-positive metabolites, determine overlap
commonMetabs <- intersect(rownames(metab), annotation$Trait)
annCols <- data.frame("id" = commonMetabs, row.names = 1,
                      "Metabolite" = annotation$Metabolite[match(commonMetabs, annotation$Trait)])
metab <- metab[commonMetabs,]
metab <- scale(t(metab)) # transpose and scale, my convention is samples x vars

# Encode SNPs from ILs
# A = Zhenshan97, B = Minghui63, H = heterozygous
ann <- genotype[,1:2]
ann$chr <- factor(ann$chr)
genotype <- genotype[,-c(1:2)] # remove chr and position.Mb. columns
genotype[genotype == "A"] <- 0
genotype[genotype == "H"] <- 1
genotype[genotype == "B"] <- 2
genotype <- t(data.matrix(genotype)) # transpose, my convention is samples x vars

# Ensure the samples match across
commonSamples <- intersect(rownames(metab),
                           rownames(genotype))
metab <- metab[commonSamples,]
genotype <- genotype[commonSamples,]

# Visualize genotype diversity
pdf("genHeat.pdf", 6, 10)
pheatmap(genotype, color = colorRampPalette(brewer.pal(n = 3, name =
                                      "Greys"))(3),
         cluster_rows = F, cluster_cols = F, border_color = NA,
         show_colnames = F, legend_breaks = c(.25, 1, 1.75),
         legend_labels = c("ZS97", "Het", "MH63"))
dev.off()

## GFLASSO
R <- cor(metab, method = "spearman")
rownames(R) <- gsub(annCols$Metabolite, pattern = "’", replacement = "")
pdf("corrMetab.pdf", 10, 6)
corrplot(R)
dev.off()
# Rewrite rownames
rownames(R) <- colnames(R)

set.seed(100)
system.time(cv <- cv_gflasso(genotype, metab, R, k = 5, times = 4, nCores = 20,
                 additionalOpts = list(delta_conv = 1e-5,
                                       iter_max = 1e5)))

# Plot cv results
pdf("cvFigure.pdf", 5, 5)
cv_plot_gflasso(cv)
dev.off()

# Fit final model with optimal pars
finalModel <- gflasso(X = genotype, Y = metab, R = R, opts = list(lambda = cv$optimal$lambda,
                                                      gamma = cv$optimal$gamma,
                                                      delta_conv = 1e-5,
                                                      iter_max = 1e5))
colnames(finalModel$B) <- colnames(metab)

pdf("results.pdf", 8, 10)
# Red = Up in Minghui63; Blue = Up in Zhenshan97
pheatmap(finalModel$B, breaks = seq(-1, 1, length.out = 100), show_rownames = F,
         annotation_row = ann, cluster_rows = F, border_color = NA)
# Plot again, but using the metabolite names
pheatmap(finalModel$B, breaks = seq(-1, 1, length.out = 100), show_rownames = F,
         annotation_row = ann, cluster_rows = F, border_color = NA,
         labels_col = gsub(annCols$Metabolite, pattern = "’", replacement = ""))
dev.off()
# Draw red circles around mQTLs from the study! Finally, save model for future predictions
# You can use predict_gflasso on new samples
save(finalModel, file = "finalModel.RData")