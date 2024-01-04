#Figure S4C- Volcano plot
res <- read.csv("CDKIMutant_vs_DMSOMutant.csv", row.names = 1)
res_plot <- res[order(res$padj), ]
res_plot$col  <- 'gray40'
res_plot$col[res_plot$log2FoldChange > 0 & res_plot$padj < 0.05] <- 'blue'
res_plot$col[res_plot$log2FoldChange < 0 & res_plot$padj < 0.05] <- 'orangered'
par(mar = c(4, 4.4, 5, 10), xpd = TRUE)
plot(res_plot$log2FoldChange, -log10(res_plot$padj), col = res_plot$col, pch = 19, cex =0.25, xlab = expression(log[2]("FC")), ylab = expression(-log[10]("FDR")), cex.lab=1.2, cex.axis=1.2, cex.main=1.0, cex.sub=1.2, frame.plot = FALSE, box(bty="l"))
legend("topright", inset = c(-0.25, 0), legend = c("Not Significant", "Upregulated", "Downregulated"), col = c("gray", "orangered", "blue"), lty = 2:4, cex = 0.9, pch = 19, bty="n")