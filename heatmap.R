load("pval.mtx2.RData")

pval.mtx2 <- pval.mtx2[, grep("kettunen", colnames(pval.mtx2))]

library(ggplot2)
library(reshape2)
library(gplots)

colnames(pval.mtx2) <-
    sapply(colnames(pval.mtx2), function(s) {a <- strsplit(s, "_")[[1]];
                                             paste(a[c(-1,-length(a))], collapse = "_")})

## heatmap.2(pmax(pval.mtx2,-24), Rowv=NULL,Colv=NULL,
##           col = rev(rainbow(20*10, start = 0/6, end = 4/6)),
##           scale="none",
##           margins=c(3,0), # ("margin.Y", "margin.X")
##           trace='none',
##           symkey=FALSE,
##           symbreaks=FALSE,
##           dendrogram='none',
##           density.info='none',
##           denscol="black",
##           keysize=0.25, "left.margin",
##           #( "bottom.margin", "left.margin", "top.margin", "left.margin" )
##           key.par=list(mar=c(3.5,0,3,0)),
##           # lmat -- added 2 lattice sections (5 and 6) for padding
##           lmat=rbind(c(5, 4, 2), c(6, 1, 3)), lhei=c(2.5, 5), lwid=c(1, 10, 1))

load("beta_samps.RData")

beta_q5 <- apply(beta_samps, 1, quantile, probs = 0.05)
beta_q95 <- apply(beta_samps, 1, quantile, probs = 0.95)
beta_q50 <- apply(beta_samps, 1, quantile, probs = 0.5)

ordered.SNPs <- rownames(beta_samps)[order(beta_q50)]

pval.mtx2 <- pval.mtx2[ordered.SNPs, ]
beta_q5 <- beta_q5[ordered.SNPs]
beta_q50 <- beta_q50[ordered.SNPs]
beta_q95 <- beta_q95[ordered.SNPs]

library(reshape2)
df <- melt(as.matrix(pval.mtx2))

traits <- c(
    "vldl_d",
    "xxl_vldl_l", "xxl_vldl_p", "xxl_vldl_pl", "xxl_vldl_tg",
    "xl_vldl_l", "xl_vldl_p", "xl_vldl_pl", "xl_vldl_tg",
    "l_vldl_c", "l_vldl_ce", "l_vldl_fc", "l_vldl_l", "l_vldl_p", "l_vldl_pl", "l_vldl_tg",
    "m_vldl_c", "m_vldl_ce", "m_vldl_fc", "m_vldl_l", "m_vldl_p", "m_vldl_pl", "m_vldl_tg",
    "s_vldl_fc", "s_vldl_l", "s_vldl_p", "s_vldl_pl", "s_vldl_tg",
    "xs_vldl_l", "xs_vldl_p", "xs_vldl_pl", "xs_vldl_tg",
    "apob", "ldl_c", "ldl_d",
    "l_ldl_c", "l_ldl_ce", "l_ldl_fc", "l_ldl_l", "l_ldl_p", "l_ldl_pl",
    "m_ldl_c", "m_ldl_ce", "m_ldl_l", "m_ldl_p", "m_ldl_pl",
    "s_ldl_c", "s_ldl_l", "s_ldl_p",
    "idl_c", "idl_fc", "idl_l", "idl_p", "idl_pl", "idl_tg",
    "apoa1", "hdl_c", "hdl_d",
    "l_hdl_c",
    "xl_hdl_c", "xl_hdl_ce", "xl_hdl_fc", "xl_hdl_l", "xl_hdl_p", "xl_hdl_pl", "xl_hdl_tg",
    "l_hdl_ce", "l_hdl_fc", "l_hdl_l", "l_hdl_p", "l_hdl_pl",
    "m_hdl_c", "m_hdl_ce", "m_hdl_fc","m_hdl_l", "m_hdl_p", "m_hdl_pl",
    "s_hdl_l", "s_hdl_p", "s_hdl_tg")
df$Var2 <- factor(df$Var2, rev(traits))

library(ggplot2)
p1 <- ggplot(df) + aes(x = Var1, y = Var2, fill = pmax(value, -12)) + geom_tile() + scale_fill_gradient2(low="red", high="white", midpoint = -1) + theme_minimal() + theme(axis.text.x = element_text(angle = 90), legend.title = element_blank(), legend.position = "top") + xlab("") + ylab("")

dfq <- data.frame(cbind(beta_q5, beta_q50, beta_q95))
dfq$SNP <- factor(rownames(dfq), ordered.SNPs)

p2 <- ggplot(dfq) + aes(x = SNP, y = beta_q50, ymin = beta_q5, ymax = beta_q95) + geom_point() + geom_errorbar() + xlab("") + ylab("") + theme_minimal() + theme(axis.text.x = element_blank())

library(cowplot)
plot_grid(p1, p2, ncol = 1, align = "v", rel_heights = c(4, 1))
