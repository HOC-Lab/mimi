library(phyloseq)
library(ggplot2)
library(edgeR)
library(tidyr)
library(cowplot)
library(gdata)
library(qvalue)

# set global colors
set.color <- as.character(c("#4DAF4A", "#984EA3",
                            "#E41A1C", "#377EB8",
                            "#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3",
                            "#E69F00", "#009E73", "#0072B2", "#D55E00", "#CC79A7"))
names(set.color) <- c("C", "V",
                      "FALSE", "TRUE",
                      "p__Actinobacteria", "p__Bacteroidetes", "p__Firmicutes", "p__Fusobacteria", "p__Proteobacteria",
                      "Birth", "1week", "1month", "3month", "6month")

t.labs <- c("Birth", "1 week", "1 month", "3 months", "6 months")
names(t.labs) <- c("Birth", "1week", "1month", "3month", "6month")

st.labs <- c("Birth", "1 wk", "1 mo", "3 mo", "6 mo")
names(st.labs) <- c("Birth", "1week", "1month", "3month", "6month")

filterPS <- function(psobj, prevpercent) {
  ps0 <- subset_taxa(psobj, !is.na(Phylum) & Family != "f__mitochondria")
  ps1 <- prune_samples(sample_sums(ps0) > 5000 & sample_data(ps0)$Timepoint != "Birth", ps0)
  ps2 <- tax_glom(ps1, "Genus", NArm = FALSE)
  tmp <- apply(otu_table(ps2), 2, function(x){sum(x > 0)})
  ps3 <- prune_taxa(tmp > prevpercent * nsamples(ps2), ps2)
  ps3.ra <- transform_sample_counts(ps3, function(x){x / sum(x)})
  tmp <- apply(otu_table(ps3.ra), 1, function(x){max(x) < 0.9})
  ps4 <- prune_samples(tmp, ps3)
  # Add Highest_Taxonomy column to tax_table
  taxonomy <- data.frame(as(tax_table(ps4),"matrix"), stringsAsFactors = FALSE)
  taxonomy$Species <- NULL
  taxonomy$Highest_Taxonomy <- ifelse(!is.na(taxonomy$Genus), taxonomy$Genus,
                                      ifelse(!is.na(taxonomy$Family), taxonomy$Family,
                                             ifelse(!is.na(taxonomy$Order), taxonomy$Order,
                                                    ifelse(!is.na(taxonomy$Class), taxonomy$Class, taxonomy$Phylum))))
  taxonomy$Highest_Taxonomy[taxonomy$Highest_Taxonomy == "g__Clostridium"] <- paste(taxonomy$Family[taxonomy$Highest_Taxonomy == "g__Clostridium"],
                                                                                    taxonomy$Genus[taxonomy$Highest_Taxonomy == "g__Clostridium"],
                                                                                    sep = ";")
  tax_table(ps4) <- tax_table(as.matrix(taxonomy))
  return(ps4)
}

ps <- readRDS("../16S analysis/phyloseq_GG_tree_nodup.RData")
filtps <- filterPS(ps, 0.05)

m <- t(as(otu_table(filtps), "matrix"))
taxonomy <- tax_table(filtps)
mtdata <- data.frame(sample_data(filtps))
mtdata$group <- as.factor(paste(mtdata$Delivery, mtdata$Timepoint, sep = ""))

d <- DGEList(counts = m, genes = taxonomy)

#########################################
#
# BETWEEN SUP AND UNSUP, ACCOUNTING FOR INTERACTION W AGE & DELIVERY
#
#########################################
#mtdata$age.month <- factor(mtdata$age.month)

mm <- model.matrix(~0 + Delivery*age.month, data = mtdata)
y <- voom(d, mm, plot = TRUE)
cor <- duplicateCorrelation(y, mm, block = mtdata$Participant.ID)$consensus
fit <- lmFit(y, mm, block = mtdata$Participant.ID, correlation = cor)
colnames(fit$coefficients) <- make.names(colnames(fit$coefficients))

contr <- makeContrasts(DeliveryV - DeliveryC, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
res1 <- topTable(tmp, sort.by = "P", n = Inf)

# wtfall plot
sigdf <- res1[res1$P.Value < 0.05,]
sigdf$Highest_Taxonomy <- factor(sigdf$Highest_Taxonomy, levels = sigdf$Highest_Taxonomy[order(sigdf$logFC, decreasing = TRUE)])
sigdf$sig <- ifelse(sigdf$adj.P.Val < 0.05, "*", "")
fig2a <- ggplot(sigdf, aes(x = Highest_Taxonomy, y = logFC, fill = Phylum)) +
  geom_col(position = "dodge") +
  geom_text(label = sigdf$sig, size = 10) +
  scale_y_continuous(expand = c(0.1, 0.05)) +
  scale_fill_manual(values = set.color, labels = c("Actinobacteria", "Bacteroidetes", "Firmicutes", "Proteobacteria")) +
  labs(y = "log2 FC \n (Vaginal - Caesarean)") +
  theme(axis.title = element_text(size = 17.5),
        axis.text.x = element_text(size = 10, angle = 60, hjust = 1),
        axis.text.y = element_text(size = 6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 0.5),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0, 10, 2, 10), "pt"))

# individual plots
pretty_labs <- as.character(c("Bifidobacterium", "Bacteroides", "Parabacteroides",
                              "Caesarean", "Vaginal"))
names(pretty_labs) <- c("g__Bifidobacterium", "g__Bacteroides", "g__Parabacteroides",
                        "C", "V")

#ps.ra <- transform_sample_counts(filtps, function(x) x / sum(x) * 100)
ps.ra <- transform_sample_counts(filtps, function(x) log2(x + 1))
long <- psmelt(ps.ra)
long <- long[long$Highest_Taxonomy %in% res1$Highest_Taxonomy[res1$adj.P.Val < 0.05],]
long$Delivery <- factor(long$Delivery, levels = c("V", "C"))
mean <- aggregate(long$Abundance, by = list(long$Timepoint, long$Delivery, long$Highest_Taxonomy), mean)
colnames(mean) <- c("Timepoint", "Delivery", "Highest_Taxonomy", "Abundance")

fig2b <- ggplot(long, aes(x = Timepoint, y = Abundance, color = Delivery)) +
  geom_line(aes(group = Participant.ID), linetype = 2, size = 0.5) +
  geom_line(data = mean, aes(group = 1), size = 1) +
  facet_grid(Highest_Taxonomy ~ Delivery, scale = "free_y",
             labeller = labeller(Highest_Taxonomy = pretty_labs, Delivery = pretty_labs)) +
  scale_color_manual(values = set.color) +
  scale_x_discrete(expand = c(0.05,0.05), labels = c("1wk", "1mo", "3mo", "6mo")) +
  ylab("log2 (count + 1)") +
  theme(axis.title = element_text(size = 17.5),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = "NA"),
        legend.position = "none",
        strip.text = element_text(size = 15, margin = margin(0, 0, 0, 0, "pt"), face = "italic"))

#########################################
#
# LIMMA VOOM - REPLICATED TIME POINTS AS FACTOR
#
#########################################

mm <- model.matrix(~0 + group, data = mtdata)
y <- voom(d, mm, plot = TRUE)
cor <- duplicateCorrelation(y, mm, block = mtdata$Participant.ID)$consensus
fit <- lmFit(y, mm, block = mtdata$Participant.ID, correlation = cor)

# dif between ages for unsupplemented infants
contr <- makeContrasts("1week" = "groupV1week - groupC1week",
                       "1month" = "groupV1month - groupC1month",
                       "3month" = "groupV3month - groupC3month",
                       "6month" = "groupV6month - groupC6month",
                       levels = coef(fit))
colnames(coef(fit))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)

res2 <- list()
for (i in 1:4){
  res2[[i]] <- topTable(tmp, coef = i, sort.by = "P", n = Inf)
  res2[[i]]$Timepoint <- colnames(tmp$coefficients)[i]
}

res2 <- do.call(rbind, res2)
res2$Timepoint <- factor(res2$Timepoint, levels = c("1week", "1month", "3month", "6month"))
res2$logp <- -log10(res2$P.Value)
res2$adj.P.value <- ifelse(res2$adj.P.Val < 0.05, "< 0.05", "> 0.05")
plabel <- data.frame(Timepoint = c("1week", "1month" , "3month", "6month"),
                     label = c(" Un-adjusted p = 0.05", "", "", ""))
res2$siglab <- ""
res2$siglab[res2$adj.P.Val < 0.05] <- res2$Highest_Taxonomy[res2$adj.P.Val < 0.05]

fig2c <- ggplot(res2, aes(x = logFC, y = logp)) +
  geom_point(size = 2, alpha = 0.8, aes(color = Phylum, shape = adj.P.value)) +
  geom_point(data = res2[res2$adj.P.Val < 0.05,], shape = 21, stroke = 1.5, size = 2, aes(fill = Phylum), show.legend = FALSE) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  facet_grid(~Timepoint, labeller = labeller(Timepoint = t.labs)) +
  scale_color_manual(values = set.color,
                     labels = c("Actinobacteria", "Bacteroidetes", "Firmicutes", "Fusobacteria", "Proteobacteria")) +
  scale_shape_manual(values = c(21, 16), guide = FALSE) +
  labs(x = "log2 Fold Change (Vaginal - Caesarean)", y = "-log10(P.Value)") +
  geom_text(data = plabel, aes(label = plabel$label, x = -Inf, y = -log10(0.05)*1.1), hjust = 0, size = 2) +
  geom_text(aes(label = siglab), hjust = 1.2, size = 2) +
  #annotate("text", label = " Un-adjusted P = 0.05", x = -Inf, y = -log10(0.05)*1.1, hjust = 0, size = 3) +
  theme(axis.title = element_text(size = 17.5),
        axis.text = element_text(size = 6),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        aspect.ratio = 1,
        legend.position = "top",
        legend.box = "vertical",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = "NA"),
        strip.text = element_text(face = "italic", size = 18.75, margin = margin(0, 0, 0, 0, "pt")),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, -8),
        plot.margin = unit(c(0, 10, 2, 10), "pt"))

# put all plots together

fig2ac <- plot_grid(fig2a, fig2c, ncol = 1, labels = c("A", "C"), label_size = 22.5, rel_heights = c(1.2, 1))

plot_grid(fig2ac, fig2b, nrow = 1, labels = c("", "B"), label_size = 22.5, rel_widths = c(1.5,1))

#ggsave("fig2.pdf", width = 11, height = 6.5, units = "in", device = "pdf")

