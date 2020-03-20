library(phyloseq)
library(ggplot2)
library(ggdendro)
library(cowplot)
library(phyloseq)
library(ggplot2)
library(vegan)
library(extrafont)

set.seed(007)

# set global colors
set.color <- as.character(c("#4DAF4A", "#984EA3",
                            "#E41A1C", "#377EB8",
                            "#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3",
                            "#E69F00", "#009E73", "#0072B2", "#D55E00", "#CC79A7"))
names(set.color) <- c("C", "V",
                      "FALSE", "TRUE",
                      "p__Actinobacteria", "p__Bacteroidetes", "p__Firmicutes", "p__Fusobacteria", "p__Proteobacteria",
                      "Birth", "1week", "1month", "3month", "6month")

minfilt <- function(psobj, prevpercent){
  # remove taxa that didn't have phylum level classification and mitochondria
  ps0 <- subset_taxa(psobj, !is.na(Phylum) & Family != "f__mitochondria")
  # place holder to calculate final sample count for prevalence threshold
  nsam <- sum(sample_sums(ps0) > 5000 & sample_data(ps0)$Timepoint != "Birth")
  # remove files without read counts
  ps1 <- prune_samples(sample_sums(ps0) > 5000, ps0)
  # agglomerate to genus level
  ps2 <- tax_glom(ps1, "Genus", NArm = FALSE)
  # prevalence cutoff
  tmp <- apply(otu_table(ps2), 2, function(x){sum(x > 0)})
  ps3 <- prune_taxa(tmp >= 0.05 * nsam, ps2)
  return(ps3)
}

ps0 <- readRDS("../16S analysis/phyloseq_GG_tree_nodup.RData")
ps <- minfilt(ps0, 0.05)
mtdata <- data.frame(sample_data(ps))

###############################################
#
# FIG 1A. MINIMUM FILTERED BARPLOT W DENDROGRAM ON TOP - WUNIFRAC DISTANCE
#
###############################################

phy <- tax_glom(ps, "Phylum")

# transform to relative abundance
phy.ra <- transform_sample_counts(phy, function(x) x / sum(x) * 100)

dendlist <- list()
cluslist <- list()

for (i in levels(mtdata$Timepoint)){
  psdf <- prune_samples(sample_data(phy.ra)$Timepoint == i, phy.ra)
  wunifrac <- distance(psdf, "wunifrac")
  row.clus <- hclust(wunifrac, method = "ward.D")
  dend <- as.dendrogram(row.clus)
  ddata <- dendro_data(dend, type = "rectangle")
  ddata <- data.frame(segment(ddata), Timepoint = i)
  if(i == "Birth"){
    ddata$x <- ddata$x * 1.05
    ddata$xend <- ddata$xend * 1.05
  }
  dendlist[[i]] <- ddata
  cluslist[[i]] <- row.clus$labels[row.clus$order]
}

denddf <- do.call(rbind, dendlist)

t.labs <- c("Birth", "1 week", "1 month", "3 months", "6 months")
names(t.labs) <- c("Birth", "1week", "1month", "3month", "6month")

p1 <- ggplot(denddf, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_segment(size = 0.2) +
  facet_grid(~Timepoint, scales = "free", space = "free", labeller = labeller(Timepoint = t.labs)) +
  scale_x_continuous(expand = c(0.025, 0.025)) +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.margin = unit(c(2, 10, 0, 10), "pt"),
        panel.spacing.x = unit(0.4, "lines"),
        strip.text = element_text(face = "italic", size = 18.75, margin = margin(0, 0, 0, 0, "pt")))

# make df to put data in
long <- psmelt(phy.ra)
long$Sample <- factor(long$Sample, levels = unlist(cluslist))

# make bottom plot
p2 <- ggplot(long, aes(x = Sample, y = Abundance)) +
  geom_bar(stat = "identity", position = "stack", color = "lightgrey", aes(fill = Phylum), size = 1e-6) +
  geom_point(data = mtdata, aes(x = file.names, y = 105, color = Supplemented.1mo)) +
  geom_point(data = mtdata, aes(x = file.names, y = 102, color = Delivery)) +
  scale_color_manual(name = "Sample Group",
                     breaks = c("C", "V", "FALSE", "TRUE"),
                     values = set.color,
                     labels = c("Caesarian", "Vaginal", "Unsupp.", "Supp.")) +
  scale_fill_manual(labels = c("Actinobacteria", "Bacteroidetes", "Firmicutes", "Fusobacteria", "Proteobacteria"),
                    values = set.color) +
  labs(x = "Relative Abundance") +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  guides(color = guide_legend(order = 1),
         fill = guide_legend(order = 2)) +
  theme(axis.title = element_text(size = 17.5),
        axis.text = element_text(size = 6),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.margin = unit(c(-13, 10, 0, 10), "pt"),
        panel.spacing.x  = unit(0.2, "lines"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, -8)) +
  facet_grid(~Timepoint, scales = "free", space = "free")

###############################################
#
# FIG 1B. MINIMUM FILTERED DIVERSITY
#
###############################################

bac <- data.frame(otu_table(ps))

diversity <- diversity(bac, index = "shannon")

# calculate p-values
pval <- data.frame(matrix(0, ncol = 2, nrow = 5))
colnames(pval) <- c("Timepoint", "P.value")
pval$Timepoint <- factor(levels(mtdata$Timepoint), levels = levels(mtdata$Timepoint))
for (i in 1:nrow(pval)){
  pval$P.value[i] <- round(t.test(diversity[mtdata$Supplemented.1mo == "TRUE" & mtdata$Timepoint == pval$Timepoint[i]],
                            diversity[mtdata$Supplemented.1mo == "FALSE" & mtdata$Timepoint == pval$Timepoint[i]])$p.value, 3)
}
pval$label <- ifelse(pval$P.value < 0.05, "sig", "ns")

fig1b <- ggplot(mtdata, aes(x = Supplemented.1mo, y = diversity)) +
  geom_boxplot(aes(fill = Supplemented.1mo), outlier.shape = NA, lwd = 0.2) +
  geom_jitter(width = 0.15) +
  scale_fill_manual(values = set.color) +
  scale_x_discrete(labels = c("Unsupp.", "Supp.")) +
  labs(x = "Supplemented by age one month",
       y = "Shannon Diversity") +
  theme(axis.title = element_text(size = 17.5),
        axis.text.x = element_text(size = 12, color = c("#E41A1C", "#377EB8")),
        axis.text.y = element_text(size = 6),
        strip.text = element_text(size = 15, margin = margin(0, 0, 0, 0, "pt"), face = "italic"),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black"),
        plot.margin = unit(c(0, 10, 2, 10), "pt")) +
  facet_grid(~Timepoint, labeller = labeller(Timepoint = t.labs)) +
  geom_text(data = pval, aes(label = label, x = 1.5, y = 2.4), size = 0.1)

#ggsave("diversity_final.pdf", plot = fig1b, width = 11, height = 8.5, units = "in", device = "pdf")

###############################################
#
# FIG 1CDE. MINIMUM FILTERED WEIGHTED UNIFRAC PCOA
#
###############################################

# convert to relative abundance for calculating bray curtis distance
ps.ra <- transform_sample_counts(ps, function(x){x / sum(x)})

# calculate PCoA coordinates using weigthed unifrac distance
ord.pcoa <- ordinate(ps.ra, method = "PCoA", distance = "wunifrac")
eig <- round(ord.pcoa$values$Relative_eig * 100, 1)

# setup df to do stat on
wunifrac <- phyloseq::distance(ps.ra, method = "wunifrac")

#pcoa plot based on different criteria
criteria <- c("Supplemented.1mo", "Timepoint", "Delivery")
pcoalist <- list()

for (i in criteria){
  # do some stat to label the plots with
  a <- adonis(as.formula(paste("wunifrac ~ ", i, sep = "")), data = mtdata)
  b <- permutest(betadisper(wunifrac, unlist(mtdata[i])))
  
  dat <- data.frame(Axis.1 = ord.pcoa$vectors[,1],
                    Axis.2 = ord.pcoa$vectors[,2],
                    crit = mtdata[,i])
  
  pcoalist[[i]] <- ggplot(dat, aes(x = Axis.1, y = Axis.2, color = crit)) +
    geom_point(size = 3, alpha = 0.8) +
    theme(axis.title = element_text(size = 17.5),
          axis.text = element_text(size = 6),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = "NA"),
          aspect.ratio = 1,
          plot.margin = unit(c(2, 0, 2, 0), "pt"),
          legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(0, 0, 0, -8)) +
    xlab(paste("Axis 1 (", eig[1], "%)", sep = "")) +
    ylab(paste("Axis 2 (", eig[2], "%)", sep = "")) +
    labs(color = i) +
    annotate("text",  x = -Inf, y = Inf,
             vjust = 1.1, hjust = -0.1,
             label = paste("adonis: p = ", round(a$aov.tab$`Pr(>F)`[1],3),
                           "\nbetadisper: p = ", round(b$tab$`Pr(>F)`[1], 3), sep = ""))
}

# minor adjustments on each pcoalist graph
pcoalist[[1]] <- pcoalist[[1]] +
  scale_color_manual(name = "Supp. by 1mo", values = set.color, labels = c("Unsupp.", "Supp."))

pcoalist[[2]] <- pcoalist[[2]] +
  scale_color_manual(name = "Time point", values = set.color, labels = c("Birth", "1 week", "1 month", "3 months", "6 months"))

pcoalist[[3]] <- pcoalist[[3]] +
  scale_color_manual(values = set.color, labels = c("Caesarian", "Vaginal"))

# put all plots together
fig1ab <- plot_grid(p1, p2, fig1b, ncol = 1, labels = c("A", "", "B"),
                    label_size = 22.5, rel_heights = c(14, 56, 35), align = "v", axis = "rl",
                    vjust = 1, hjust = 0)

fig1cde <- plot_grid(pcoalist[[1]],
                     pcoalist[[2]],
                     pcoalist[[3]],
                     ncol = 1, align = "v",
                     labels = c("C", "D", "E"), label_size = 22.5, vjust = 1, hjust = 0.5)

plot_grid(fig1ab, fig1cde, nrow = 1, rel_widths = c(7, 3))

#ggsave("fig1.pdf", width = 13, height = 7, units = "in", device = "pdf")
