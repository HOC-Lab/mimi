library(qvalue)
library(ade4)
library(ggplot2)
library(gdata)
library(vegan)
library(lmerTest)
library(cowplot)
library(tidyr)
library(phyloseq)

set.seed(070)

# set global colors
set.color <- as.character(c("#4DAF4A", "#984EA3",
                            "#E41A1C", "#377EB8",
                            "#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3",
                            "#E69F00", "#009E73", "#0072B2", "#D55E00", "#CC79A7"))
names(set.color) <- c("C", "V",
                      "FALSE", "TRUE",
                      "p__Actinobacteria", "p__Bacteroidetes", "p__Firmicutes", "p__Fusobacteria", "p__Proteobacteria",
                      "Birth", "1week", "1month", "3month", "6month")

####################################
#
# LOAD DATA
#
####################################

imm <- read.csv("../PBMC analysis/MIMI cell phenotype.csv")
# clean up mtdata
mtdata <- read.xls("../MIMIStudy_rawfromValerie.xlsx", sheet = "for R_updated190710")
# reorder mtdata to match imm df
imtdata <- merge(data.frame(Participant.ID = as.factor(imm$Infant.ID),
                           Collection.date = as.Date(imm$Date, format = "%m/%d/%Y"),
                           Timepoint = factor(gsub(" ", "", imm$Timepoint), levels = c("Birth", "1month", "3month", "6month")),
                           age.month = ifelse(imm$Timepoint == "Birth", 0,
                                              ifelse(imm$Timepoint == "1 month", 1,
                                                     ifelse(imm$Timepoint == "3 month", 3, 6)))),
                mtdata, by = "Participant.ID")

# clean up mtdata
imtdata$Birthday <- as.Date(imtdata$Birthday)
imtdata$Influenza.date <- as.Date(imtdata$Influenza.date)
rownames(imtdata) <- paste(imtdata$Participant.ID, imtdata$Timepoint, sep = ".")

# clean up imm data
rownames(imm) <- rownames(imtdata)
imm[,1:4] <- NULL

# add total monocyte column
imm$Total.Monocyte <- imm$Monocytes.CD16nCD14p + imm$Monocytes.CD16pCD14p + imm$Monocytes.CD16pCD14n
imm$Total.Monocyte.CD80p <- imm$Monocytes.CD16nCD14p.CD80p + imm$Monocytes.CD16pCD14p.CD80p + imm$Monocytes.CD16pCD14n.CD80p
imm$Monocytes.CD16nCD14n <- NULL
imm$Monocytes.CD16nCD14n.CD80p <- NULL

####################################
#
# PCA
#
####################################

pca <- dudi.pca(imm, nf = 5, scannf = FALSE)
axis <- pca$li
loading <- pca$c1

# dif between sup
eu <- dist(imm)
a1 <- adonis(eu ~ Supplemented.1mo, data = imtdata, method = "eu")
b1 <- permutest(betadisper(eu, imtdata$Supplemented.1mo))

##using data table from ade4
fig3a <- ggplot(axis, aes(x = Axis1, y = Axis2, col = imtdata$Supplemented.1mo))+
  geom_point(size = 2, alpha = 0.8)+
  #stat_ellipse() +
  xlab(paste("PC1 (", round(pca$eig/sum(pca$eig)*100, 1)[1], "%)", sep = "")) +
  ylab(paste("PC2 (", round(pca$eig/sum(pca$eig)*100, 1)[2], "%)", sep = "")) +
  theme(axis.title = element_text(size = 17.5),
        axis.text = element_text(size = 6),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = "NA"),
        aspect.ratio = 1,
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, -8),
        plot.margin = unit(c(2, 2, 2, 10), "pt")) +
  annotate("text",  x = -Inf, y = Inf,
           vjust = 1.1, hjust = -0.1, size = 2,
           label = paste(" adonis: p = ", round(a1$aov.tab$`Pr(>F)`[1],3),
                         "\nbetadisper: p = ", round(b1$tab$`Pr(>F)`[1], 3), sep = "")) +
  scale_color_manual(name = "Supp. by 1mo", values = set.color, labels = c("Unsupp.", "Supp."))

# dif between timepoints
a2 <- adonis(eu ~ Timepoint, data = imtdata, method = "eu")
b2 <- permutest(betadisper(eu, imtdata$Timepoint))

##using data table from ade4
fig3b <- ggplot(axis, aes(x = Axis1, y = Axis2, col = imtdata$Timepoint))+
  geom_point(size = 2, alpha = 0.8)+
  #stat_ellipse() +
  xlab(paste("PC1 (", round(pca$eig/sum(pca$eig)*100, 1)[1], "%)", sep = "")) +
  ylab(paste("PC2 (", round(pca$eig/sum(pca$eig)*100, 1)[2], "%)", sep = "")) +
  theme(axis.title = element_text(size = 17.5),
        axis.text = element_text(size = 6),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = "NA"),
        aspect.ratio = 1,
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, -8),
        plot.margin = unit(c(0, 2, 2, 10), "pt")) +
  annotate("text",  x = -Inf, y = Inf,
           vjust = 1.1, hjust = -0.1, size = 2,
           label = paste(" adonis: p = ", round(a2$aov.tab$`Pr(>F)`[1],3),
                         "\nbetadisper: p = ", round(b2$tab$`Pr(>F)`[1], 3), sep = "")) +
  scale_color_manual(name = "Time point", values = set.color, labels = c("Birth", "1 week", "1 month", "3 months", "6 months"))

####################################
#
# MIXED MODEL - SUP 1 MO
#
####################################

out <- data.frame(matrix(0, ncol = 4, nrow = ncol(imm)))
colnames(out) <- c("imm", "SupTRUE", "age.month", "SupTRUE_age.month")
out$imm <- colnames(imm)

for (i in 1:nrow(out)){
  mod <- anova(lmer(imm[,i] ~ imtdata$Supplemented.1mo * imtdata$age.month + (1|imtdata$Participant.ID)))
  out[i, 2:4] <- mod[1:3, 6]
}

# int term: CD4.HLADRp is significant - but seemed to be driven by one infant
# qvalue for age.month p val
out$age.month.q <- qvalue(out$age.month)$qvalues
length(which(out$age.month < 0.05))
length(which(out$age.month.q < 0.05))
length(which(out$age.month.q < 0.03))

#sig <- out$imm[which(out$age.month.q < 0.05)]
sig <- out$imm[which(out$age.month.q < 0.03)]

long <- cbind(imtdata, imm[,sig])
long <- gather(long, "imm", "percentage", sig)
#mean <- aggregate(long$percentage, by = list(long$age.month, long$Supplemented.1mo, long$imm), mean)
#colnames(mean) <- c("age.month", "Supplemented.1mo", "imm", "percentage")

pretty_labs <- c("Unsupplemented", "Supplemented",
                 "CD4+HLA-DR+ \nT cell", "CD4+HLA-DR+ \nem T cell", "CD4+HLA-DR+ \nmem T cell", "CD8+ \nmem T cell", "CD16+CD14- \nmonocytes")
names(pretty_labs) <- c("FALSE", "TRUE",
                        "CD4.HLADRp", "emCD4.HLADRp", "memCD4.HLADRp", "memCD8", "Monocytes.CD16pCD14n")

#ggplot(imtdata, aes(x = age.month, y = imm$CD8)) +
#  geom_smooth(method = "lm", se = FALSE, size = 1) +
#  geom_line(aes(group = Participant.ID), linetype = 2, size = 0.5) +
#  facet_wrap(~Supplemented.1mo)

fig3c <- ggplot(long, aes(x = age.month, y = percentage, color = Supplemented.1mo)) +
  geom_line(aes(group = Participant.ID), linetype = 2, size = 0.5) +
  #geom_line(data = mean, aes(group = 1), size = 1) +
  geom_smooth(method = "lm", se = FALSE, size = 1) +
  labs(x = "Age (month)", y = "Population Percentage") +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  theme(axis.title = element_text(size = 17.5),
        axis.text = element_text(size = 6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = "NA"),
        legend.position = "none",
        strip.text = element_text(size = 10, margin = margin(0, 0, 0, 0, "pt"), face = "italic"),
        plot.margin = unit(c(2, 10, 2, 10), "pt")) +
  scale_color_manual(values = set.color) +
  facet_grid(imm~Supplemented.1mo, scale = "free", labeller = labeller(Supplemented.1mo = pretty_labs, imm = pretty_labs))

####################################
#
# VACCINE RESPONSE
#
####################################

igg <- read.xls("../ELISA/IgG ELISA.xlsx", sheet = "TT IgG", stringsAsFactors = FALSE)
igm <- read.xls("../ELISA/IgM ELISA.xlsx", sheet = "TT IgM", stringsAsFactors = FALSE)

# merge IgM and IgG data, make sure the order is still the same
TT <- dplyr::right_join(igm[,c("Sample", "TT.IgM..IU.mL.")], igg, by = "Sample")

dat <- merge(data.frame(Participant.ID = as.factor(TT$Participant.ID),
                        Collection.date = as.Date(as.character(TT$Collection.date)),
                        Timepoint = factor(gsub(" ", "", TT$Timepoint), levels = c("Birth", "1month", "3month", "6month")),
                        age.month = ifelse(TT$Timepoint == "Birth", 0,
                                           ifelse(TT$Timepoint == "1 month", 1,
                                                  ifelse(TT$Timepoint == "3 month", 3, 6)))),
             mtdata, by = "Participant.ID")

# clean up mtdata
dat$Birthday <- as.Date(dat$Birthday)
dat$Influenza.date <- as.Date(dat$Influenza.date)
rownames(dat) <- paste(dat$Participant.ID, dat$Timepoint, sep = ".")
dat$igm <- TT$TT.IgM..IU.mL.
dat$igg <- TT$TT.IgG..IU.mL.

# remove the one infant that didn't get vaccine at 2 months
dat <- dat[dat$MultiVaccine.2month == TRUE,]

fig3d <- ggplot(dat, aes(x = age.month, y = igg, color = Supplemented.1mo)) +
  geom_rect(aes(xmin = 2, xmax = Inf, ymin = -Inf, ymax = Inf), fill = "goldenrod1", color = NA, alpha = 0.02) +
  geom_path(aes(group = Participant.ID), linetype = 2, size = 0.5) +
  geom_smooth(method = "lm", se = FALSE, size = 1) +
  labs(x = "Age (month)", y = "Tetanus Toxoid \nIgG (IU/mL)") +
  theme(axis.title = element_text(size = 17.5),
        axis.text = element_text(size = 6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = "NA"),
        legend.position = "none",
        #aspect.ratio = 1,
        strip.text = element_text(size = 15, margin = margin(0, 0, 0, 0, "pt"), face = "italic"),
        plot.margin = unit(c(0, 10, 2, 10), "pt")) +
  scale_color_manual(values = set.color) +
  facet_wrap(~Supplemented.1mo, labeller = labeller(Supplemented.1mo = pretty_labs)) +
  annotate("text", label = " Post-vaccination", x = 2, y = -Inf, hjust = 0, vjust = -0.5, size = 2)

fig3e <- ggplot(dat, aes(x = age.month, y = igm, color = Supplemented.1mo)) +
  geom_rect(aes(xmin = 2, xmax = Inf, ymin = -Inf, ymax = Inf), fill = "goldenrod1", color = NA, alpha = 0.02) +
  geom_path(aes(group = Participant.ID), linetype = 2, size = 0.5) +
  geom_smooth(method = "lm", se = FALSE, size = 1) +
  labs(x = "Age (month)", y = "Tetanus Toxoid \nIgM (IU/mL)") +
  theme(axis.title = element_text(size = 17.5),
        axis.text = element_text(size = 6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = "NA"),
        legend.position = "none",
        #aspect.ratio = 1,
        strip.text = element_text(size = 15, margin = margin(0, 0, 0, 0, "pt"), face = "italic"),
        plot.margin = unit(c(0, 10, 2, 10), "pt")) +
  scale_color_manual(values = set.color) +
  facet_wrap(~Supplemented.1mo, labeller = labeller(Supplemented.1mo = pretty_labs)) +
  annotate("text", label = " Post-vaccination", x = 2, y = -Inf, hjust = 0, vjust = -0.5, size = 2)

####################################
#
# VACCINE RESPONSE COR W BIF
#
####################################

#filterPSg <- function(psobj, prevpercent) {
#  ps0 <- subset_taxa(psobj, !is.na(Phylum))
#  tmp <- sample_sums(ps0) > 5000
#  ps1 <- prune_samples(tmp, ps0)
#  ps2 <- tax_glom(ps1, "Genus", NArm = TRUE)
#  # maybe get some changes in between
#  prevdf <- apply(otu_table(ps2), 2, function(x){sum(x > 0)})
#  prevalenceThreshold <-  prevpercent * nsamples(ps2)
#  keepTaxa <- names(prevdf)[(prevdf >= prevalenceThreshold)]
#  ps3 <- prune_taxa(keepTaxa, ps2)
#  ps3.ra <- transform_sample_counts(ps3, function(x){x / sum(x)})
#  tmp <- apply(otu_table(ps3.ra), 1, function (x) {max(x)}) < 0.9
#  ps4 <- prune_samples(tmp, ps3)
#  return(ps4)
#}

#ps <- readRDS("../16S analysis/phyloseq_GG_tree_nodup.RData")
#filtpsg <- filterPSg(ps, 0.05)
#filtpsg.ra <- transform_sample_counts(filtpsg, function(x){x / sum(x)})
#bif <- subset_taxa(filtpsg.ra, Genus == "g__Bifidobacterium")
#bif <- data.frame(otu_table(bif))

#all <- cbind(dat[intersect(rownames(dat), rownames(bif)),], bif[intersect(rownames(dat), rownames(bif)),])
#colnames(all) <- c(colnames(dat), "bif_genus")

# keep 3 month and 6 month samples only
#all2 <- all[all$age.month %in% c(3, 6),]
#ttcor <- cor.test(all2$TT, all2$bif_genus, method = "spearman")
#ttcor$p.value
#ttcor$estimate

#fig3e <- ggplot(all[all$age.month %in% c(3, 6),], aes(x = bif_genus, y = TT)) +
#  geom_point(size = 3) +
#  labs(x = "Relative Abundance of Bifidobacterium Genus (%)", y = "Tetanus Toxoid IgG (IU/mL)") +
#  theme(text = element_text(size = 15),
#        panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank(),
#        panel.background = element_blank(),
#        panel.border = element_rect(fill = "NA"),
#        aspect.ratio = 1) +
#  annotate("text", x = -Inf, y = Inf,
#           vjust = 1.1, hjust = 0,
#           label = paste0(" rho = ", round(ttcor$estimate, 3),
#                         "\n p-value = ", round(ttcor$p.value, 3)))

####################################
#
# PUT EVERYTHING TOGETHER
#
####################################

fig3ab <- plot_grid(fig3a, fig3b, ncol = 1, align = "v", labels = c("A", "B"), label_size = 22.5)
fig3abc <- plot_grid(fig3ab, fig3c, nrow = 1, labels = c("", "C"), label_size = 22.5)
fig3de <- plot_grid(fig3d, fig3e, nrow = 1, align = "h", axis = "tb", labels = c("D", "E"), label_size = 22.5, vjust = 0)

plot_grid(fig3abc, fig3de, ncol = 1, rel_heights = c(2.5, 1))

#ggsave("fig3_v3.pdf", width = 8.5, height  = 8, device = "pdf", units = "in")

