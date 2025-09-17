###############################################
##### Guo et al., 2025 - Tannin induction #####
###############################################

library(nortest)
library(sf)
library(ggplot2) 
library(knitr)
library(tidyverse)
library(brms)
library(ape)
library(coda)
library(modelr)
library(gridExtra)
library(pBrackets)
library(RColorBrewer)
library(performance)
library(phytools)
library(kableExtra)
library(tidybayes)
library(formattable)
library(grid)
library(Taxonstand)
library(dplyr)
library(ggplot2)
library(bayesplot)
library(metafor)
library(ggbeeswarm)
library(pander)
library(raster)
library(maps)
library(mapdata)
library(rworldmap)
library(readxl)
library(utils)
library(reshape2)
library(ggpubr)
library(stringr)
library(ggstance)
library(ggridges)
library(sp)
library(geodata)
library(multcompView)
library(multcomp)
library(rotl)
library(psych)
library(meta)
library(devtools)
library(processx)
library(metaAidR)
library(gghalves)
library(ggtree)
library(sf)
library(leaflet)
library(ggspatial)
library(clipr)
library(gstat)
library(knitr)
library(tidyverse)
library(brms)
library(ape)
library(coda)
library(modelr)
library(gridExtra)
library(pBrackets)
library(RColorBrewer)
library(performance)
library(phytools)
library(kableExtra)
library(tidybayes)
library(formattable)
library(grid)
library(Taxonstand)
library(base)
library(bayesplot)
library(metafor)
library(ggbeeswarm)
library(pander)
library(raster)
library(maps)
library(mapdata)
library(rworldmap)
library(readxl)
library(utils)
library(reshape2)
library(ggpubr)
library(stringr)
library(ggstance)
library(ggridges) 
library(sp)
library(raster)
library(geodata)
library(car)
library(ggplot2)
library(dplyr)
library(broom)
library(extrafont)
library(tidyverse)
library(tidyselect)
library(lme4)
library(emmeans)
library(broom)
library(cowplot)
library(reshape2)
library(ggcorrplot)
library(MASS)
library(geodata)
library(ggbeeswarm)
library(ggsignif)
library(emmeans)
library(sf)
library(dplyr)
library(tidyverse)
library(agricolae)
library(agricolae)
library(V.PhyloMaker2)
library(ape)
library(ggtree)
library(dplyr)
library(tidyr)
library(treeio)
library(tidytree)
library(scales)
library(taxize)
library(ggExtra)
library(readxl)
library(metafor)
library(patchwork)
library(metafor)

### Figure - tree ###
setwd("/Users/yaolin/Desktop/My papers/Manuscripts/2025 - Induced tannin/Submission")
data <- read_excel("Raw data - Induced tannin - 0322.xlsx")
write.csv(data, file = "data.csv", row.names = FALSE)
data <- read.csv("data.csv")

num_levels <- n_distinct(data$Plant_species, data$Plant_species)
print(num_levels)

unique_publications <- unique(data$Publication)
print(unique_publications)
count_unique <- length(unique_publications)
print(count_unique)

unique_Plant_species <- unique(data$Plant_species)
print(unique_Plant_species)
count_unique <- length(unique_Plant_species)
print(count_unique)

data <- data[,c("lnRR",  
                "rvar",
                "Family",
                "Plant_species",
                "Plant_type")]
print(unique_species, n = 100)

# Figure - tree #
length(unique(data$'Plant_species'))

unique_species <- data %>%
  distinct(data$'Plant_species')
print(unique_species, n = 100)

data <- data %>%
  group_by(Plant_species, Family) %>%
  mutate(lnRR = mean(lnRR, na.rm = TRUE)) %>%
  slice(1) %>%
  ungroup()

print(data, n = 500)

data_tree <- data[, c("lnRR", "Plant_species", "Family", "Plant_type")] %>%
  mutate(Plant_species = str_replace_all(Plant_species, "﻿", "")) %>%
  separate(Plant_species, into = c("genus", "species"), sep = " ", remove = FALSE) %>%
  filter(species != "spp.")

species_list <- data_tree %>% select(species, genus, Family)
phylo_tree_result <- phylo.maker(sp.list = species_list)
phylo_tree <- phylo_tree_result$scenario.3

tree_data <- phylo_tree
tree_data$tip.label <- gsub("_", " ", tree_data$tip.label)

data_tree <- data_tree %>% rename(label = species)
data_tree <- data_tree[ , !duplicated(names(data_tree))]
unique(tree_data$tip.label)
unique(data_tree$label)

p <- ggtree(tree_data, layout = "circular") + geom_tree()
p$data$label <- tolower(p$data$label)
data_tree$label <- tolower(data_tree$label)

p$data <- p$data %>% left_join(data_tree, by = "label")
p$data <- p$data %>% mutate(species_short = word(label, -1)) 
p$data$species_short <- str_to_title(p$data$species_short)
filtered_data <- p$data %>% filter(species_short %in% tree_data$tip.label)

p <- p %<+% filtered_data
p <- p +
  geom_tippoint(aes(color = Plant_type, size = lnRR)) +
  geom_tiplab(aes(label = Plant_species), size = 3, hjust = -0.1, offset = 1) +
  theme(
    plot.margin = unit(c(3, 3, 3, 3), "cm"),
    text = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.position = c(2.2, 0.2),
    legend.justification = c(2.2, 0.2)
  ) +
  scale_color_manual(values = c(
    "Woody"     = "#f16c23",
    "Non-woody" = "#2b6a99"
  )) +
  scale_size_continuous(range = c(0, 5))
p
ggsave('./Figure 1a.pdf', p, height = 150, width = 300, units = c("mm"))






### Test phylo relationship ###
data_tree <- data_tree %>% rename(species = label)
species_list <- data_tree %>% select(species, genus, Family)
phylo_tree_result <- phylo.maker(sp.list = species_list)
phylo_tree <- phylo_tree_result$scenario.3


setdiff(data$Plant_species, phylo_tree$tip.label)
data$Plant_species <- gsub("\ufeff", "", data$Plant_species)

matched_list <- phylo_tree_result$species.list
genus_species_names <- paste0(matched_list$genus, "_", matched_list$species)
phylo_tree$tip.label <- genus_species_names
data$Plant_species <- gsub(" ", "_", data$Plant_species)
common <- intersect(phylo_tree$tip.label, data$Plant_species)
cat("Matched species count:", length(common), "\n")
print(common)
str(data)
head(data)
sum(is.na(data$Plant_species)) 
sum(is.na(data$lnRR))

if (is.null(phylo_tree$edge.length)) {
  cat("Tree has no branch lengths. Setting them all to 1.\n")
  phylo_tree$edge.length <- rep(1, nrow(phylo_tree$edge))
}
head(phylo_tree$edge.length)

effect_size_species <- data %>%
  group_by(Plant_species) %>%
  summarize(effect_size = mean(lnRR, na.rm = TRUE)) %>%
  as.data.frame()

effect_size_species <- effect_size_species[!is.na(effect_size_species$effect_size), ]
rownames(effect_size_species) <- effect_size_species$Plant_species

effect_size_species$Plant_species <- gsub(" ", "_", effect_size_species$Plant_species)
rownames(effect_size_species)     <- gsub(" ", "_", rownames(effect_size_species))

cat("Tree tip labels example:\n")
print(head(phylo_tree$tip.label, 10))
cat("Data rownames example:\n")
print(head(rownames(effect_size_species), 10))

common_species <- intersect(phylo_tree$tip.label, rownames(effect_size_species))
cat("Number of common species between tree and data:", length(common_species), "\n")

effect_size_species <- effect_size_species[common_species, , drop = FALSE]
effect_size_species <- effect_size_species[match(phylo_tree$tip.label, rownames(effect_size_species)), ]
effect_size_species <- effect_size_species[!is.na(rownames(effect_size_species)), ]
names(effect_size_species$effect_size) <- rownames(effect_size_species)

lambda_result_all <- phylosig(
  tree  = phylo_tree,
  x     = effect_size_species$effect_size,
  method= "lambda",
  test  = TRUE
)
print(lambda_result_all)

# Woody
data_woody <- subset(data, Plant_type == "Woody")
effect_size_species_woody <- data_woody %>%
  group_by(Plant_species) %>%
  summarize(effect_size = mean(lnRR, na.rm = TRUE)) %>%
  as.data.frame()

effect_size_species_woody <- effect_size_species_woody[!is.na(effect_size_species_woody$effect_size), ]
rownames(effect_size_species_woody) <- effect_size_species_woody$Plant_species
effect_size_species_woody$Plant_species <- gsub(" ", "_", effect_size_species_woody$Plant_species)
rownames(effect_size_species_woody)     <- gsub(" ", "_", rownames(effect_size_species_woody))
common_species_woody <- intersect(phylo_tree$tip.label, rownames(effect_size_species_woody))
cat("Overlap species count (Woody):", length(common_species_woody), "\n")
effect_size_species_woody <- effect_size_species_woody[common_species_woody, , drop=FALSE]
extra_species_woody <- setdiff(phylo_tree$tip.label, common_species_woody)
pruned_tree_woody   <- drop.tip(phylo_tree, extra_species_woody)
effect_size_species_woody <- effect_size_species_woody[match(pruned_tree_woody$tip.label,
                                                             rownames(effect_size_species_woody)), , drop=FALSE]
effect_size_species_woody <- effect_size_species_woody[!is.na(rownames(effect_size_species_woody)), ]

x_woody <- effect_size_species_woody$effect_size
names(x_woody) <- rownames(effect_size_species_woody)
length(x_woody)
length(pruned_tree_woody$tip.label)
all.equal(names(x_woody), pruned_tree_woody$tip.label)

lambda_result_woody <- phylosig(
  tree  = pruned_tree_woody,
  x     = x_woody,
  method= "lambda",
  test  = TRUE
)
print(lambda_result_woody)

# Non-woody
data_nonwoody <- subset(data, Plant_type == "Non-woody")
effect_size_species_nonwoody <- data_nonwoody %>%
  dplyr::group_by(Plant_species) %>%
  dplyr::summarize(effect_size = mean(lnRR, na.rm = TRUE)) %>%
  as.data.frame()

effect_size_species_nonwoody <- effect_size_species_nonwoody[!is.na(effect_size_species_nonwoody$effect_size), ]
rownames(effect_size_species_nonwoody) <- effect_size_species_nonwoody$Plant_species
effect_size_species_nonwoody$Plant_species <- gsub(" ", "_", effect_size_species_nonwoody$Plant_species)
rownames(effect_size_species_nonwoody)     <- gsub(" ", "_", rownames(effect_size_species_nonwoody))
common_species_nonwoody <- intersect(phylo_tree$tip.label, rownames(effect_size_species_nonwoody))
cat("Non-woody overlap with tree:", length(common_species_nonwoody), "\n")
effect_size_species_nonwoody <- effect_size_species_nonwoody[common_species_nonwoody, , drop=FALSE]
extra_species_nonwoody <- setdiff(phylo_tree$tip.label, common_species_nonwoody)
pruned_tree_nonwoody   <- drop.tip(phylo_tree, extra_species_nonwoody)
effect_size_species_nonwoody <- effect_size_species_nonwoody[match(pruned_tree_nonwoody$tip.label,
                                                                   rownames(effect_size_species_nonwoody)), ,
                                                             drop=FALSE]
effect_size_species_nonwoody <- effect_size_species_nonwoody[!is.na(rownames(effect_size_species_nonwoody)), ]
x_nonwoody <- effect_size_species_nonwoody$effect_size
names(x_nonwoody) <- rownames(effect_size_species_nonwoody)
cat("Number of tips in pruned_tree_nonwoody:", length(pruned_tree_nonwoody$tip.label), "\n")
cat("Length of x_nonwoody:", length(x_nonwoody), "\n")
all.equal(names(x_nonwoody), pruned_tree_nonwoody$tip.label)
lambda_result_nonwoody <- phylosig(
  tree  = pruned_tree_nonwoody,
  x     = x_nonwoody,
  method= "lambda",
  test  = TRUE
)
print(lambda_result_nonwoody)










### Figure 1b ###
setwd("/Users/yaolin/Desktop/My papers/Manuscripts/2025 - Induced tannin/Submission - 0830")
data <- read_excel("Raw data - Induced tannin - 0322.xlsx")
write.csv(data, file = "data.csv", row.names = FALSE)
data <- read.csv("data.csv")

rm1 <- rma.mv(lnRR, rvar,
              mods = ~ 1,
              random = list( ~ 1 | Publication/Effect_ID, ~ 1 | Plant_species),
              method = "REML",
              data = data)
rm1

rm2 <- rma.mv(lnRR, rvar,
              mods = ~ paste(Plant_type) - 1,
              random = list( ~ 1 | Publication/Effect_ID, ~ 1 | Plant_species),
              method = "REML", 
              data = data)
rm2

coef(rm2)
L <- rbind("Woody - Non-woody" = c(-1, 1))
anova(rm2, L = L) 

sample.sizes <- as.data.frame(table(data$Plant_type))
sample.sizes

dat1 <- data.frame(coef = coef(rm1), se = sqrt(diag(vcov(rm1))), Type = c("Total"))
dat1$lower <- dat1$coef - 1.96 * dat1$se
dat1$upper <- dat1$coef + 1.96 * dat1$se
dat1
dat2 <- data.frame(coef = coef(rm2), se = sqrt(diag(vcov(rm2))), Type = c("Non-woody", "Woody"))
dat2$lower <- dat2$coef - 1.96 * dat2$se
dat2$upper <- dat2$coef + 1.96 * dat2$se
dat2
dat <- rbind(dat1, dat2)
dat

data_total <- data
data_total$Plant_type <- factor("Total", levels = c(levels(data$Plant_type), "Total"))
data_total <- rbind(data, data_total)

p <- ggplot() +
  geom_point(data = data_total, aes(x = lnRR, y = Plant_type, color = Plant_type),
             alpha = 0.5, position = position_jitter(height = 0.15, width = 0), size = 3, shape = 16) +
  geom_errorbarh(data = dat, aes(xmin = lower, xmax = upper, y = Type), 
                 width = 1, size = 1, colour = "grey10", height = 0.06, linewidth = 1.25) +
  geom_point(data = dat, aes(x = coef, y = Type, color = Type),
             shape = 21, fill = "white", color = "black", size = 3, stroke = 1) +
  geom_vline(xintercept = 0, linetype = 2, colour = "steelblue", size = 1) +
  scale_color_manual(values = c("Total" = "#bfbfbf",
                                "Woody" = "#f16c23",
                                "Non-woody" = "#2b6a99")) +
  scale_size_continuous(range = c(2, 6)) +
  theme_classic() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.position = "None") +
  labs(x = expression(paste("Effect size (ln ", italic("RR"), ")")), y = NULL) +
  scale_y_discrete(limits = c("Non-woody", 
                              "Woody", 
                              "Total")) +
  scale_x_continuous(limits = c(-0.6, 0.6), labels = scales::number_format(accuracy = 0.1)) +
  coord_fixed(ratio = 0.6)
p
ggsave('./Figure 1b1.pdf', p, height = 125, width = 200, units = c("mm"))



### Figure 2 ###
# Figure 2a
rm3 <- rma.mv(lnRR, rvar,
              mods = ~ 0 + factor(Tannin_type),
              random = list( ~ 1 | Publication/Effect_ID, ~ 1 | Plant_species),
              method = "REML", 
              data = data)
rm3

coef(rm3)
L <- rbind("Hydrolysable − Condensed" = c(-1,  1,  0),
           "Total − Condensed"        = c(-1,  0,  1),
           "Hydrolysable − Total"     = c( 0,  1, -1))
anova(rm3, L = L)

rm4 <- rma.mv(lnRR, rvar,
              mods = ~ 0 + paste(Tannin_type, Plant_type),
              random = list( ~ 1 | Publication, ~ 1 | Plant_species),
              method = "REML", 
              data = data)
rm4

bn <- names(coef(rm4))
p  <- length(bn)
ix <- function(pat) which(grepl(pat, bn))
mk <- function(plus, minus){
  v <- rep(0, p); v[plus] <- 1; v[minus] <- -1; v
}

i_CW <- ix("Condensed.*Woody")
i_CN <- ix("Condensed.*Non-woody")
i_TW <- ix("Total.*Woody")
i_TN <- ix("Total.*Non-woody")
i_HW <- ix("Hydrolysable.*Woody")
i_HN <- ix("Hydrolysable.*Non-woody")
L <- rbind(if(length(i_CW)==1 && length(i_CN)==1) "Woody − Non-woody | Condensed"    = mk(i_CW, i_CN),
           if(length(i_TW)==1 && length(i_TN)==1) "Woody − Non-woody | Total"        = mk(i_TW, i_TN),
           if(length(i_HW)==1 && length(i_HN)==1) "Woody − Non-woody | Hydrolysable" = mk(i_HW, i_HN))
anova(rm4, L = L)

sample.sizes <- as.data.frame(table(data$Tannin_type))
sample.sizes
sample.sizes <- as.data.frame(table(data$Plant_type, data$Tannin_type))
sample.sizes

dat1 <- data.frame(coef = coef(rm3), se = sqrt(diag(vcov(rm3))), Type = c("Condensed tannins", 
                                                                          "Hydrolysable tannins", 
                                                                          "Total tannins"))
dat1$lower <- dat1$coef - 1.96 * dat1$se
dat1$upper <- dat1$coef + 1.96 * dat1$se
dat2 <- data.frame(coef = coef(rm4), se = sqrt(diag(vcov(rm4))), Type = c("Condensed tannins - Non-woody", 
                                                                          "Condensed tannins - Woody",
                                                                          "Hydrolysable tannins - Woody", 
                                                                          "Total tannins - Non-woody",
                                                                          "Total tannins - Woody"))
dat2$lower <- dat2$coef - 1.96 * dat2$se
dat2$upper <- dat2$coef + 1.96 * dat2$se
dat <- rbind(dat1, dat2)
dat

p1 <- ggplot(dat, aes(x = coef, y = Type)) +
  geom_errorbarh(mapping = aes(xmin = lower, xmax = upper), 
                 width = 0, size=1, colour = "grey10", height = 0.2) +
  geom_point(size = 2, shape = 21, color = "grey10", stroke = 1, 
             mapping = aes(x = coef, y = Type, colour = Type, fill = Type)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "steelblue", size = 1) +
  theme_classic() +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        axis.text.x = element_text(size = 12, hjust = 0.55),
        axis.text.y = element_text(size = 12, hjust = 1),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  scale_fill_manual(values = c( "Condensed tannins - Non-woody" = "#2b6a99",                       
                                "Total tannins - Non-woody" = "#2b6a99",                        
                                "Hydrolysable tannins - Woody" = "#f16c23",                  
                                "Condensed tannins - Woody" = "#f16c23",
                                "Total tannins - Woody"  = "#f16c23",
                                "Hydrolysable tannins" = "#bfbfbf",
                                "Condensed tannins"  = "#bfbfbf",
                                "Total tannins" = "#bfbfbf"),
                    breaks = c("Condensed tannins - Non-woody",                       
                               "Total tannins - Non-woody",                        
                               "Hydrolysable tannins - Woody",                  
                               "Condensed tannins - Woody",
                               "Total tannins - Woody",
                               "Hydrolysable tannins",
                               "Condensed tannins",
                               "Total tannins")) +
  scale_y_discrete(limits = c("Condensed tannins - Non-woody",                       
                              "Total tannins - Non-woody",                        
                              "Hydrolysable tannins - Woody",                  
                              "Condensed tannins - Woody",
                              "Total tannins - Woody",
                              "Hydrolysable tannins",
                              "Condensed tannins",
                              "Total tannins")) +
  xlim(c(-1, 1)) +
  xlab("Effect size (ln"~italic(RR)~")") +
  ylab("Tannin type and plant growth form") +
  coord_fixed(ratio = 0.4)
p1
ggsave('./Figure 2a.pdf', p1, height = 125, width = 200, units = c("mm"))



# Figure 2b
rm5 <- rma.mv(lnRR, rvar,
              mods = ~ paste(Leaf_type) - 1,
              random = list( ~ 1 | Publication/Effect_ID, ~ 1 | Plant_species),
              method = "REML", 
              data = data)
rm5

bn <- names(coef(rm5)); p <- length(bn)
i_M <- which(bn == "paste(Leaf_type)Mixed leaves")
i_T <- which(bn == "paste(Leaf_type)Treated leaves")
i_U <- which(bn == "paste(Leaf_type)Untreated leaves")
mk <- function(plus, minus){ v <- numeric(p); v[plus] <- 1; v[minus] <- -1; v }
L <- rbind("Treated − Untreated" = mk(i_T, i_U),
           "Mixed − Treated"     = mk(i_M, i_T),
           "Mixed − Untreated"   = mk(i_M, i_U))
anova(rm5, L = L)

rm6 <- rma.mv(lnRR, rvar,
              mods = ~ paste(Leaf_type, Plant_type) - 1,
              random = list( ~ 1 | Publication/Effect_ID, ~ 1 | Plant_species),
              method = "REML", 
              data = data)
rm6

bn <- names(coef(rm6))
p  <- length(bn)
ix <- function(pat) which(grepl(pat, bn))
mk <- function(plus, minus){ v <- rep(0, p); v[plus] <- 1; v[minus] <- -1; v }
i_MW <- ix("Mixed.*Woody")
i_MN <- ix("Mixed.*Non-woody")
i_TW <- ix("Treated.*Woody")
i_TN <- ix("Treated.*Non-woody")
i_UW <- ix("Untreated.*Woody")
i_UN <- ix("Untreated.*Non-woody")
L <- rbind(if(length(i_MW)==1 && length(i_MN)==1) "Woody − Non-woody | Mixed"     = mk(i_MW, i_MN),
           if(length(i_TW)==1 && length(i_TN)==1) "Woody − Non-woody | Treated"   = mk(i_TW, i_TN),
           if(length(i_UW)==1 && length(i_UN)==1) "Woody − Non-woody | Untreated" = mk(i_UW, i_UN))
anova(rm6, L = L)

sample.sizes <- as.data.frame(table(data$Leaf_type))
sample.sizes
sample.sizes <- as.data.frame(table(data$Plant_type, data$Leaf_type))
sample.sizes

dat1 <- data.frame(coef = coef(rm5), se = sqrt(diag(vcov(rm5))), Type = c("Mixed leaves",
                                                                          "Treated leaves", 
                                                                          "Untreated leaves"))
dat1$lower <- dat1$coef - 1.96 * dat1$se
dat1$upper <- dat1$coef + 1.96 * dat1$se
dat1
dat2 <- data.frame(coef = coef(rm6), se = sqrt(diag(vcov(rm6))), Type = c("Mixed leaves - Non-woody",
                                                                          "Mixed leaves - Woody",
                                                                          "Treated leaves - Non-woody", 
                                                                          "Treated leaves - Woody", 
                                                                          "Untreated leaves - Non-woody",
                                                                          "Untreated leaves - Woody"))
dat2$lower <- dat2$coef - 1.96 * dat2$se
dat2$upper <- dat2$coef + 1.96 * dat2$se
dat2
dat <- rbind(dat1, dat2)
dat

p2 <- ggplot(dat, aes(x = coef, y = Type)) +
  geom_errorbarh(mapping = aes(xmin = lower, xmax = upper), 
                 width = 0, size=1, colour = "grey10", height = 0.2) +
  geom_point(size = 2, shape = 21, color = "grey10", stroke = 1, 
             mapping = aes(x = coef, y = Type, colour = Type, fill = Type)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "steelblue", size = 1) +
  theme_classic() +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        axis.text.x = element_text(size = 12, hjust = 0.55),
        axis.text.y = element_text(size = 12, hjust = 1),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  scale_fill_manual(values = c( "Untreated leaves - Non-woody" = "#2b6a99",              
                                "Treated leaves - Non-woody" = "#2b6a99",                        
                                "Mixed leaves - Non-woody" = "#2b6a99",                       
                                "Untreated leaves - Woody" = "#f16c23",
                                "Treated leaves - Woody" = "#f16c23", 
                                "Mixed leaves - Woody" = "#f16c23",
                                "Untreated leaves" = "#bfbfbf",
                                "Treated leaves" = "#bfbfbf",
                                "Mixed leaves" = "#bfbfbf"),
                    breaks = c("Untreated leaves - Non-woody",              
                               "Treated leaves - Non-woody",                        
                               "Mixed leaves - Non-woody",                       
                               "Untreated leaves - Woody",
                               "Treated leaves - Woody", 
                               "Mixed leaves - Woody",
                               "Untreated leaves",
                               "Treated leaves",
                               "Mixed leaves")) +
  scale_y_discrete(limits = c("Untreated leaves - Non-woody",              
                              "Treated leaves - Non-woody",                        
                              "Mixed leaves - Non-woody",                       
                              "Untreated leaves - Woody",
                              "Treated leaves - Woody", 
                              "Mixed leaves - Woody",
                              "Untreated leaves",
                              "Treated leaves",
                              "Mixed leaves")) +
  xlim(c(-1, 1)) +
  xlab("Effect size (ln"~italic(RR)~")") +
  ylab("Leaf status and plant growth form") +
  coord_fixed(ratio = 0.3555555555556)
p2
ggsave('./Figure 2b.pdf', p2, height = 125, width = 200, units = c("mm"))

# Figure 2c
rm7 <- rma.mv(lnRR, rvar,
              mods = ~ paste(Leaf_maturity) - 1,
              random = list( ~ 1 | Publication/Effect_ID, ~ 1 | Plant_species),
              method = "REML", 
              data = data)
rm7

bn <- names(coef(rm7)); p <- length(bn)
i_M <- which(bn == "paste(Leaf_maturity)Mixed leaves")
i_P <- which(bn == "paste(Leaf_maturity)Preexisting leaves")
i_R <- which(bn == "paste(Leaf_maturity)Re-growth leaves")
mk <- function(plus, minus){ v <- numeric(p); v[plus] <- 1; v[minus] <- -1; v }
L <- rbind("Re-growth − Preexisting" = mk(i_R, i_P),
           "Preexisting − Mixed"     = mk(i_P, i_M),
           "Re-growth − Mixed"       = mk(i_R, i_M))
anova(rm7, L = L)

rm8 <- rma.mv(lnRR, rvar,
               mods = ~ paste(Leaf_maturity, Plant_type) - 1,
              random = list( ~ 1 | Publication/Effect_ID, ~ 1 | Plant_species),
               method = "REML", 
               data = data)
rm8

bn <- names(coef(rm8)); p <- length(bn)
ix <- function(pat) which(grepl(pat, bn))
mk <- function(plus, minus){ v <- rep(0, p); v[plus] <- 1; v[minus] <- -1; v }

i_MW <- ix("Mixed.*Woody")
i_MN <- ix("Mixed.*Non-woody")
i_PW <- ix("Preexisting.*Woody")
i_PN <- ix("Preexisting.*Non-woody")
i_RW <- ix("Re[- ]?growth.*Woody")
i_RN <- ix("Re[- ]?growth.*Non-woody")
L <- rbind(if(length(i_MW)==1 && length(i_MN)==1) "Woody − Non-woody | Mixed"       = mk(i_MW, i_MN),
           if(length(i_PW)==1 && length(i_PN)==1) "Woody − Non-woody | Preexisting" = mk(i_PW, i_PN))
anova(rm8, L = L)

sample.sizes <- as.data.frame(table(data$Leaf_maturity))
sample.sizes
sample.sizes <- as.data.frame(table(data$Plant_type, data$Leaf_maturity))
sample.sizes

dat1 <- data.frame(coef = coef(rm7), se = sqrt(diag(vcov(rm7))), Type = c("Mixed leaves", 
                                                                          "Preexisting leaves", 
                                                                          "Re-growth leaves"))
dat1$lower <- dat1$coef - 1.96 * dat1$se
dat1$upper <- dat1$coef + 1.96 * dat1$se
dat1
dat2 <- data.frame(coef = coef(rm8), se = sqrt(diag(vcov(rm8))), Type = c("Mixed leaves - Non-woody",
                                                                          "Mixed leaves - Woody", 
                                                                          "Preexisting leaves - Non-woody", 
                                                                          "Preexisting leaves - Woody",
                                                                          "Re-growth leaves - Woody"))
dat2$lower <- dat2$coef - 1.96 * dat2$se
dat2$upper <- dat2$coef + 1.96 * dat2$se
dat2
dat <- rbind(dat1, dat2)
dat

p3 <- ggplot(dat, aes(x = coef, y = Type)) +
  geom_errorbarh(mapping = aes(xmin = lower, xmax = upper), 
                 width = 0, size=1, colour = "grey10", height = 0.2) +
  geom_point(size = 2, shape = 21, color = "grey10", stroke = 1, 
             mapping = aes(x = coef, y = Type, colour = Type, fill = Type)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "steelblue", size = 1) +
  theme_classic() +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        axis.text.x = element_text(size = 12, hjust = 0.55),
        axis.text.y = element_text(size = 12, hjust = 1),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  scale_fill_manual(values = c("Preexisting leaves - Non-woody" = "#2b6a99", 
                               "Mixed leaves - Non-woody" = "#2b6a99",
                               "Re-growth leaves - Woody" = "#f16c23",
                               "Preexisting leaves - Woody" = "#f16c23", 
                               "Mixed leaves - Woody" = "#f16c23",
                               "Re-growth leaves" = "#bfbfbf",
                               "Preexisting leaves" = "#bfbfbf",
                               "Mixed leaves" = "#bfbfbf"),
                    breaks = c("Preexisting leaves - Non-woody", 
                               "Mixed leaves - Non-woody",
                               "Re-growth leaves - Woody",
                               "Preexisting leaves - Woody", 
                               "Mixed leaves - Woody",
                               "Re-growth leaves",
                               "Preexisting leaves",
                               "Mixed leaves")) +
  scale_y_discrete(limits = c("Preexisting leaves - Non-woody", 
                              "Mixed leaves - Non-woody",
                              "Re-growth leaves - Woody",
                              "Preexisting leaves - Woody", 
                              "Mixed leaves - Woody",
                              "Re-growth leaves",
                              "Preexisting leaves",
                              "Mixed leaves")) +
  xlim(c(-1, 1)) +
  xlab("Effect size (ln"~italic(RR)~")") +
  ylab("Leaf developmental stage and plant growth form") +
  coord_fixed(ratio = 0.4)
p3
ggsave('./Figure 2c.pdf', p3, height = 125, width = 200, units = c("mm"))

p4 <- p1 + p2 + p3
p4
ggsave('./Figure 2 - 0829.pdf', p4, height = 125, width = 600, units = c("mm"))












































### Figure 3 ###

# Figure 3a #
rm9 <- rma.mv(lnRR, rvar,
              mods = ~ paste(Treatment_type) - 1,
              random = list( ~ 1 | Publication/Effect_ID, ~ 1 | Plant_species),
              method = "REML", 
              data = data)
rm9

bn <- names(coef(rm9)); p <- length(bn)
i_A <- which(bn == "paste(Treatment_type)Artifical")
i_C <- which(bn == "paste(Treatment_type)Chemical")
i_H <- which(bn == "paste(Treatment_type)Herbivore")
i_J <- which(bn == "paste(Treatment_type)Joint")
mk <- function(plus, minus){ v <- numeric(p); v[plus] <- 1; v[minus] <- -1; v }
L <- rbind("Artifical − Chemical"  = mk(i_A, i_C),
           "Artifical − Herbivore" = mk(i_A, i_H),
           "Artifical − Joint"     = mk(i_A, i_J),
           "Herbivore − Chemical"  = mk(i_H, i_C),
           "Herbivore − Joint"     = mk(i_H, i_J),
           "Joint − Chemical"      = mk(i_J, i_C))
anova(rm9, L = L)

rm10 <- rma.mv(lnRR, rvar,
               mods = ~ paste(Treatment_type, Plant_type) - 1,
               random = list( ~ 1 | Publication/Effect_ID, ~ 1 | Plant_species),
               method = "REML", 
               data = data)
rm10

bn <- names(coef(rm10)); p <- length(bn)
pick <- function(method, form) {
  which(bn == paste0("paste(Treatment_type, Plant_type)", method, " ", form))
}

i_AN <- pick("Artifical", "Non-woody")
i_AW <- pick("Artifical", "Woody")
i_CN <- pick("Chemical",  "Non-woody")
i_CW <- pick("Chemical",  "Woody")
i_HN <- pick("Herbivore", "Non-woody")
i_HW <- pick("Herbivore", "Woody")
i_JN <- pick("Joint",     "Non-woody")
i_JW <- pick("Joint",     "Woody")

mk <- function(plus, minus){ v <- numeric(p); v[plus] <- 1; v[minus] <- -1; v }
L <- rbind("Woody − Non-woody | Artificial" = mk(i_AW, i_AN),
           "Woody − Non-woody | Chemical"   = mk(i_CW, i_CN),
           "Woody − Non-woody | Herbivore"  = mk(i_HW, i_HN),
           "Woody − Non-woody | Joint"      = mk(i_JW, i_JN))

anova(rm10, L = L)

sample.sizes <- as.data.frame(table(data$Treatment_type))
sample.sizes
sample.sizes <- as.data.frame(table(data$Plant_type, data$Treatment_type))
sample.sizes

dat1 <- data.frame(coef = coef(rm9), se = sqrt(diag(vcov(rm9))), Type = c("Artifical treatment", 
                                                                          "Chemical treatment",
                                                                          "Herbivore treatment", 
                                                                          "Joint treatment"))
dat1$lower <- dat1$coef - 1.96 * dat1$se
dat1$upper <- dat1$coef + 1.96 * dat1$se
dat1
dat2 <- data.frame(coef = coef(rm10), se = sqrt(diag(vcov(rm10))), Type = c("Artifical treatment - Non-woody", 
                                                                            "Artifical treatment - Woody",
                                                                            "Chemical treatment - Non-woody",
                                                                            "Chemical treatment - Woody",
                                                                            "Herbivore treatment - Non-woody",
                                                                            "Herbivore treatment - Woody", 
                                                                            "Joint treatment - Non-woody",
                                                                            "Joint treatment - Woody"))
dat2$lower <- dat2$coef - 1.96 * dat2$se
dat2$upper <- dat2$coef + 1.96 * dat2$se
dat2
dat <- rbind(dat1, dat2)
dat

p1 <- ggplot(dat, aes(x = coef, y = Type)) +
  geom_errorbarh(mapping = aes(xmin = lower, xmax = upper), 
                 width = 0, size=1, colour = "grey10", height = 0.2) +
  geom_point(size = 2, shape = 21, color = "grey10", stroke = 1, 
             mapping = aes(x = coef, y = Type, colour = Type, fill = Type)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "steelblue", size = 1) +
  theme_classic() +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        axis.text.x = element_text(size = 12, hjust = 0.55),
        axis.text.y = element_text(size = 12, hjust = 1),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  scale_fill_manual(values = c("Joint treatment - Non-woody" = "#2b6a99", 
                               "Chemical treatment - Non-woody" = "#2b6a99",
                               "Herbivore treatment - Non-woody" = "#2b6a99", 
                               "Artifical treatment - Non-woody" = "#2b6a99",
                               "Joint treatment - Woody" = "#f16c23",
                               "Chemical treatment - Woody" = "#f16c23",
                               "Herbivore treatment - Woody" = "#f16c23",
                               "Artifical treatment - Woody" = "#f16c23",
                               "Joint treatment" = "#bfbfbf",
                               "Chemical treatment" = "#bfbfbf",
                               "Herbivore treatment" = "#bfbfbf",
                               "Artifical treatment" = "#bfbfbf"),
                    breaks = c("Joint treatment - Non-woody",
                               "Chemical treatment - Non-woody",
                               "Herbivore treatment - Non-woody",
                               "Artifical treatment - Non-woody",
                               "Joint treatment - Woody",
                               "Chemical treatment - Woody",
                               "Herbivore treatment - Woody",
                               "Artifical treatment - Woody",
                               "Joint treatment",
                               "Chemical treatment",
                               "Herbivore treatment",
                               "Artifical treatment")) +
  scale_y_discrete(limits = c("Joint treatment - Non-woody",
                              "Chemical treatment - Non-woody",
                              "Herbivore treatment - Non-woody",
                              "Artifical treatment - Non-woody",
                              "Joint treatment - Woody",
                              "Chemical treatment - Woody",
                              "Herbivore treatment - Woody",
                              "Artifical treatment - Woody",
                              "Joint treatment",
                              "Chemical treatment",
                              "Herbivore treatment",
                              "Artifical treatment")) +
  xlim(c(-1, 1)) +
  xlab("Effect size (ln"~italic(RR)~")") +
  ylab("Induction treatment and plant growth form") +
  coord_fixed(ratio = 0.5)
p1
ggsave('./Figure 3a.pdf', p1, height = 260, width = 175, units = c("mm"))



# Figure 3b #
data1 <- data %>%
  filter(Treatment_type == "Artifical") %>%
  droplevels()                              

rm11 <- rma.mv(lnRR, rvar,
               mods = ~ paste(Treatment_typeII) - 1,
               random = list( ~ 1 | Publication/Effect_ID, ~ 1 | Plant_species),
               method = "REML", 
               data = data1)
rm11

bn <- names(coef(rm11)); p <- length(bn)
i_D <- which(bn == "paste(Treatment_typeII)Defoliation")
i_Pr <- which(bn == "paste(Treatment_typeII)Pruning")
i_Pu <- which(bn == "paste(Treatment_typeII)Punching")
mk <- function(plus, minus){ v <- numeric(p); v[plus] <- 1; v[minus] <- -1; v }
L <- rbind("Defoliation − Pruning"   = mk(i_D,  i_Pr),
           "Defoliation − Punching"  = mk(i_D,  i_Pu),
           "Pruning − Punching"      = mk(i_Pr, i_Pu))
anova(rm11, L = L)

rm12 <- rma.mv(lnRR, rvar,
               mods = ~ paste(Treatment_typeII, Plant_type) - 1,
               random = list( ~ 1 | Publication/Effect_ID, ~ 1 | Plant_species),
               method = "REML", 
               data = data1)
rm12

bn <- names(coef(rm12)); p <- length(bn)
pick <- function(method, form) {
  which(bn == paste0("paste(Treatment_typeII, Plant_type)", method, " ", form))
}

i_DN <- pick("Defoliation", "Non-woody")
i_DW <- pick("Defoliation", "Woody")
i_PuN <- pick("Punching", "Non-woody")
i_PuW <- pick("Punching", "Woody")
mk <- function(plus, minus){ v <- numeric(p); v[plus] <- 1; v[minus] <- -1; v }
L <- rbind("Woody − Non-woody | Defoliation" = mk(i_DW,  i_DN),
           "Woody − Non-woody | Punching"    = mk(i_PuW, i_PuN))
anova(rm12, L = L)

sample.sizes <- as.data.frame(table(data1$Treatment_typeII))
sample.sizes
sample.sizes <- as.data.frame(table(data1$Plant_type, data1$Treatment_typeII))
sample.sizes

dat1 <- data.frame(coef = coef(rm11), se = sqrt(diag(vcov(rm11))), Type = c("Defoliation", 
                                                                            "Pruning",
                                                                            "Punching"))
dat1$lower <- dat1$coef - 1.96 * dat1$se
dat1$upper <- dat1$coef + 1.96 * dat1$se
dat1

dat2 <- data.frame(coef = coef(rm12), se = sqrt(diag(vcov(rm12))), Type = c("Defoliation - Non-woody", 
                                                                            "Defoliation - Woody",
                                                                            "Pruning - Woody",
                                                                            "Punching - Non-woody",
                                                                            "Punching - Woody"))
dat2$lower <- dat2$coef - 1.96 * dat2$se
dat2$upper <- dat2$coef + 1.96 * dat2$se
dat2
dat <- rbind(dat1, dat2)
dat

p2 <- ggplot(dat, aes(x = coef, y = Type)) +
  geom_errorbarh(mapping = aes(xmin = lower, xmax = upper), 
                 width = 0, size=1, colour = "grey10", height = 0.2) +
  geom_point(size = 2, shape = 21, color = "grey10", stroke = 1, 
             mapping = aes(x = coef, y = Type, colour = Type, fill = Type)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "steelblue", size = 1) +
  theme_classic() +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        axis.text.x = element_text(size = 12, hjust = 0.55),
        axis.text.y = element_text(size = 12, hjust = 1),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  scale_fill_manual(values = c("Punching - Non-woody" = "#2b6a99",
                               "Defoliation - Non-woody" = "#2b6a99",
                               "Punching - Woody" = "#f16c23",
                               "Pruning - Woody" = "#f16c23",
                               "Defoliation - Woody" = "#f16c23",
                               "Punching" = "#bfbfbf",
                               "Pruning" = "#bfbfbf",
                               "Defoliation"  = "#bfbfbf"),
                    breaks = c("Punching - Non-woody",
                               "Defoliation - Non-woody",
                               "Punching - Woody",
                               "Pruning - Woody",
                               "Defoliation - Woody",
                               "Punching",
                               "Pruning",
                               "Defoliation")) +
  scale_y_discrete(limits = c("Punching - Non-woody",
                              "Defoliation - Non-woody",
                              "Punching - Woody",
                              "Pruning - Woody",
                              "Defoliation - Woody",
                              "Punching",
                              "Pruning",
                              "Defoliation")) +
  xlim(c(-1.5, 1.5)) +
  xlab("Effect size (ln"~italic(RR)~")") +
  ylab("Specific induction triggers and plant growth form") +
  coord_fixed(ratio = 0.6)
p2

ggsave('./Figure 3b.pdf', p2, height = 125, width = 200, units = c("mm"))



# Figure 3c
data1 <- data %>%
  filter(Treatment_type == "Herbivore") %>%
  droplevels()                              

rm13 <- rma.mv(lnRR, rvar,
               mods = ~ paste(Treatment_typeII) - 1,
               random = list( ~ 1 | Publication/Effect_ID, ~ 1 | Plant_species),
               method = "REML", 
               data = data1)
rm13

bn <- names(coef(rm13))
i_G <- which(bn == "paste(Treatment_typeII)Generalist")
i_S <- which(bn == "paste(Treatment_typeII)Specialist")
L <- rbind("Specialist − Generalist" = {v <- numeric(length(bn)); v[i_S] <- 1; v[i_G] <- -1; v})
anova(rm13, L = L)

rm14 <- rma.mv(lnRR, rvar,
               mods = ~ paste(Treatment_typeII, Plant_type) - 1,
               random = list( ~ 1 | Publication/Effect_ID, ~ 1 | Plant_species),
               method = "REML", 
               data = data1)
rm14

bn <- names(coef(rm14)); p <- length(bn)
pick <- function(level, form) {
  which(bn == paste0("paste(Treatment_typeII, Plant_type)", level, " ", form))
}
i_GN <- pick("Generalist", "Non-woody")
i_GW <- pick("Generalist", "Woody")
i_SN <- pick("Specialist", "Non-woody")
i_SW <- pick("Specialist", "Woody")
mk <- function(plus, minus){ v <- numeric(p); v[plus] <- 1; v[minus] <- -1; v }
L <- rbind("Woody − Non-woody | Generalist" = mk(i_GW, i_GN),
           "Woody − Non-woody | Specialist" = mk(i_SW, i_SN))
anova(rm14, L = L)

sample.sizes <- as.data.frame(table(data1$Treatment_typeII))
sample.sizes
sample.sizes <- as.data.frame(table(data1$Plant_type, data1$Treatment_typeII))
sample.sizes

dat1 <- data.frame(coef = coef(rm13), se = sqrt(diag(vcov(rm13))), Type = c("Generalist", "Specialist"))
dat1$lower <- dat1$coef - 1.96 * dat1$se
dat1$upper <- dat1$coef + 1.96 * dat1$se
dat1
dat2 <- data.frame(coef = coef(rm14), se = sqrt(diag(vcov(rm14))), Type = c("Generalist - Non-woody", 
                                                                            "Generalist - Woody",
                                                                            "Specialist - Non-woody",
                                                                            "Specialist - Woody"))
dat2$lower <- dat2$coef - 1.96 * dat2$se
dat2$upper <- dat2$coef + 1.96 * dat2$se
dat2
dat <- rbind(dat1, dat2)
dat

p3 <- ggplot(dat, aes(x = coef, y = Type)) +
  geom_errorbarh(mapping = aes(xmin = lower, xmax = upper), 
                 width = 0, size=1, colour = "grey10", height = 0.2) +
  geom_point(size = 2, shape = 21, color = "grey10", stroke = 1, 
             mapping = aes(x = coef, y = Type, colour = Type, fill = Type)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "steelblue", size = 1) +
  theme_classic() +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        axis.text.x = element_text(size = 12, hjust = 0.55),
        axis.text.y = element_text(size = 12, hjust = 1),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  scale_fill_manual(values = c("Specialist - Non-woody" = "#2b6a99",
                               "Generalist - Non-woody" = "#2b6a99",
                               "Specialist - Woody" = "#f16c23", 
                               "Generalist - Woody" = "#f16c23",
                               "Specialist" = "#bfbfbf",
                               "Generalist" = "#bfbfbf"),
                    breaks = c("Specialist - Non-woody",
                               "Generalist - Non-woody",
                               "Specialist - Woody", 
                               "Generalist - Woody",
                               "Specialist",
                               "Generalist")) +
  scale_y_discrete(limits = c("Specialist - Non-woody",
                              "Generalist - Non-woody",
                              "Specialist - Woody", 
                              "Generalist - Woody",
                              "Specialist",
                              "Generalist")) +
  xlim(c(-1, 1)) +
  xlab("Effect size (ln"~italic(RR)~")") +
  ylab("Herbivore feeding guild and plant growth form") +
  coord_fixed(ratio = 0.53333333333)
p3
ggsave('./Figure 3c.pdf', p3, height = 125, width = 200, units = c("mm"))



# Figure 3d
data1 <- data %>%
  filter(Treatment_type == "Chemical") %>%
  droplevels()                              

rm15 <- rma.mv(lnRR, rvar,
               mods = ~ paste(Treatment_typeII) - 1,
               random = list( ~ 1 | Publication/Effect_ID, ~ 1 | Plant_species),
               method = "REML", 
               data = data1)
rm15

bn <- names(coef(rm15))
i_JA <- which(bn == "paste(Treatment_typeII)JA")
i_SA <- which(bn == "paste(Treatment_typeII)SA")
L <- rbind("JA − SA" = {v <- numeric(length(bn)); v[i_JA] <- 1; v[i_SA] <- -1; v})
anova(rm15, L = L)

rm16 <- rma.mv(lnRR, rvar,
               mods = ~ paste(Treatment_typeII, Plant_type) - 1,
               random = list( ~ 1 | Publication/Effect_ID, ~ 1 | Plant_species),
               method = "REML", 
               data = data1)
rm16

bn <- names(coef(rm16)); p <- length(bn)
pick <- function(level, form) {
  which(bn == paste0("paste(Treatment_typeII, Plant_type)", level, " ", form))
}
i_JN <- pick("JA", "Non-woody")
i_JW <- pick("JA", "Woody")
i_SN <- pick("SA", "Non-woody")
i_SW <- pick("SA", "Woody")
mk <- function(plus, minus){ v <- numeric(p); v[plus] <- 1; v[minus] <- -1; v }
L <- rbind("Woody − Non-woody | JA" = mk(i_JW, i_JN),
           "Woody − Non-woody | SA" = mk(i_SW, i_SN))
anova(rm16, L = L)

sample.sizes <- as.data.frame(table(data1$Treatment_typeII))
sample.sizes
sample.sizes <- as.data.frame(table(data1$Plant_type, data1$Treatment_typeII))
sample.sizes

dat1 <- data.frame(coef = coef(rm15), se = sqrt(diag(vcov(rm15))), Type = c("JA", "SA"))
dat1$lower <- dat1$coef - 1.96 * dat1$se
dat1$upper <- dat1$coef + 1.96 * dat1$se
dat1
dat2 <- data.frame(coef = coef(rm16), se = sqrt(diag(vcov(rm16))), Type = c("JA - Non-woody", 
                                                                            "JA - Woody",
                                                                            "SA - Non-woody",
                                                                            "SA - Woody"))
dat2$lower <- dat2$coef - 1.96 * dat2$se
dat2$upper <- dat2$coef + 1.96 * dat2$se
dat2
dat <- rbind(dat1, dat2)
dat

p4 <- ggplot(dat, aes(x = coef, y = Type)) +
  geom_errorbarh(mapping = aes(xmin = lower, xmax = upper), 
                 width = 0, size=1, colour = "grey10", height = 0.2) +
  geom_point(size = 2, shape = 21, color = "grey10", stroke = 1, 
             mapping = aes(x = coef, y = Type, colour = Type, fill = Type)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "steelblue", size = 1) +
  theme_classic() +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        axis.text.x = element_text(size = 12, hjust = 0.55),
        axis.text.y = element_text(size = 12, hjust = 1),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  scale_fill_manual(values = c("SA - Non-woody" = "#2b6a99",
                               "JA - Non-woody" = "#2b6a99",
                               "SA - Woody" = "#f16c23", 
                               "JA - Woody" = "#f16c23",
                               "SA" = "#bfbfbf",
                               "JA" = "#bfbfbf"),
                    breaks = c("SA - Non-woody",
                               "JA - Non-woody",
                               "SA - Woody", 
                               "JA - Woody",
                               "SA",
                               "JA")) +
  scale_y_discrete(limits = c("SA - Non-woody",
                              "JA - Non-woody",
                              "SA - Woody", 
                              "JA - Woody",
                              "SA",
                              "JA")) +
  xlim(c(-1, 1)) +
  xlab("Effect size (ln"~italic(RR)~")") +
  ylab("Herbivore feeding guild and plant growth form") +
  coord_fixed(ratio = 0.53333333333)
p4
ggsave('./Figure 3d.pdf', p4, height = 125, width = 200, units = c("mm"))



# Figure 3e
data1 <- data %>%
  filter(Treatment_type == "Joint") %>%
  droplevels()                              

rm17 <- rma.mv(lnRR, rvar,
               mods = ~ paste(Treatment_typeII) - 1,
               random = list( ~ 1 | Publication/Effect_ID, ~ 1 | Plant_species),
               method = "REML", 
               data = data1)
rm17

bn <- names(coef(rm17)); p <- length(bn)
w  <- function(label) which(bn == paste0("paste(Treatment_typeII)", label))
i_DH <- w("Defoliation + Herbivore")
i_DJ <- w("Defoliation + JA")
i_DS <- w("Defoliation + Salivation")
mk <- function(plus, minus){ v <- numeric(p); v[plus] <- 1; v[minus] <- -1; v }
L <- rbind("Defol+Herb − Defol+JA"        = mk(i_DH, i_DJ),
           "Defol+Herb − Defol+Salivation"= mk(i_DH, i_DS),
           "Defol+JA   − Defol+Salivation"= mk(i_DJ, i_DS))
anova(rm17, L = L)

rm18 <- rma.mv(lnRR, rvar,
               mods = ~ paste(Treatment_typeII, Plant_type) - 1,
               random = list( ~ 1 | Publication/Effect_ID, ~ 1 | Plant_species),
               method = "REML", 
               data = data1)
rm18

bn <- names(coef(rm18)); p <- length(bn)
pick <- function(level, form) {
  which(bn == paste0("paste(Treatment_typeII, Plant_type)", level, " ", form))
}
i_JN <- pick("Defoliation + JA",         "Non-woody")
i_JW <- pick("Defoliation + JA",         "Woody")
i_SN <- pick("Defoliation + Salivation", "Non-woody")
i_SW <- pick("Defoliation + Salivation", "Woody")
i_HN <- pick("Defoliation + Herbivore",  "Non-woody")
i_HW <- pick("Defoliation + Herbivore",  "Woody")
mk <- function(plus, minus){ v <- numeric(p); v[plus] <- 1; v[minus] <- -1; v }
L_list <- list()
if (length(i_JW)==1 && length(i_JN)==1)
  L_list[["Woody − Non-woody | Defol + JA"]] <- mk(i_JW, i_JN)
if (length(i_SW)==1 && length(i_SN)==1)
  L_list[["Woody − Non-woody | Defol + Salivation"]] <- mk(i_SW, i_SN)
if (length(i_HW)==1 && length(i_HN)==1)
  L_list[["Woody − Non-woody | Defol + Herbivore"]] <- mk(i_HW, i_HN)  # 只在两格都在时才比较
stopifnot(length(L_list) > 0)
L <- do.call(rbind, L_list)
anova(rm18, L = L)

dsample.sizes <- as.data.frame(table(data1$Treatment_typeII))
sample.sizes
sample.sizes <- as.data.frame(table(data1$Plant_type, data1$Treatment_typeII))
sample.sizes

dat1 <- data.frame(coef = coef(rm17), se = sqrt(diag(vcov(rm17))), Type = c("Defoliation + Herbivore", 
                                                                            "Defoliation + JA",
                                                                            "Defoliation + Salivation"))
dat1$lower <- dat1$coef - 1.96 * dat1$se
dat1$upper <- dat1$coef + 1.96 * dat1$se
dat1
dat2 <- data.frame(coef = coef(rm18), se = sqrt(diag(vcov(rm18))), Type = c("Defoliation + Herbivore Woody", 
                                                                            "Defoliation + JA Non-woody",
                                                                            "Defoliation + JA Woody",
                                                                            "Defoliation + Salivation Non-woody",
                                                                            "Defoliation + Salivation Woody"))
dat2$lower <- dat2$coef - 1.96 * dat2$se
dat2$upper <- dat2$coef + 1.96 * dat2$se
dat2
dat <- rbind(dat1, dat2)
dat

p5 <- ggplot(dat, aes(x = coef, y = Type)) +
  geom_errorbarh(mapping = aes(xmin = lower, xmax = upper), 
                 width = 0, size=1, colour = "grey10", height = 0.2) +
  geom_point(size = 2, shape = 21, color = "grey10", stroke = 1, 
             mapping = aes(x = coef, y = Type, colour = Type, fill = Type)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "steelblue", size = 1) +
  theme_classic() +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        axis.text.x = element_text(size = 12, hjust = 0.55),
        axis.text.y = element_text(size = 12, hjust = 1),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  scale_fill_manual(values = c("Defoliation + JA Non-woody" = "#2b6a99",
                               "Defoliation + Salivation Non-woody" = "#2b6a99",
                               "Defoliation + JA Woody" = "#f16c23",
                               "Defoliation + Salivation Woody" = "#f16c23",
                               "Defoliation + Herbivore Woody" = "#f16c23", 
                               "Defoliation + JA" = "#bfbfbf",
                               "Defoliation + Salivation" = "#bfbfbf",
                               "Defoliation + Herbivore" = "#bfbfbf"),
                    breaks = c("Defoliation + JA Non-woody",
                               "Defoliation + Salivation Non-woody",
                               "Defoliation + JA Woody",
                               "Defoliation + Salivation Woody",
                               "Defoliation + Herbivore Woody", 
                               "Defoliation + JA",
                               "Defoliation + Salivation",
                               "Defoliation + Herbivore")) +
  scale_y_discrete(limits = c("Defoliation + JA Non-woody",
                              "Defoliation + Salivation Non-woody",
                              "Defoliation + JA Woody",
                              "Defoliation + Salivation Woody",
                              "Defoliation + Herbivore Woody", 
                              "Defoliation + JA",
                              "Defoliation + Salivation",
                              "Defoliation + Herbivore")) +
  xlim(c(-1, 1)) +
  xlab("Effect size (ln"~italic(RR)~")") +
  ylab("Herbivore feeding guild and plant growth form") +
  coord_fixed(ratio = 0.4)
p5
ggsave('./Figure 3e.pdf', p5, height = 125, width = 200, units = c("mm"))

p6 <- p2 + p3 + p4 + p5
p6
ggsave('./Figure 3.pdf', p6, height = 250, width = 400, units = c("mm"))

rm9
rm10
rm11
rm12
rm13
rm14
rm15
rm16
rm17
rm18






### Figure 4 ###
# Figure 4a #
setwd("/Users/yaolin/Desktop/My papers/Manuscripts/2024 - Sun - Tannin induction")
data <- read_excel("RawData - Tannin.xlsx")
write.csv(data, file = "data.csv", row.names = FALSE)
data <- read.csv("data.csv")

fit <- lm(formula = lnRR ~ Herbivory_intensity - 1, data = data)
summary(fit)
fit <- lm(formula = lnRR ~ Herbivory_intensity - 1, data = subset(data, Plant_type == "Woody"))
summary(fit)
fit <- lm(formula = lnRR ~ Herbivory_intensity - 1, data = subset(data, Plant_type == "Non-woody"))
summary(fit)

p1 <- ggplot(data, aes(x = Herbivory_intensity, y = lnRR)) +
  geom_point(aes(color = Plant_type), 
             size = 4, 
             shape = 16, 
             stroke = 0.25, 
             alpha = 0.6) +
  geom_smooth(aes(color = "Total", fill = "Total", linetype = "Total"),
              method = "lm", se = TRUE, alpha = 0.15,
              formula = y ~ x - 1) +
  geom_smooth(aes(color = Plant_type, fill = Plant_type, linetype = Plant_type),
              method = "lm", se = TRUE, alpha = 0.15,
              formula = y ~ x - 1) +
  scale_color_manual(values = c("Woody"     = "#f16c23",
                                "Non-woody" = "#2b6a99",
                                "Total" = "#696969")) +
  scale_fill_manual(values = c("Woody"     = "#f16c23",
                               "Non-woody" = "#2b6a99",
                               "Total" = "#696969"), 
                    guide = "none") +
  scale_linetype_manual(values = c("Total" = "solid",
                                   "Woody" = "solid",
                                   "Non-woody" = "dashed"), 
                        guide = "none") +
  guides(color = guide_legend(title = "Plant_type")) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x        = element_text(size = 12, hjust = 0.5),
    axis.title.x       = element_text(size = 12),
    legend.position    = "none",       
    legend.justification = c(0, 1),
    axis.line.x        = element_line(color = "black", linewidth = 0.5),
    axis.line.y        = element_line(color = "black", linewidth = 0.5),
    axis.ticks         = element_line(color = "black", linewidth = 0.5),
    axis.title         = element_text(size = 12),
    axis.text          = element_text(size = 12)
  ) +
  xlab("Herbivory intensity (%)") +
  ylab("Effect size (ln"~italic(RR)~")") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  geom_text(
    data = data.frame(x = 4, y = -1.25),
    aes(x = x, y = y),
    label = "italic(R)[Total]^2 == '0.137' * ',' ~ italic(p) < 0.001",
    parse = TRUE,
    color = "#696969",
    size  = 4,
    inherit.aes = FALSE,
    hjust = 0) +
  geom_text(
    data = data.frame(x = 4, y = -1.5),
    aes(x = x, y = y),
    label = "italic(R)[Woody]^2 == '0.274' * ',' ~ italic(p) < 0.001",
    parse = TRUE,
    color = "#f16c23",
    size  = 4,
    inherit.aes = FALSE,
    hjust = 0) +
  geom_text(
    data = data.frame(x = 4, y = -1.75),
    aes(x = x, y = y),
    label = "italic(R)[Non-woody]^2 == '0.001' * ',' ~ italic(p) == 0.778",
    parse = TRUE,
    color = "#2b6a99",
    size  = 4,
    inherit.aes = FALSE,
    hjust = 0)
p1
p1 <- ggMarginal(
  p1,
  type        = "density",
  margins     = "both",
  groupColour = TRUE,
  groupFill   = TRUE,
  alpha       = 0.4
)
print(p1)
ggsave("Figure 4a.pdf", p1, height = 150, width = 400, units = "mm")



#Figure 4b
setwd("/Users/yaolin/Desktop/My papers/Manuscripts/2024 - Sun - Tannin induction")
data <- read_excel("RawData - Tannin.xlsx")
write.csv(data, file = "data.csv", row.names = FALSE)
data <- read.csv("data.csv")

fit <- lm(formula = lnRR ~ Treatment_interval, data = data)
summary(fit)
fit <- lm(formula = lnRR ~ Treatment_interval, data = subset(data, Plant_type == "Woody"))
summary(fit)
fit <- lm(formula = lnRR ~ Treatment_interval, data = subset(data, Plant_type == "Non-woody"))
summary(fit)

p2 <- ggplot(data, aes(x = Treatment_interval, y = lnRR)) +
  geom_point(aes(color = Plant_type), 
             size = 4, 
             shape = 16, 
             stroke = 0.25, 
             alpha = 0.6) +
  geom_smooth(aes(color = "Total", fill = "Total", linetype = "Total"),
              method = "lm", se = TRUE, alpha = 0.15,
              formula = y ~ x) +
  geom_smooth(aes(color = Plant_type, fill = Plant_type, linetype = Plant_type),
              method = "lm", se = TRUE, alpha = 0.15,
              formula = y ~ x) +
  scale_color_manual(values = c("Woody"     = "#f16c23",
                                "Non-woody" = "#2b6a99",
                                "Total" = "#696969")) +
  scale_fill_manual(values = c("Woody"     = "#f16c23",
                               "Non-woody" = "#2b6a99",
                               "Total" = "#696969"), 
                    guide = "none") +
  scale_linetype_manual(values = c("Total" = "dashed",
                                   "Woody" = "dashed",
                                   "Non-woody" = "solid"), 
                        guide = "none") +
  guides(color = guide_legend(title = "Plant_type")) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x        = element_text(size = 12, hjust = 0.5),
    axis.title.x       = element_text(size = 12),
    legend.position    = "none",       
    legend.justification = c(0, 1),
    axis.line.x        = element_line(color = "black", linewidth = 0.5),
    axis.line.y        = element_line(color = "black", linewidth = 0.5),
    axis.ticks         = element_line(color = "black", linewidth = 0.5),
    axis.title         = element_text(size = 12),
    axis.text          = element_text(size = 12)
  ) +
  xlab("Treatment interval (days)") +
  ylab("Effect size (ln"~italic(RR)~")") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  geom_text(
    data = data.frame(x = 4, y = -1.25),
    aes(x = x, y = y),
    label = "italic(R)[Total]^2 == '0.002' * ',' ~ italic(p) == 0.382",
    parse = TRUE,
    color = "#696969",
    size  = 4,
    inherit.aes = FALSE,
    hjust = 0) +
  geom_text(
    data = data.frame(x = 4, y = -1.5),
    aes(x = x, y = y),
    label = "italic(R)[Non-woody]^2 == '0.000' * ',' ~ italic(p) == 0.834",
    parse = TRUE,
    color = "#f16c23",
    size  = 4,
    inherit.aes = FALSE,
    hjust = 0) +
  geom_text(
    data = data.frame(x = 4, y = -1.75),
    aes(x = x, y = y),
    label = "italic(R)[Non-woody]^2 == '0.078' * ',' ~ italic(p) == 0.004",
    parse = TRUE,
    color = "#2b6a99",
    size  = 4,
    inherit.aes = FALSE,
    hjust = 0)
p2
p2 <- ggMarginal(
  p2,
  type        = "density",
  margins     = "both",
  groupColour = TRUE,
  groupFill   = TRUE,
  alpha       = 0.4
)

print(p2)
ggsave("Figure 4b.pdf", p2, height = 150, width = 400, units = "mm")



#Figure 4c
setwd("/Users/yaolin/Desktop/My papers/Manuscripts/2024 - Sun - Tannin induction")
data <- read_excel("RawData - Tannin.xlsx")
write.csv(data, file = "data.csv", row.names = FALSE)
data <- read.csv("data.csv")

fit <- lm(formula = lnRR ~ Measurement_interval, data = data)
summary(fit)
fit <- lm(formula = lnRR ~ Measurement_interval, data = subset(data, Plant_type == "Woody"))
summary(fit)
fit <- lm(formula = lnRR ~ Measurement_interval, data = subset(data, Plant_type == "Non-woody"))
summary(fit)

p3 <- ggplot(data, aes(x = Measurement_interval, y = lnRR)) +
  geom_point(aes(color = Plant_type), 
             size = 4, 
             shape = 16, 
             stroke = 0.25, 
             alpha = 0.6) +
  geom_smooth(aes(color = "Total", fill = "Total", linetype = "Total"),
              method = "lm", se = TRUE, alpha = 0.15,
              formula = y ~ x) +
  geom_smooth(aes(color = Plant_type, fill = Plant_type, linetype = Plant_type),
              method = "lm", se = TRUE, alpha = 0.15,
              formula = y ~ x) +
  scale_color_manual(values = c("Woody"     = "#f16c23",
                                "Non-woody" = "#2b6a99",
                                "Total" = "#696969")) +
  scale_fill_manual(values = c("Woody"     = "#f16c23",
                               "Non-woody" = "#2b6a99",
                               "Total" = "#696969"), 
                    guide = "none") +
  scale_linetype_manual(values = c("Total" = "dashed",
                                   "Woody" = "dashed",
                                   "Non-woody" = "dashed"), 
                        guide = "none") +
  guides(color = guide_legend(title = "Plant_type")) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x        = element_text(size = 12, hjust = 0.5),
    axis.title.x       = element_text(size = 12),
    legend.position    = "none",       
    legend.justification = c(0, 1),
    axis.line.x        = element_line(color = "black", linewidth = 0.5),
    axis.line.y        = element_line(color = "black", linewidth = 0.5),
    axis.ticks         = element_line(color = "black", linewidth = 0.5),
    axis.title         = element_text(size = 12),
    axis.text          = element_text(size = 12)
  ) +
  xlab("Measurement interval (days)") +
  ylab("Effect size (ln"~italic(RR)~")") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  geom_text(
    data = data.frame(x = 4, y = -1.25),
    aes(x = x, y = y),
    label = "italic(R)[Total]^2 == '0.002' * ',' ~ italic(p) == 0.384",
    parse = TRUE,
    color = "#696969",
    size  = 4,
    inherit.aes = FALSE,
    hjust = 0) +
  geom_text(
    data = data.frame(x = 4, y = -1.5),
    aes(x = x, y = y),
    label = "italic(R)[Non-woody]^2 == '0.001' * ',' ~ italic(p) == 0.700",
    parse = TRUE,
    color = "#f16c23",
    size  = 4,
    inherit.aes = FALSE,
    hjust = 0) +
  geom_text(
    data = data.frame(x = 4, y = -1.75),
    aes(x = x, y = y),
    label = "italic(R)[Non-woody]^2 == '0.009' * ',' ~ italic(p) == 0.315",
    parse = TRUE,
    color = "#2b6a99",
    size  = 4,
    inherit.aes = FALSE,
    hjust = 0)
p3
p3 <- ggMarginal(
  p3,
  type        = "density",
  margins     = "both",
  groupColour = TRUE,
  groupFill   = TRUE,
  alpha       = 0.4
)
print(p3)
ggsave("Figure 4c.pdf", p3, height = 150, width = 400, units = "mm")


p4 <- grid.arrange(p1, p2, p3, ncol = 1)
p4
ggsave('./Figure 4.pdf', p4, height = 400, width = 150, units = c("mm"))










### Figure 5 Test the bias ###
# Figure 5a
setwd("/Users/yaolin/Desktop/My papers/Manuscripts/2024 - Sun - Tannin induction")
data <- read_excel("RawData - Tannin.xlsx")
write.csv(data, file = "data.csv", row.names = FALSE)
data <- read.csv("data.csv")

standard.model <- rma(lnRR, rvar, 
                      mods = ~ 1, 
                      data = data)
regtest(standard.model)

forest.model <- rma.mv(lnRR, rvar, 
                       mods = ~ 1,
                       random = list(~ 1 | Plant_species),
                       method = "REML",
                       data = data)
forest.model 

# Obtain residuals
resstandards <- rstandard.rma.mv(forest.model, type = "response")
# Obtain grand mean effect size 
grand.mean <- as.numeric(forest.model$b) 
# Create new df with residuals replacing raw
df.forest.model <- data
df.forest.model$lnRR <- resstandards$resid + grand.mean 
df.forest.model$sei <- resstandards$se
# Funnel plot for all outcome classes
make.funnel <- function(dataset, model) {
  apatheme <- theme_bw() +  
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border     = element_blank(),
          axis.line        = element_line(),
          text             = element_text(family = 'Times'),
          legend.position  = 'none')
  estimate <- as.numeric(model$b)
  SE       <- model$se
  se.seq <- seq(
    from = 0, 
    to   = max(sqrt(dataset$rvar), na.rm = TRUE),
    by   = 0.001
  )
  
  dfCI <- data.frame(
    ll95     = estimate - (1.96 * se.seq),
    ul95     = estimate + (1.96 * se.seq),
    ll99     = estimate - (3.29 * se.seq),
    ul99     = estimate + (3.29 * se.seq),
    se.seq   = se.seq,
    meanll95 = estimate - (1.96 * SE),
    meanul95 = estimate + (1.96 * SE)
  )
  
  ggplot(dataset, aes(x = sqrt(rvar), y = lnRR)) +
    geom_point(color = '#696969', shape = 16, size = 4, alpha = 0.6) +
    xlab("Standard error") +
    ylab("Effect size (ln"~italic(RR)~")") +
    geom_line(aes(x = se.seq, y = ll95), data = dfCI, linetype = 'dotted') +
    geom_line(aes(x = se.seq, y = ul95), data = dfCI, linetype = 'dotted') +
    geom_line(aes(x = se.seq, y = ll99), data = dfCI, linetype = 'dashed') +
    geom_line(aes(x = se.seq, y = ul99), data = dfCI, linetype = 'dashed') +
    geom_segment(
      aes(x = min(se.seq), y = meanll95, xend = max(se.seq), yend = meanll95),
      data    = dfCI,
      colour  = "steelblue",
      linetype= 'dotdash',
      size    = 0.75
    ) +
    geom_segment(
      aes(x = min(se.seq), y = meanul95, xend = max(se.seq), yend = meanul95),
      data    = dfCI,
      colour  = "steelblue",
      linetype= 'dotdash',
      size    = 0.75
    ) +
    scale_x_reverse() +
    coord_flip() +
    theme_bw() +
    theme(
      panel.spacing      = unit(0.5, "lines"),
      panel.border       = element_blank(),
      text               = element_text(size = 12),
      axis.line          = element_line(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      legend.text        = element_text(size = 12),
      legend.title       = element_text(size = 12, face = "bold"),
      axis.title.x       = element_text(hjust = 0.5, size = 12),
      axis.title.y       = element_text(size = 12),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12)
    ) +
    geom_text(
      data = data.frame(x = 0 , y = -4.5),
      aes(x = x, y = y),
      label = "italic(z) == '0.124' * ',' ~ italic(p) == 0.902",
      parse = TRUE,
      color = "#696969",
      size  = 4,
      inherit.aes = FALSE,
      hjust = 0) +
    geom_text(data = data.frame(x = 0 , y = 4.5, label = "A"),
              aes(x = x, y = y, label = label),
              color = "black", fontface = "bold", size = 6,
              inherit.aes = FALSE)
    
}

funnel.plot <- make.funnel(data, forest.model)
print(funnel.plot)
ggsave("Figure 5a.pdf", funnel.plot, height = 125, width = 100, units = "mm")



# Figure 5b # 
# the impact by IF
fitIF <- lm(lnRR ~ Impact_factor, data = data)
summary(fitIF)

IF.plot <- ggplot(data = data, aes(x = Impact_factor, y = lnRR)) +
  geom_point(color = '#696969', shape = 16, size = 4, alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = 2, colour = "steelblue", size = 0.75) +
  scale_x_log10(limits = c(-5, 15), breaks = c(0, 1, 2, 5, 10, 15)) +
  geom_smooth(method  = "lm", se = TRUE, color = "#696969", fill = "#696969", alpha   = 0.2, linetype = "dashed") +
  labs(size = 'Weight (%)', y = "Effect size (ln"~italic(RR)~")", x = 'Journal impact factor') + 
  guides(size = F) +
  scale_x_continuous(breaks = seq(0, 15, by = 3)) +
  theme_bw() +
  theme(panel.spacing = unit(0.5, "lines"),
        panel.border= element_blank(),
        text = element_text(size = 12),
        axis.line=element_line(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.text = element_text(size = 12),
        legend.title=element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)
        ) +
  geom_text(
    data = data.frame(x = 0 , y = 2),
    aes(x = x, y = y),
    label = "italic(R)^2 == '0.004' * ',' ~ italic(p) == 0.151",
    parse = TRUE,
    color = "#696969",
    size  = 4,
    inherit.aes = FALSE,
    hjust = 0) +
  geom_text(data = data.frame(x = 11 , y = 2, label = "B"),
            aes(x = x, y = y, label = label),
            color = "black", fontface = "bold", size = 6,
            inherit.aes = FALSE) +
  scale_y_continuous(labels = number_format(accuracy = 0.1))


IF.plot
ggsave("Figure 5b.pdf", IF.plot, height = 125, width = 100, units = "mm")


# Figure 5c # the impact by published year
fityear <- lm(lnRR ~ Year, data = data)
summary(fityear)

time.plot <- data %>% 
  ggplot(aes(x = Year, y = lnRR))+
  geom_jitter(fill ='grey', alpha=.75, size = 4, shape = 16, colour ='#696969')+
  geom_hline(yintercept = 0, linetype = 2, colour = "steelblue", size = 0.75) +
  geom_smooth(method  = "lm", se = TRUE, color = "#696969", fill = "#696969", alpha   = 0.2) +
  guides(size = F) +
  labs(y = "Effect size (ln"~italic(RR)~")", x = 'Year of publication')+
  theme_bw()+
  theme(panel.spacing = unit(0.5, "lines"),
        panel.border= element_blank(),
        axis.line=element_line(),
        text = element_text(size = 12),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position    = "none",
        panel.grid.minor.x = element_blank(),
        legend.text = element_text(size = 12),
        legend.title=element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +
  geom_text(
    data = data.frame(x = 1985, y = 2),
    aes(x = x, y = y),
    label = "italic(R)^2 == '0.030' * ',' ~ italic(p) < 0.001",
    parse = TRUE,
    color = "#696969",
    size  = 4,
    inherit.aes = FALSE,
    hjust = 0) +
  geom_text(data = data.frame(x = 2030 , y = 2, label = "C"),
            aes(x = x, y = y, label = label),
            color = "black", fontface = "bold", size = 6,
            inherit.aes = FALSE) +
  scale_y_continuous(labels = number_format(accuracy = 0.1))

time.plot
ggsave("Figure 4c.pdf", IF.plot, height = 125, width = 100, units = "mm")

p <- grid.arrange(funnel.plot, IF.plot, time.plot, ncol = 1)
p
ggsave("Figure 5.pdf", p, height = 325, width = 125, units = "mm")









rm1
rm2
rm3
rm4
rm5
rm6

rm7
rm8

rm9
rm10

rm11
rm12

rm13
rm14

rm15
rm16

rm17
rm18
