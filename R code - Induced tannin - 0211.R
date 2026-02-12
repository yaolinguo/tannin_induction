##############################################
##### Sun & Guo, 2026 - Tannin induction #####
##############################################

library(nortest)
library(sf)
library(ggplot2) 
library(knitr)
library(tidyverse)
library(brms)
library(ape)
library(readxl)
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
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggtree)
library(V.PhyloMaker2)
library(grid)

# ==================================================================
# Figure 1a: Circular phylogeny with tip size mapped from lnRR_wmean
# ==================================================================

setwd("/Users/yaolin/Desktop/My papers/Manuscripts/2025 - Induced tannin/Submission - Plant Diversity - 0209/Data & Code")
data <- read_excel("Raw data - Induced tannin - 0211.xlsx")
write.csv(data, file = "data.csv", row.names = FALSE)
data <- read.csv("data.csv")

setwd("/Users/yaolin/Desktop/My papers/Manuscripts/2025 - Induced tannin/Submission - Plant Diversity - 0209/Data & Code")
infile  <- "Raw data - Induced tannin - 0211.xlsx"
out_pdf <- "Figure 1a.pdf"

Mode1 <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}

wmean_safe <- function(x, rvar) {
  w <- 1 / rvar
  w[!is.finite(w) | w <= 0] <- NA_real_
  if (all(is.na(w))) return(mean(x, na.rm = TRUE))
  weighted.mean(x, w = w, na.rm = TRUE)
}

# Read + clean data
dat0 <- readxl::read_excel(infile) %>% as.data.frame()
needed <- c("lnRR", "rvar", "Family", "Plant_species", "Plant_type")
miss <- setdiff(needed, names(dat0))
if (length(miss) > 0) stop("Missing columns: ", paste(miss, collapse = ", "))

dat <- dat0 %>%
  transmute(lnRR = suppressWarnings(as.numeric(lnRR)),
            rvar = suppressWarnings(as.numeric(rvar)),
            Family = str_squish(as.character(Family)),
            Plant_species = str_squish(str_replace_all(as.character(Plant_species), "\ufeff", "")),
            Plant_type = str_squish(as.character(Plant_type))) %>%
  separate(Plant_species, into = c("genus", "species_ep"), sep = "\\s+", remove = FALSE,
           extra = "drop", fill = "right") %>%
  filter(!is.na(genus), !is.na(species_ep)) %>%
  filter(!tolower(species_ep) %in% c("sp", "sp.", "spp", "spp.")) %>%
  mutate(
    genus = str_to_title(genus),
    species_ep = str_to_lower(species_ep),
    Plant_species_binom = paste(genus, species_ep),                      # "Genus species"
    label = str_to_lower(paste(genus, species_ep, sep = "_")),           # "genus_species"
    Plant_species_display = paste(genus, species_ep),
    Plant_type = case_when(
      str_to_lower(Plant_type) %in% c("woody") ~ "Woody",
      str_to_lower(Plant_type) %in% c("non-woody", "nonwoody", "non woody", "herbaceous") ~ "Non-woody",
      TRUE ~ Plant_type
    )
  ) %>%
  filter(is.finite(lnRR))

# Species-level aggregation (inverse-variance weighted mean)
data_tree <- dat %>%
  group_by(label) %>%
  summarise(genus = first(genus),
            species_ep = first(species_ep),
            Plant_species_binom = first(Plant_species_binom),
            Plant_species = first(Plant_species_display),
            Family = Mode1(Family),
            Plant_type = Mode1(Plant_type),
            n_es = n(),
            lnRR_mean  = mean(lnRR, na.rm = TRUE),
            lnRR_wmean = wmean_safe(lnRR, rvar),
            .groups = "drop")
stopifnot(!anyDuplicated(data_tree$label))

# EXACT & stable size mapping
lnRR_min <- -0.5
lnRR_max <-  1.0
size_min <-  1.0  # mm
size_max <-  6.0  # mm

size_from_lnRR <- function(x) {
  x2 <- pmin(pmax(x, lnRR_min), lnRR_max)  # clamp
  size_min + (size_max - size_min) * (x2 - lnRR_min) / (lnRR_max - lnRR_min)
}

data_tree <- data_tree %>%
  mutate(
    pt_size = size_from_lnRR(lnRR_wmean)   # <- this is the ACTUAL point size in mm
  )

# Build phylogeny
sp_list <- data_tree %>%
  distinct(genus, species_ep, Family, Plant_species_binom) %>%
  filter(!is.na(Family), Family != "") %>%
  transmute(species = Plant_species_binom,  # "Genus species"
            genus   = genus,
            family  = Family)
phylo_res <- V.PhyloMaker2::phylo.maker(sp.list = sp_list)
tree_data <- phylo_res$scenario.3

tree_data$tip.label <- str_to_lower(str_replace_all(tree_data$tip.label, " ", "_"))
stopifnot(!anyDuplicated(tree_data$tip.label))

data_attach <- data_tree %>% filter(label %in% tree_data$tip.label)
stopifnot(!anyDuplicated(data_attach$label))

# Plot
legend_labs <- c(-0.5, 0, 0.5, 1)
legend_brks <- size_from_lnRR(legend_labs)  # breaks must be in mm (pt_size scale)
p <- ggtree(tree_data, layout = "circular") %<+% data_attach +
  geom_tree() +
  geom_tippoint(aes(color = Plant_type, size = pt_size)) +
  geom_tiplab(aes(label = Plant_species), size = 3, hjust = -0.1, offset = 1) +
  scale_color_manual(values = c("Woody" = "#f16c23", "Non-woody" = "#2b6a99")) +
  scale_size_identity(name   = "lnRR (inv-var wmean; more negative = smaller)",
                      breaks = legend_brks,
                      labels = legend_labs,
                      guide  = "legend") +
  theme(plot.margin = unit(c(3, 3, 3, 3), "cm"),
        text = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = c(2.2, 0.2),
        legend.justification = c(2.2, 0.2))
print(p)
ggsave(out_pdf, p, height = 150, width = 300, units = "mm")
cat("Saved:", out_pdf, "\n")

# =======================================================================
# Figure 1b: overall effect size and contrast between woody and non-woody
# =======================================================================

setwd("/Users/yaolin/Desktop/My papers/Manuscripts/2025 - Induced tannin/Submission - Plant Diversity - 0209/Data & Code")
data <- read_excel("Raw data - Induced tannin - 0211.xlsx")
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

p <- ggplot(dat, aes(x = coef, y = Type)) +
  geom_errorbarh(mapping = aes(xmin = lower, xmax = upper), 
                 width = 0.075, size = 1.25, colour = "grey10", height = 0.2) +
  geom_point(size = 2.75, shape = 21, color = "grey10", stroke = 1.25, 
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
  scale_fill_manual(values = c("Total"     = "#bfbfbf",
                               "Woody"     = "#f16c23",
                               "Non-woody" = "#2b6a99"),
                    breaks = c("Non-woody", "Woody", "Total")) +
  scale_y_discrete(limits = c("Non-woody", "Woody", "Total")) +
  scale_x_continuous(limits = c(-0.6, 0.6),
                     labels = scales::number_format(accuracy = 0.1)) +
  xlab(expression(paste("Effect size (ln ", italic("RR"), ")"))) +
  ylab(NULL) +
  coord_fixed(ratio = 0.6)
p
ggsave("./Figure 1b.pdf", p, height = 125, width = 200, units = "mm")

# =====================================================================
# Figure 2: tannin group, leaf cohort, leaf damage status and life span
# =====================================================================

setwd("/Users/yaolin/Desktop/My papers/Manuscripts/2025 - Induced tannin/Submission - Plant Diversity - 0209/Data & Code")
data <- read_excel("Raw data - Induced tannin - 0211.xlsx")
write.csv(data, file = "data.csv", row.names = FALSE)
data <- read.csv("data.csv")

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
              random = list( ~ 1 | Publication/Effect_ID, ~ 1 | Plant_species),
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
                 width = 0.15, size = 1.25, colour = "grey10", height = 0.5) +
  geom_point(size = 2.75, shape = 21, color = "grey10", stroke = 1.25, 
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
                 width = 0.15, size = 1.25, colour = "grey10", height = 0.2) +
  geom_point(size = 2.75, shape = 21, color = "grey10", stroke = 1.25, 
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
                 width = 0.15, size = 1.25, colour = "grey10", height = 0.2) +
  geom_point(size = 2.75, shape = 21, color = "grey10", stroke = 1.25, 
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

# Figure S1
setwd("/Users/yaolin/Desktop/My papers/Manuscripts/2025 - Induced tannin/Submission - Plant Diversity - 0209/Data & Code")
data <- read_excel("Raw data - Induced tannin - 0211.xlsx")

data$Plant_species <- gsub("\ufeff", "", data$Plant_species)
data$Life_span     <- trimws(data$Life_span)
data$Plant_type    <- trimws(data$Plant_type)
data$Life_span <- ifelse(grepl("Annual", data$Life_span, ignore.case = TRUE), "Annual", "Perennial")
data$Life_span <- factor(data$Life_span, levels = c("Annual", "Perennial"))

rm9 <- rma.mv(lnRR, rvar,
                 mods = ~ 0 + factor(Life_span),
                 random = list(~ 1 | Publication/Effect_ID, ~ 1 | Plant_species),
                 method = "REML",
                 data = data)
rm9
coef(rm9)

L <- rbind("Perennial − Annual" = c(-1, 1))
anova(rm9, L = L)

rm10 <- rma.mv(lnRR, rvar,
                 mods = ~ 0 + paste(Life_span, Plant_type),
                 random = list(~ 1 | Publication/Effect_ID, ~ 1 | Plant_species),
                 method = "REML",
                 data = data)
rm10

bn <- names(coef(rm10))
p  <- length(bn)
ix <- function(pat) which(grepl(pat, bn))
mk <- function(plus, minus){
  v <- rep(0, p); v[plus] <- 1; v[minus] <- -1; v
}

i_AN <- ix("^paste\\(Life_span, Plant_type\\)Annual.*Non-woody|Annual.*Non-woody")
i_AW <- ix("^paste\\(Life_span, Plant_type\\)Annual.*Woody|Annual.*Woody")
i_PN <- ix("^paste\\(Life_span, Plant_type\\)Perennial.*Non-woody|Perennial.*Non-woody")
i_PW <- ix("^paste\\(Life_span, Plant_type\\)Perennial.*Woody|Perennial.*Woody")

L2 <- rbind(
  if(length(i_PN)==1 && length(i_AN)==1) "Perennial − Annual | Non-woody" = mk(i_PN, i_AN),
  if(length(i_PW)==1 && length(i_AW)==1) "Perennial − Annual | Woody"     = mk(i_PW, i_AW),
  if(length(i_AW)==1 && length(i_AN)==1) "Woody − Non-woody | Annual"     = mk(i_AW, i_AN),
  if(length(i_PW)==1 && length(i_PN)==1) "Woody − Non-woody | Perennial"  = mk(i_PW, i_PN)
)
anova(rm10, L = L2)

as.data.frame(table(data$Life_span))
as.data.frame(table(data$Plant_type, data$Life_span))

dat1 <- data.frame(coef = coef(rm9),
                   se   = sqrt(diag(vcov(rm9))),
                   Type = gsub("^factor\\(Life_span\\)", "", names(coef(rm9))))

dat1$Type  <- ifelse(dat1$Type == "Annual", "Annual", "Perennial")
dat1$lower <- dat1$coef - 1.96 * dat1$se
dat1$upper <- dat1$coef + 1.96 * dat1$se

bn2 <- names(coef(rm10))
Type2 <- gsub("^paste\\(Life_span, Plant_type\\)", "", bn2)
Type2 <- trimws(Type2)
Type2 <- gsub("  ", " ", Type2)
Type2 <- gsub("Annual Non-woody",    "Annual - Non-woody",    Type2)
Type2 <- gsub("Annual Woody",        "Annual - Woody",        Type2)
Type2 <- gsub("Perennial Non-woody", "Perennial - Non-woody", Type2)
Type2 <- gsub("Perennial Woody",     "Perennial - Woody",     Type2)

dat2 <- data.frame(coef = coef(rm10),
                   se   = sqrt(diag(vcov(rm10))),
                   Type = Type2)
  
dat2$lower <- dat2$coef - 1.96 * dat2$se
dat2$upper <- dat2$coef + 1.96 * dat2$se

dat <- rbind(dat2, dat1)
fill_map <- c("Annual - Non-woody"    = "#2b6a99",
              "Perennial - Non-woody" = "#2b6a99",
              "Annual - Woody"        = "#f16c23",
              "Perennial - Woody"     = "#f16c23",
              "Annual"                = "#bfbfbf",
              "Perennial"             = "#bfbfbf")
y_order <- c("Annual - Non-woody", "Perennial - Non-woody",
             "Annual - Woody", "Perennial - Woody",
             "Annual", "Perennial")
y_order <- y_order[y_order %in% dat$Type]

p4 <- ggplot(dat, aes(x = coef, y = Type)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper),
                 width = 0.15, size = 1.25, colour = "grey10", height = 0.5) +
  geom_point(size = 2.75, shape = 21, color = "grey10", stroke = 1.25,
             aes(fill = Type)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "steelblue", size = 1) +
  theme_classic() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size = 12, hjust = 0.55),
        axis.text.y = element_text(size = 12, hjust = 1),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  scale_fill_manual(values = fill_map, breaks = y_order) +
  scale_y_discrete(limits = y_order) +
  xlim(c(-1, 1)) +
  xlab("Effect size (ln"~italic(RR)~")") +
  ylab("Life span and plant growth form") +
  coord_fixed(ratio = 0.64)

p4
ggsave("./Figure S1.pdf", p_ls, height = 125, width = 200, units = "mm")

p4 <- p1 + p2 + p3
p4
ggsave('./Figure 2 - 0829.pdf', p4, height = 250, width = 600, units = c("mm"))

# ==========================
# Figure 3: induction regime
# ==========================

# Figure 3a #
setwd("/Users/yaolin/Desktop/My papers/Manuscripts/2025 - Induced tannin/Submission - Plant Diversity - 0209/Data & Code")
data <- read_excel("Raw data - Induced tannin - 0211.xlsx")
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
  theme(panel.grid.major.x = element_blank(),
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
        axis.text          = element_text(size = 12)) +
  xlab("Herbivory intensity (%)") +
  ylab("Effect size (ln"~italic(RR)~")") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  geom_text(data = data.frame(x = 4, y = -1.25),
            aes(x = x, y = y),
            label = "italic(R)[Total]^2 == '0.137' * ',' ~ italic(p) < 0.001",
            parse = TRUE,
            color = "#696969",
            size  = 4,
            inherit.aes = FALSE,
            hjust = 0) +
  geom_text(data = data.frame(x = 4, y = -1.5),
            aes(x = x, y = y),
            label = "italic(R)[Woody]^2 == '0.274' * ',' ~ italic(p) < 0.001",
            parse = TRUE,
            color = "#f16c23",
            size  = 4,
            inherit.aes = FALSE,
            hjust = 0) +
  geom_text(data = data.frame(x = 4, y = -1.75),
            aes(x = x, y = y),
            label = "italic(R)[Non-woody]^2 == '0.001' * ',' ~ italic(p) == 0.778",
            parse = TRUE,
            color = "#2b6a99",
            size  = 4,
            inherit.aes = FALSE,
            hjust = 0)
p1
p1 <- ggMarginal(p1,
                 type        = "density",
                 margins     = "both",
                 groupColour = TRUE,
                 groupFill   = TRUE,
                 alpha       = 0.4)
print(p1)
ggsave("Figure 3a.pdf", p1, height = 150, width = 400, units = "mm")



#Figure 3b
setwd("/Users/yaolin/Desktop/My papers/Manuscripts/2025 - Induced tannin/Submission - Plant Diversity - 0209/Data & Code")
data <- read_excel("Raw data - Induced tannin - 0211.xlsx")
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
  theme(panel.grid.major.x = element_blank(),
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
        axis.text          = element_text(size = 12)) +
  xlab("Treatment interval (days)") +
  ylab("Effect size (ln"~italic(RR)~")") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  geom_text(data = data.frame(x = 4, y = -1.25),
            aes(x = x, y = y),
            label = "italic(R)[Total]^2 == '0.002' * ',' ~ italic(p) == 0.382",
            parse = TRUE,
            color = "#696969",
            size  = 4,
            inherit.aes = FALSE,
            hjust = 0) +
  geom_text(data = data.frame(x = 4, y = -1.5),
            aes(x = x, y = y),
            label = "italic(R)[Non-woody]^2 == '0.000' * ',' ~ italic(p) == 0.834",
            parse = TRUE,
            color = "#f16c23",
            size  = 4,
            inherit.aes = FALSE,
            hjust = 0) +
  geom_text(data = data.frame(x = 4, y = -1.75),
            aes(x = x, y = y),
            label = "italic(R)[Non-woody]^2 == '0.078' * ',' ~ italic(p) == 0.004",
            parse = TRUE,
            color = "#2b6a99",
            size  = 4,
            inherit.aes = FALSE,
            hjust = 0)
p2
p2 <- ggMarginal(p2,
                 type        = "density",
                 margins     = "both",
                 groupColour = TRUE,
                 groupFill   = TRUE,
                 alpha       = 0.4)
print(p2)
ggsave("Figure 3b.pdf", p2, height = 150, width = 400, units = "mm")



#Figure 3c
setwd("/Users/yaolin/Desktop/My papers/Manuscripts/2025 - Induced tannin/Submission - Plant Diversity - 0209/Data & Code")
data <- read_excel("Raw data - Induced tannin - 0211.xlsx")
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
  theme(panel.grid.major.x = element_blank(),
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
        axis.text          = element_text(size = 12)) +
  xlab("Measurement interval (days)") +
  ylab("Effect size (ln"~italic(RR)~")") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  geom_text(data = data.frame(x = 4, y = -1.25),
            aes(x = x, y = y),
            label = "italic(R)[Total]^2 == '0.002' * ',' ~ italic(p) == 0.384",
            parse = TRUE,
            color = "#696969",
            size  = 4,
            inherit.aes = FALSE,
            hjust = 0) +
  geom_text(data = data.frame(x = 4, y = -1.5),
            aes(x = x, y = y),
            label = "italic(R)[Non-woody]^2 == '0.001' * ',' ~ italic(p) == 0.700",
            parse = TRUE,
            color = "#f16c23",
            size  = 4,
            inherit.aes = FALSE,
            hjust = 0) +
  geom_text(data = data.frame(x = 4, y = -1.75),
            aes(x = x, y = y),
            label = "italic(R)[Non-woody]^2 == '0.009' * ',' ~ italic(p) == 0.315",
            parse = TRUE,
            color = "#2b6a99",
            size  = 4,
            inherit.aes = FALSE,
            hjust = 0)
p3
p3 <- ggMarginal(p3,
                 type        = "density",
                 margins     = "both",
                 groupColour = TRUE,
                 groupFill   = TRUE,
                 alpha       = 0.4)
print(p3)
ggsave("Figure 3c.pdf", p3, height = 150, width = 400, units = "mm")

p4 <- grid.arrange(p1, p2, p3, ncol = 1)
p4
ggsave('./Figure 3 - 1226.pdf', p4, height = 400, width = 150, units = c("mm"))

# ==========================
# Figure 4: induction method
# ==========================

# Figure 4a #
setwd("/Users/yaolin/Desktop/My papers/Manuscripts/2025 - Induced tannin/Submission - Plant Diversity - 0209/Data & Code")
data <- read_excel("Raw data - Induced tannin - 0211.xlsx")
write.csv(data, file = "data.csv", row.names = FALSE)
data <- read.csv("data.csv")

rm11 <- rma.mv(lnRR, rvar,
              mods = ~ paste(Treatment_type) - 1,
              random = list( ~ 1 | Publication/Effect_ID, ~ 1 | Plant_species),
              method = "REML", 
              data = data)
rm11

bn <- names(coef(rm11)); p <- length(bn)
i_A <- which(bn == "paste(Treatment_type)Artificial")
i_C <- which(bn == "paste(Treatment_type)Chemical")
i_H <- which(bn == "paste(Treatment_type)Herbivore")
i_J <- which(bn == "paste(Treatment_type)Joint")
mk <- function(plus, minus){ v <- numeric(p); v[plus] <- 1; v[minus] <- -1; v }
L <- rbind("Artificial − Chemical"  = mk(i_A, i_C),
           "Artificial − Herbivore" = mk(i_A, i_H),
           "Artificial − Joint"     = mk(i_A, i_J),
           "Herbivore − Chemical"  = mk(i_H, i_C),
           "Herbivore − Joint"     = mk(i_H, i_J),
           "Joint − Chemical"      = mk(i_J, i_C))
anova(rm11, L = L)

rm12 <- rma.mv(lnRR, rvar,
               mods = ~ paste(Treatment_type, Plant_type) - 1,
               random = list( ~ 1 | Publication/Effect_ID, ~ 1 | Plant_species),
               method = "REML", 
               data = data)
rm12

bn <- names(coef(rm12)); p <- length(bn)
pick <- function(method, form) {
  which(bn == paste0("paste(Treatment_type, Plant_type)", method, " ", form))
}

i_AN <- pick("Artificial", "Non-woody")
i_AW <- pick("Artificial", "Woody")
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
anova(rm12, L = L)

sample.sizes <- as.data.frame(table(data$Treatment_type))
sample.sizes
sample.sizes <- as.data.frame(table(data$Plant_type, data$Treatment_type))
sample.sizes

dat1 <- data.frame(coef = coef(rm11), se = sqrt(diag(vcov(rm11))), Type = c("Artificial treatment", 
                                                                            "Chemical treatment",
                                                                            "Herbivore treatment", 
                                                                            "Joint treatment"))
dat1$lower <- dat1$coef - 1.96 * dat1$se
dat1$upper <- dat1$coef + 1.96 * dat1$se
dat1
dat2 <- data.frame(coef = coef(rm12), se = sqrt(diag(vcov(rm12))), Type = c("Artificial treatment - Non-woody", 
                                                                            "Artificial treatment - Woody",
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
                 width = 0.15, size = 1.25, colour = "grey10", height = 0.2) +
  geom_point(size = 2.75, shape = 21, color = "grey10", stroke = 1.25, 
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
                               "Artificial treatment - Non-woody" = "#2b6a99",
                               "Joint treatment - Woody" = "#f16c23",
                               "Chemical treatment - Woody" = "#f16c23",
                               "Herbivore treatment - Woody" = "#f16c23",
                               "Artificial treatment - Woody" = "#f16c23",
                               "Joint treatment" = "#bfbfbf",
                               "Chemical treatment" = "#bfbfbf",
                               "Herbivore treatment" = "#bfbfbf",
                               "Artificial treatment" = "#bfbfbf"),
                    breaks = c("Joint treatment - Non-woody",
                               "Chemical treatment - Non-woody",
                               "Herbivore treatment - Non-woody",
                               "Artificial treatment - Non-woody",
                               "Joint treatment - Woody",
                               "Chemical treatment - Woody",
                               "Herbivore treatment - Woody",
                               "Artificial treatment - Woody",
                               "Joint treatment",
                               "Chemical treatment",
                               "Herbivore treatment",
                               "Artificial treatment")) +
  scale_y_discrete(limits = c("Joint treatment - Non-woody",
                              "Chemical treatment - Non-woody",
                              "Herbivore treatment - Non-woody",
                              "Artificial treatment - Non-woody",
                              "Joint treatment - Woody",
                              "Chemical treatment - Woody",
                              "Herbivore treatment - Woody",
                              "Artificial treatment - Woody",
                              "Joint treatment",
                              "Chemical treatment",
                              "Herbivore treatment",
                              "Artificial treatment")) +
  xlim(c(-1, 1)) +
  xlab("Effect size (ln"~italic(RR)~")") +
  ylab("Induction treatment and plant growth form") +
  coord_fixed(ratio = 0.5)
p1
ggsave('./Figure 4a.pdf', p1, height = 260, width = 175, units = c("mm"))



# Figure 4b #
setwd("/Users/yaolin/Desktop/My papers/Manuscripts/2025 - Induced tannin/Submission - Plant Diversity - 0209/Data & Code")
data <- read_excel("Raw data - Induced tannin - 0211.xlsx")
write.csv(data, file = "data.csv", row.names = FALSE)
data <- read.csv("data.csv")

data1 <- data %>%
  filter(Treatment_type == "Artificial") %>%
  droplevels()                              

rm13 <- rma.mv(lnRR, rvar,
               mods = ~ paste(Treatment_typeII) - 1,
               random = list( ~ 1 | Publication/Effect_ID, ~ 1 | Plant_species),
               method = "REML", 
               data = data1)
rm13

bn <- names(coef(rm13)); p <- length(bn)
i_D <- which(bn == "paste(Treatment_typeII)Defoliation")
i_Pr <- which(bn == "paste(Treatment_typeII)Pruning")
i_Pu <- which(bn == "paste(Treatment_typeII)Punching")
mk <- function(plus, minus){ v <- numeric(p); v[plus] <- 1; v[minus] <- -1; v }
L <- rbind("Defoliation − Pruning"   = mk(i_D,  i_Pr),
           "Defoliation − Punching"  = mk(i_D,  i_Pu),
           "Pruning − Punching"      = mk(i_Pr, i_Pu))
anova(rm13, L = L)

rm14 <- rma.mv(lnRR, rvar,
               mods = ~ paste(Treatment_typeII, Plant_type) - 1,
               random = list( ~ 1 | Publication/Effect_ID, ~ 1 | Plant_species),
               method = "REML", 
               data = data1)
rm14

bn <- names(coef(rm14)); p <- length(bn)
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
anova(rm14, L = L)

sample.sizes <- as.data.frame(table(data1$Treatment_typeII))
sample.sizes
sample.sizes <- as.data.frame(table(data1$Plant_type, data1$Treatment_typeII))
sample.sizes

dat1 <- data.frame(coef = coef(rm13), se = sqrt(diag(vcov(rm13))), Type = c("Defoliation", 
                                                                            "Pruning",
                                                                            "Punching"))
dat1$lower <- dat1$coef - 1.96 * dat1$se
dat1$upper <- dat1$coef + 1.96 * dat1$se
dat1

dat2 <- data.frame(coef = coef(rm14), se = sqrt(diag(vcov(rm14))), Type = c("Defoliation - Non-woody", 
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
                 width = 0.15, size = 1.25, colour = "grey10", height = 0.2) +
  geom_point(size = 2.75, shape = 21, color = "grey10", stroke = 1.25, 
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
ggsave('./Figure 4b.pdf', p2, height = 125, width = 200, units = c("mm"))

# Figure 4c
setwd("/Users/yaolin/Desktop/My papers/Manuscripts/2025 - Induced tannin/Submission - Plant Diversity - 0209/Data & Code")
data <- read_excel("Raw data - Induced tannin - 0211.xlsx")
write.csv(data, file = "data.csv", row.names = FALSE)
data <- read.csv("data.csv")
data1 <- data %>%
  filter(Treatment_type == "Herbivore") %>%
  droplevels()                              

rm15 <- rma.mv(lnRR, rvar,
               mods = ~ paste(Treatment_typeII) - 1,
               random = list( ~ 1 | Publication/Effect_ID, ~ 1 | Plant_species),
               method = "REML", 
               data = data1)
rm15

bn <- names(coef(rm15))
i_G <- which(bn == "paste(Treatment_typeII)Generalist")
i_S <- which(bn == "paste(Treatment_typeII)Specialist")
L <- rbind("Specialist − Generalist" = {v <- numeric(length(bn)); v[i_S] <- 1; v[i_G] <- -1; v})
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
i_GN <- pick("Generalist", "Non-woody")
i_GW <- pick("Generalist", "Woody")
i_SN <- pick("Specialist", "Non-woody")
i_SW <- pick("Specialist", "Woody")
mk <- function(plus, minus){ v <- numeric(p); v[plus] <- 1; v[minus] <- -1; v }
L <- rbind("Woody − Non-woody | Generalist" = mk(i_GW, i_GN),
           "Woody − Non-woody | Specialist" = mk(i_SW, i_SN))
anova(rm16, L = L)

sample.sizes <- as.data.frame(table(data1$Treatment_typeII))
sample.sizes
sample.sizes <- as.data.frame(table(data1$Plant_type, data1$Treatment_typeII))
sample.sizes

dat1 <- data.frame(coef = coef(rm15), se = sqrt(diag(vcov(rm15))), Type = c("Generalist", "Specialist"))
dat1$lower <- dat1$coef - 1.96 * dat1$se
dat1$upper <- dat1$coef + 1.96 * dat1$se
dat1
dat2 <- data.frame(coef = coef(rm16), se = sqrt(diag(vcov(rm16))), Type = c("Generalist - Non-woody", 
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
                 width = 0.15, size = 1.25, colour = "grey10", height = 0.2) +
  geom_point(size = 2.75, shape = 21, color = "grey10", stroke = 1.25, 
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
ggsave('./Figure 4c.pdf', p3, height = 125, width = 200, units = c("mm"))

# Figure 4d
setwd("/Users/yaolin/Desktop/My papers/Manuscripts/2025 - Induced tannin/Submission - Plant Diversity - 0209/Data & Code")
data <- read_excel("Raw data - Induced tannin - 0211.xlsx")
write.csv(data, file = "data.csv", row.names = FALSE)
data <- read.csv("data.csv")

data1 <- data %>%
  filter(Treatment_type == "Chemical") %>%
  droplevels()                              

rm17 <- rma.mv(lnRR, rvar,
               mods = ~ paste(Treatment_typeII) - 1,
               random = list( ~ 1 | Publication/Effect_ID, ~ 1 | Plant_species),
               method = "REML", 
               data = data1)
rm17

bn <- names(coef(rm17))
i_JA <- which(bn == "paste(Treatment_typeII)JA")
i_SA <- which(bn == "paste(Treatment_typeII)SA")
L <- rbind("JA − SA" = {v <- numeric(length(bn)); v[i_JA] <- 1; v[i_SA] <- -1; v})
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
i_JN <- pick("JA", "Non-woody")
i_JW <- pick("JA", "Woody")
i_SN <- pick("SA", "Non-woody")
i_SW <- pick("SA", "Woody")
mk <- function(plus, minus){ v <- numeric(p); v[plus] <- 1; v[minus] <- -1; v }
L <- rbind("Woody − Non-woody | JA" = mk(i_JW, i_JN),
           "Woody − Non-woody | SA" = mk(i_SW, i_SN))
anova(rm18, L = L)

sample.sizes <- as.data.frame(table(data1$Treatment_typeII))
sample.sizes
sample.sizes <- as.data.frame(table(data1$Plant_type, data1$Treatment_typeII))
sample.sizes

dat1 <- data.frame(coef = coef(rm17), se = sqrt(diag(vcov(rm17))), Type = c("JA", "SA"))
dat1$lower <- dat1$coef - 1.96 * dat1$se
dat1$upper <- dat1$coef + 1.96 * dat1$se
dat1
dat2 <- data.frame(coef = coef(rm18), se = sqrt(diag(vcov(rm18))), Type = c("JA - Non-woody", 
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
                 width = 0.15, size = 1.25, colour = "grey10", height = 0.2) +
  geom_point(size = 2.75, shape = 21, color = "grey10", stroke = 1.25, 
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
ggsave('./Figure 4d.pdf', p4, height = 125, width = 200, units = c("mm"))

# Figure 4e
setwd("/Users/yaolin/Desktop/My papers/Manuscripts/2025 - Induced tannin/Submission - Plant Diversity - 0209/Data & Code")
data <- read_excel("Raw data - Induced tannin - 0211.xlsx")
write.csv(data, file = "data.csv", row.names = FALSE)
data <- read.csv("data.csv")

data1 <- data %>%
  filter(Treatment_type == "Joint") %>%
  droplevels()                              

rm19 <- rma.mv(lnRR, rvar,
               mods = ~ paste(Treatment_typeII) - 1,
               random = list( ~ 1 | Publication/Effect_ID, ~ 1 | Plant_species),
               method = "REML", 
               data = data1)
rm19

bn <- names(coef(rm19)); p <- length(bn)
w  <- function(label) which(bn == paste0("paste(Treatment_typeII)", label))
i_DH <- w("Defoliation + Herbivore")
i_DJ <- w("Defoliation + JA")
i_DS <- w("Defoliation + Salivation")
mk <- function(plus, minus){ v <- numeric(p); v[plus] <- 1; v[minus] <- -1; v }
L <- rbind("Defol+Herb − Defol+JA"        = mk(i_DH, i_DJ),
           "Defol+Herb − Defol+Salivation"= mk(i_DH, i_DS),
           "Defol+JA   − Defol+Salivation"= mk(i_DJ, i_DS))
anova(rm19, L = L)

rm20 <- rma.mv(lnRR, rvar,
               mods = ~ paste(Treatment_typeII, Plant_type) - 1,
               random = list( ~ 1 | Publication/Effect_ID, ~ 1 | Plant_species),
               method = "REML", 
               data = data1)
rm20

bn <- names(coef(rm20)); p <- length(bn)
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
anova(rm20, L = L)

sample.sizes <- as.data.frame(table(data1$Treatment_typeII))
sample.sizes
sample.sizes <- as.data.frame(table(data1$Plant_type, data1$Treatment_typeII))
sample.sizes

dat1 <- data.frame(coef = coef(rm19), se = sqrt(diag(vcov(rm19))), Type = c("Defoliation + Herbivore", 
                                                                            "Defoliation + JA",
                                                                            "Defoliation + Salivation"))
dat1$lower <- dat1$coef - 1.96 * dat1$se
dat1$upper <- dat1$coef + 1.96 * dat1$se
dat1
dat2 <- data.frame(coef = coef(rm20), se = sqrt(diag(vcov(rm20))), Type = c("Defoliation + Herbivore Woody", 
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
                 width = 0.15, size = 1.25, colour = "grey10", height = 0.2) +
  geom_point(size = 2.75, shape = 21, color = "grey10", stroke = 1.25, 
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
ggsave('./Figure 4e.pdf', p5, height = 125, width = 200, units = c("mm"))

p6 <- p2 + p3 + p4 + p5
p6
ggsave('./Figure 4.pdf', p6, height = 250, width = 400, units = c("mm"))

# =======================
# Figure 5: Test the bias
# =======================

# Figure 5a
setwd("/Users/yaolin/Desktop/My papers/Manuscripts/2025 - Induced tannin/Submission - Plant Diversity - 0209/Data & Code")
data <- read_excel("Raw data - Induced tannin - 0211.xlsx")
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
  se.seq <- seq(from = 0, 
                to   = max(sqrt(dataset$rvar), na.rm = TRUE),
                by   = 0.001)
  dfCI <- data.frame(ll95     = estimate - (1.96 * se.seq),
                     ul95     = estimate + (1.96 * se.seq),
                     ll99     = estimate - (3.29 * se.seq),
                     ul99     = estimate + (3.29 * se.seq),
                     se.seq   = se.seq,
                     meanll95 = estimate - (1.96 * SE),
                     meanul95 = estimate + (1.96 * SE))
  ggplot(dataset, aes(x = sqrt(rvar), y = lnRR)) +
    geom_point(color = '#696969', shape = 16, size = 4, alpha = 0.6) +
    xlab("Standard error") +
    ylab("Effect size (ln"~italic(RR)~")") +
    geom_line(aes(x = se.seq, y = ll95), data = dfCI, linetype = 'dotted') +
    geom_line(aes(x = se.seq, y = ul95), data = dfCI, linetype = 'dotted') +
    geom_line(aes(x = se.seq, y = ll99), data = dfCI, linetype = 'dashed') +
    geom_line(aes(x = se.seq, y = ul99), data = dfCI, linetype = 'dashed') +
    geom_segment(aes(x = min(se.seq), y = meanll95, xend = max(se.seq), yend = meanll95),
                 data    = dfCI,
                 colour  = "steelblue",
                 linetype= 'dotdash',
                 size    = 0.75) +
    geom_segment(aes(x = min(se.seq), y = meanul95, xend = max(se.seq), yend = meanul95),
                 data    = dfCI,
                 colour  = "steelblue",
                 linetype= 'dotdash',
                 size    = 0.75) +
    scale_x_reverse() +
    coord_flip() +
    theme_bw() +
    theme(panel.spacing      = unit(0.5, "lines"),
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
          axis.text.y = element_text(size = 12)) +
    geom_text(data = data.frame(x = 0 , y = -4.5),
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
  geom_smooth(method = "lm", se = TRUE,
              color = "#696969", fill = "#696969", alpha = 0.2, linetype = "dashed") +
  labs(y = "Effect size (ln"~italic(RR)~")", x = 'Journal impact factor') +
  scale_x_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 10)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme_bw() +
  theme(panel.spacing = unit(0.5, "lines"),
        panel.border= element_blank(),
        text = element_text(size = 12),
        axis.line=element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  geom_text(data = data.frame(x = 0, y = 2),
            aes(x = x, y = y),
            label = "italic(R)^2 == '0.004' * ',' ~ italic(p) == 0.151",
            parse = TRUE,
            color = "#696969", size = 4, inherit.aes = FALSE, hjust = 0) +
  geom_text(data = data.frame(x = 11, y = 2, label = "B"),
            aes(x = x, y = y, label = label),
            color = "black", fontface = "bold", size = 6, inherit.aes = FALSE)
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
  geom_text(data = data.frame(x = 1985, y = 2),
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
ggsave("Figure 5c.pdf", IF.plot, height = 125, width = 100, units = "mm")

p <- grid.arrange(funnel.plot, IF.plot, time.plot, ncol = 1)
p
ggsave("Figure 5.pdf", p, height = 325, width = 125, units = "mm")

#===================================
## Table S generator (rm3 ... rm20)
## Wald contrasts: plus - minus
#===================================

library(metafor)
library(dplyr)
library(tibble)

.pick_unique <- function(x, pattern, ignore.case = FALSE) {
  hit <- grep(pattern, x, value = TRUE, ignore.case = ignore.case, perl = TRUE)
  if (length(hit) == 1) return(hit)
  return(NA_character_)  # 0 or >1 => NA (safer)
}

.make_L_by_name <- function(fit, plus_name, minus_name) {
  bnames <- names(coef(fit))
  if (is.na(plus_name) || is.na(minus_name)) return(NA)
  
  if (!(plus_name %in% bnames) || !(minus_name %in% bnames)) return(NA)
  
  L <- rep(0, length(bnames)); names(L) <- bnames
  L[plus_name]  <-  1
  L[minus_name] <- -1
  L
}

.make_L_by_pattern <- function(fit, plus_pat, minus_pat, ignore.case = FALSE) {
  bnames <- names(coef(fit))
  plus  <- .pick_unique(bnames, plus_pat, ignore.case = ignore.case)
  minus <- .pick_unique(bnames, minus_pat, ignore.case = ignore.case)
  .make_L_by_name(fit, plus, minus)
}

.wald_from_L <- function(fit, L) {
  if (all(is.na(L))) {
    return(tibble(Estimate=NA_real_, SE=NA_real_, z=NA_real_, p=NA_real_,
                  CI_low=NA_real_, CI_high=NA_real_))
  }
  b <- coef(fit)
  V <- vcov(fit)
  
  est <- as.numeric(sum(L * b))
  se  <- sqrt(as.numeric(t(L) %*% V %*% L))
  z   <- est / se
  p   <- 2 * pnorm(-abs(z))
  ci  <- est + c(-1, 1) * qnorm(0.975) * se
  
  tibble(Estimate=est, SE=se, z=z, p=p, CI_low=ci[1], CI_high=ci[2])
}

.run_specs <- function(fit, model_label, specs) {
  out <- lapply(seq_len(nrow(specs)), function(i){
    L   <- .make_L_by_pattern(fit, specs$plus_pat[i], specs$minus_pat[i])
    res <- .wald_from_L(fit, L)
    res$Model <- model_label
    res$Test  <- specs$Test[i]
    res
  })
  bind_rows(out) %>% select(Model, Test, everything())
}

# Contrast specs (match your coef naming)
# rm3
spec_rm3 <- tibble(Test = c("Hydrolysable − Condensed",
                            "Total − Condensed",
                            "Hydrolysable − Total"),
                   plus_pat = c("^factor\\(Tannin_type\\)Hydro",
                                "^factor\\(Tannin_type\\)Total",
                                "^factor\\(Tannin_type\\)Hydro"),
                   minus_pat= c("^factor\\(Tannin_type\\)Condens",
                                "^factor\\(Tannin_type\\)Condens",
                                "^factor\\(Tannin_type\\)Total"))
  
# rm4
spec_rm4 <- tibble(Test = c("Condensed: Woody − Non-woody",
                            "Hydrolysable: Woody − Non-woody",
                            "Total: Woody − Non-woody"),
                   plus_pat = c("^paste\\(Tannin_type,\\s*Plant_type\\).*Condens.*\\sWoody$",
                                "^paste\\(Tannin_type,\\s*Plant_type\\).*Hydro.*\\sWoody$",
                                "^paste\\(Tannin_type,\\s*Plant_type\\).*Total.*\\sWoody$"),
                   minus_pat= c("^paste\\(Tannin_type,\\s*Plant_type\\).*Condens.*\\sNon-woody$",
                                "^paste\\(Tannin_type,\\s*Plant_type\\).*Hydro.*\\sNon-woody$",
                                "^paste\\(Tannin_type,\\s*Plant_type\\).*Total.*\\sNon-woody$"))
  
# rm5
spec_rm5 <- tibble(Test = c("Treated − Untreated",
                            "Mixed − Treated",
                            "Mixed − Untreated"),
                   plus_pat = c("^paste\\(Leaf_type\\)Treated leaves$",
                                "^paste\\(Leaf_type\\)Mixed leaves$",
                                "^paste\\(Leaf_type\\)Mixed leaves$"),
                   minus_pat= c("^paste\\(Leaf_type\\)Untreated leaves$",
                                "^paste\\(Leaf_type\\)Treated leaves$",
                                "^paste\\(Leaf_type\\)Untreated leaves$"))

# rm6
spec_rm6 <- tibble(Test = c("Mixed: Woody − Non-woody",
                            "Treated: Woody − Non-woody",
                            "Untreated: Woody − Non-woody"),
                   plus_pat = c("^paste\\(Leaf_type,\\s*Plant_type\\)Mixed leaves\\s+Woody$",
                                "^paste\\(Leaf_type,\\s*Plant_type\\)Treated leaves\\s+Woody$",
                                "^paste\\(Leaf_type,\\s*Plant_type\\)Untreated leaves\\s+Woody$"),
                   minus_pat= c("^paste\\(Leaf_type,\\s*Plant_type\\)Mixed leaves\\s+Non-woody$",
                                "^paste\\(Leaf_type,\\s*Plant_type\\)Treated leaves\\s+Non-woody$",
                                "^paste\\(Leaf_type,\\s*Plant_type\\)Untreated leaves\\s+Non-woody$"))
  
# rm7
spec_rm7 <- tibble(Test = c("Re-growth − Preexisting",
                            "Preexisting − Mixed",
                            "Re-growth − Mixed"),
                   plus_pat = c("^paste\\(Leaf_maturity\\)Re[- ]?growth leaves$",
                                "^paste\\(Leaf_maturity\\)Preexisting leaves$",
                                "^paste\\(Leaf_maturity\\)Re[- ]?growth leaves$"),
                   minus_pat= c("^paste\\(Leaf_maturity\\)Preexisting leaves$",
                                "^paste\\(Leaf_maturity\\)Mixed leaves$",
                                "^paste\\(Leaf_maturity\\)Mixed leaves$"))

# rm8
spec_rm8 <- tibble(Test = c("Mixed: Woody − Non-woody",
                            "Preexisting: Woody − Non-woody",
                            "Re-growth: Woody − Non-woody"),
                   plus_pat = c("^paste\\(Leaf_maturity,\\s*Plant_type\\)Mixed leaves\\s+Woody$",
                                "^paste\\(Leaf_maturity,\\s*Plant_type\\)Preexisting leaves\\s+Woody$",
                                "^paste\\(Leaf_maturity,\\s*Plant_type\\)Re[- ]?growth leaves\\s+Woody$"),
                   minus_pat= c("^paste\\(Leaf_maturity,\\s*Plant_type\\)Mixed leaves\\s+Non-woody$",
                                "^paste\\(Leaf_maturity,\\s*Plant_type\\)Preexisting leaves\\s+Non-woody$",
                                "^paste\\(Leaf_maturity,\\s*Plant_type\\)Re[- ]?growth leaves\\s+Non-woody$"))
  
# rm9
spec_rm9 <- tibble(Test     = "Perennial − Annual",
                   plus_pat = "^factor\\(Life_span\\)Perennial$",
                   minus_pat= "^factor\\(Life_span\\)Annual$")
  
# rm10
spec_rm10 <- tibble(Test = c( "Perennial plants: Woody vs Non-woody",
                              "Annual plants: Woody vs Non-woody"),
                    plus_pat = c("^paste\\(Life_span,\\s*Plant_type\\)Perennial\\s+Woody$",
                                 "^paste\\(Life_span,\\s*Plant_type\\)Annual\\s+Woody$"),
                    minus_pat = c("^paste\\(Life_span,\\s*Plant_type\\)Perennial\\s+Non-woody$",
                                  "^paste\\(Life_span,\\s*Plant_type\\)Annual\\s+Non-woody$"))

# rm11
spec_rm11 <- tibble(Test = c("Artificial − Chemical",
                             "Artificial − Herbivore",
                             "Artificial − Joint",
                             "Herbivore − Chemical",
                             "Herbivore − Joint",
                             "Joint − Chemical"),
                    plus_pat = c("^paste\\(Treatment_type\\)(Artifical|Artificial)$",
                                 "^paste\\(Treatment_type\\)(Artifical|Artificial)$",
                                 "^paste\\(Treatment_type\\)(Artifical|Artificial)$",
                                 "^paste\\(Treatment_type\\)Herbivore$",
                                 "^paste\\(Treatment_type\\)Herbivore$",
                                 "^paste\\(Treatment_type\\)Joint$"),
                    minus_pat= c("^paste\\(Treatment_type\\)Chemical$",
                                 "^paste\\(Treatment_type\\)Herbivore$",
                                 "^paste\\(Treatment_type\\)Joint$",
                                 "^paste\\(Treatment_type\\)Chemical$",
                                 "^paste\\(Treatment_type\\)Joint$",
                                 "^paste\\(Treatment_type\\)Chemical$"))

# rm12
spec_rm12 <- tibble(Test = c("Artificial: Woody − Non-woody",
                             "Chemical: Woody − Non-woody",
                             "Herbivore: Woody − Non-woody",
                             "Joint: Woody − Non-woody"),
                    plus_pat = c("^paste\\(Treatment_type,\\s*Plant_type\\)(Artifical|Artificial)\\s+Woody$",
                                 "^paste\\(Treatment_type,\\s*Plant_type\\)Chemical\\s+Woody$",
                                 "^paste\\(Treatment_type,\\s*Plant_type\\)Herbivore\\s+Woody$",
                                 "^paste\\(Treatment_type,\\s*Plant_type\\)Joint\\s+Woody$"),
                    minus_pat= c("^paste\\(Treatment_type,\\s*Plant_type\\)(Artifical|Artificial)\\s+Non-woody$",
                                 "^paste\\(Treatment_type,\\s*Plant_type\\)Chemical\\s+Non-woody$",
                                 "^paste\\(Treatment_type,\\s*Plant_type\\)Herbivore\\s+Non-woody$",
                                 "^paste\\(Treatment_type,\\s*Plant_type\\)Joint\\s+Non-woody$"))

# rm13
spec_rm13 <- tibble(Test = c("Defoliation − Pruning",
                             "Defoliation − Punching",
                             "Pruning − Punching"),
                    plus_pat = c("^paste\\(Treatment_typeII\\)Defoliation$",
                                 "^paste\\(Treatment_typeII\\)Defoliation$",
                                 "^paste\\(Treatment_typeII\\)Pruning$"),
                    minus_pat= c("^paste\\(Treatment_typeII\\)Pruning$",
                                 "^paste\\(Treatment_typeII\\)Punching$",
                                 "^paste\\(Treatment_typeII\\)Punching$"))
  
# rm14
spec_rm14 <- tibble(Test = c("Defoliation: Woody − Non-woody",
                             "Pruning: Woody − Non-woody",
                             "Punching: Woody − Non-woody"),
                    plus_pat = c("^paste\\(Treatment_typeII,\\s*Plant_type\\)Defoliation\\s+Woody$",
                                 "^paste\\(Treatment_typeII,\\s*Plant_type\\)Pruning\\s+Woody$",
                                 "^paste\\(Treatment_typeII,\\s*Plant_type\\)Punching\\s+Woody$"),
                    minus_pat= c("^paste\\(Treatment_typeII,\\s*Plant_type\\)Defoliation\\s+Non-woody$",
                                 "^paste\\(Treatment_typeII,\\s*Plant_type\\)Pruning\\s+Non-woody$",
                                 "^paste\\(Treatment_typeII,\\s*Plant_type\\)Punching\\s+Non-woody$"))

# rm15
spec_rm15 <- tibble(Test     = "Specialist − Generalist",
                    plus_pat = "^paste\\(Treatment_typeII\\)Specialist$",
                    minus_pat= "^paste\\(Treatment_typeII\\)Generalist$")

# rm16
spec_rm16 <- tibble(Test = c("Generalist: Woody − Non-woody",
                             "Specialist: Woody − Non-woody"),
                    plus_pat = c("^paste\\(Treatment_typeII,\\s*Plant_type\\)Generalist\\s+Woody$",
                                 "^paste\\(Treatment_typeII,\\s*Plant_type\\)Specialist\\s+Woody$"),
                    minus_pat= c("^paste\\(Treatment_typeII,\\s*Plant_type\\)Generalist\\s+Non-woody$",
                                 "^paste\\(Treatment_typeII,\\s*Plant_type\\)Specialist\\s+Non-woody$"))

# rm17
spec_rm17 <- tibble(Test     = "JA − SA",
                    plus_pat = "^paste\\(Treatment_typeII\\)JA$",
                    minus_pat= "^paste\\(Treatment_typeII\\)SA$")

# rm18
spec_rm18 <- tibble(Test = c("JA: Woody − Non-woody",
                             "SA: Woody − Non-woody"),
                    plus_pat = c("^paste\\(Treatment_typeII,\\s*Plant_type\\)JA\\s+Woody$",
                                 "^paste\\(Treatment_typeII,\\s*Plant_type\\)SA\\s+Woody$"),
                    minus_pat= c("^paste\\(Treatment_typeII,\\s*Plant_type\\)JA\\s+Non-woody$",
                                 "^paste\\(Treatment_typeII,\\s*Plant_type\\)SA\\s+Non-woody$"))
  
# rm19
spec_rm19 <- tibble(Test = c("Defol+Herb − Defol+JA",
                             "Defol+Herb − Defol+Salivation",
                             "Defol+JA − Defol+Salivation"),
                    plus_pat = c("^paste\\(Treatment_typeII\\)Defoliation \\+ Herbivore$",
                                 "^paste\\(Treatment_typeII\\)Defoliation \\+ Herbivore$",
                                 "^paste\\(Treatment_typeII\\)Defoliation \\+ JA$"),
                    minus_pat= c("^paste\\(Treatment_typeII\\)Defoliation \\+ JA$",
                                 "^paste\\(Treatment_typeII\\)Defoliation \\+ Salivation$",
                                 "^paste\\(Treatment_typeII\\)Defoliation \\+ Salivation$"))
  
# rm20
spec_rm20 <- tibble(Test = c("Defoliation + Herbivore: Woody − Non-woody",
                             "Defoliation + JA: Woody − Non-woody",
                             "Defoliation + Salivation: Woody − Non-woody"),
                    plus_pat = c("^paste\\(Treatment_typeII,\\s*Plant_type\\)Defoliation \\+ Herbivore\\s+Woody$",
                                 "^paste\\(Treatment_typeII,\\s*Plant_type\\)Defoliation \\+ JA\\s+Woody$",
                                 "^paste\\(Treatment_typeII,\\s*Plant_type\\)Defoliation \\+ Salivation\\s+Woody$"),
                    minus_pat= c("^paste\\(Treatment_typeII,\\s*Plant_type\\)Defoliation \\+ Herbivore\\s+Non-woody$",
                                 "^paste\\(Treatment_typeII,\\s*Plant_type\\)Defoliation \\+ JA\\s+Non-woody$",
                                 "^paste\\(Treatment_typeII,\\s*Plant_type\\)Defoliation \\+ Salivation\\s+Non-woody$"))

# Run all models and bind
Table_S2 <- bind_rows(.run_specs(rm3,  "Model 3: lnRR ~ Tannin type", spec_rm3),
                      .run_specs(rm4,  "Model 4: lnRR ~ Tannin type × Plant type", spec_rm4),
                      .run_specs(rm5,  "Model 5: lnRR ~ Leaf type", spec_rm5),
                      .run_specs(rm6,  "Model 6: lnRR ~ Leaf type × Plant type", spec_rm6),
                      .run_specs(rm7,  "Model 7: lnRR ~ Leaf maturity", spec_rm7),
                      .run_specs(rm8,  "Model 8: lnRR ~ Leaf maturity × Plant type", spec_rm8),
                      .run_specs(rm9,  "Model 9: lnRR ~ Life span", spec_rm9),
                      .run_specs(rm10, "Model 10: lnRR ~ Life span × Plant type", spec_rm10),
                      .run_specs(rm11, "Model 11: lnRR ~ Treatment type", spec_rm11),
                      .run_specs(rm12, "Model 12: lnRR ~ Treatment type × Plant type", spec_rm12),
                      .run_specs(rm13, "Model 13: lnRR ~ Artificial triggers", spec_rm13),
                      .run_specs(rm14, "Model 14: lnRR ~ Artificial triggers × Plant type", spec_rm14),
                      .run_specs(rm15, "Model 15: lnRR ~ Herbivore guild", spec_rm15),
                      .run_specs(rm16, "Model 16: lnRR ~ Herbivore guild × Plant type", spec_rm16),
                      .run_specs(rm17, "Model 17: lnRR ~ Chemical elicitors", spec_rm17),
                      .run_specs(rm18, "Model 18: lnRR ~ Chemical elicitors × Plant type", spec_rm18),
                      .run_specs(rm19, "Model 19: lnRR ~ Joint triggers", spec_rm19),
                      .run_specs(rm20, "Model 20: lnRR ~ Joint triggers × Plant type", spec_rm20)) %>%
  mutate(bold_95CI_excludes0 = ifelse(is.na(CI_low), NA, (CI_low > 0) | (CI_high < 0)),
         Estimate = round(Estimate, 3),
         SE       = round(SE, 3),
         z        = round(z, 3),
         p        = ifelse(is.na(p), NA_real_, p),
         p        = round(p, 3),
         CI_low   = round(CI_low, 3),
         CI_high  = round(CI_high, 3))
print(Table_S2, n = Inf)
write.csv(Table_S2, "Table_pairwise_contrasts_rm3_rm20.csv", row.names = FALSE)

#============================
# Export Table (rm3 ... rm20)
#============================

library(metafor)
library(dplyr)

models <- list(rm3=rm3, rm4=rm4, 
               rm5=rm5, rm6=rm6, 
               rm7=rm7, rm8=rm8,
               rm9=rm9, rm10=rm10,
               rm11=rm11, rm12=rm12,
               rm13=rm13, rm14=rm14,
               rm15=rm15, rm16=rm16,
               rm17=rm17, rm18=rm18,
               rm19=rm19, rm20=rm20)

model_labels <- c(rm3 ="Model 3: Effect size (lnRR) ~ Tannin type",
                  rm4 ="Model 4: Effect size (lnRR) ~ Tannin type × Plant type",
                  rm5 ="Model 5: Effect size (lnRR) ~ Leaf type",
                  rm6 ="Model 6: Effect size (lnRR) ~ Leaf type × Plant type",
                  rm7 ="Model 7: Effect size (lnRR) ~ Leaf maturity",
                  rm8 ="Model 8: Effect size (lnRR) ~ Leaf maturity × Plant type",
                  rm9 ="Model 9: Effect size (lnRR) ~ Life span",
                  rm10="Model 10: Effect size (lnRR) ~ Life span × Plant type",
                  rm11="Model 11: Effect size (lnRR) ~ Treatment type",
                  rm12="Model 12: Effect size (lnRR) ~ Treatment type × Plant type",
                  rm13="Model 13: Effect size (lnRR) ~ Artificial triggers",
                  rm14="Model 14: Effect size (lnRR) ~ Artificial triggers × Plant type",
                  rm15="Model 15: Effect size (lnRR) ~ Herbivore guild",
                  rm16="Model 16: Effect size (lnRR) ~ Herbivore guild × Plant type",
                  rm17="Model 17: Effect size (lnRR) ~ Chemical elicitors",
                  rm18="Model 18: Effect size (lnRR) ~ Chemical elicitors × Plant type",
                  rm19="Model 19: Effect size (lnRR) ~ Joint triggers",
                  rm20="Model 20: Effect size (lnRR) ~ Joint triggers × Plant type")

fmt_p <- function(p) {
  ifelse(is.na(p), NA_character_,
         ifelse(p < 0.001, "<0.001",
                formatC(p, digits = 3, format = "f")))
}

clean_term <- function(x) {
  x2 <- gsub("^paste\\([^\\)]*\\)", "", x)
  x2 <- gsub("^factor\\([^\\)]*\\)", "", x2)
  x2 <- gsub("^0 \\+ ", "", x2)
  x2 <- gsub("^\\s+|\\s+$", "", x2)
  x2 <- ifelse(x2 %in% c("intrcpt", "(Intercept)", "Intercept"), "Full model", x2)
  x2
}

get_df_QE <- function(fit) {
  k <- tryCatch(as.numeric(fit$k), error=function(e) NA_real_)
  p <- tryCatch(as.numeric(fit$p), error=function(e) NA_real_)
  if (!is.na(k) && !is.na(p)) return(k - p)
  NA_real_
}

variance_line <- function(fit) {
  tau <- tryCatch(as.numeric(fit$sigma2), error=function(e) numeric(0))
  QE  <- tryCatch(as.numeric(fit$QE),  error=function(e) NA_real_)
  QEp <- tryCatch(as.numeric(fit$QEp), error=function(e) NA_real_)
  dfQ <- get_df_QE(fit)
  
  tau_txt <- if (length(tau) == 0) {
    "τ² = NA"
  } else {
    paste(sprintf("τ%d² = %.3f", seq_along(tau), tau), collapse = ", ")
  }
  
  qe_txt <- paste0("QE = ", ifelse(is.na(QE), "NA", formatC(QE, digits=3, format="f")),
                   ", df = ", ifelse(is.na(dfQ), "NA", formatC(dfQ, digits=0, format="f")),
                   ", p ", ifelse(is.na(QEp), "= NA", paste0(ifelse(QEp < 0.001, "< ", "= "), fmt_p(QEp))))
  
  paste0("(", tau_txt, "; ", qe_txt, ")")
}

coef_table <- function(fit) {
  term <- tryCatch(rownames(fit$b), error=function(e) NULL)
  est  <- tryCatch(as.numeric(fit$b), error=function(e) NA_real_)
  se   <- tryCatch(as.numeric(fit$se), error=function(e) rep(NA_real_, length(est)))
  z    <- tryCatch(as.numeric(fit$zval), error=function(e) rep(NA_real_, length(est)))
  p    <- tryCatch(as.numeric(fit$pval), error=function(e) rep(NA_real_, length(est)))
  lci  <- tryCatch(as.numeric(fit$ci.lb), error=function(e) rep(NA_real_, length(est)))
  uci  <- tryCatch(as.numeric(fit$ci.ub), error=function(e) rep(NA_real_, length(est)))
  
  if (is.null(term)) term <- paste0("b", seq_along(est))
  
  tibble(Moderator_variables = clean_term(term),
         Estimate = est,
         SE = se,
         z = z,
         p = p,
         LCI = lci,
         UCI = uci) %>%
    mutate(sig_CI = ifelse(is.na(LCI) | is.na(UCI), NA, (LCI > 0) | (UCI < 0)),
           Estimate = round(Estimate, 3),
           SE       = round(SE, 3),
           z        = round(z, 3),
           p_txt    = fmt_p(p),
           LCI      = round(LCI, 3),
           UCI      = round(UCI, 3)) %>%
    select(Moderator_variables, Estimate, SE, z, p_txt, LCI, UCI, sig_CI)
}

build_TableS1 <- function(models, labels) {
  out <- list()
  
  for (nm in names(models)) {
    fit <- models[[nm]]
    model_title <- ifelse(!is.na(labels[nm]), labels[nm], nm)
    
    r1 <- tibble(Moderator_variables = model_title,
                 Estimate = NA_real_, SE = NA_real_, z = NA_real_, p_txt = NA_character_,
                 LCI = NA_real_, UCI = NA_real_, sig_CI = NA)
    r2 <- tibble(Moderator_variables = variance_line(fit),
                 Estimate = NA_real_, SE = NA_real_, z = NA_real_, p_txt = NA_character_,
                 LCI = NA_real_, UCI = NA_real_, sig_CI = NA)
    r3 <- coef_table(fit)
    out[[nm]] <- bind_rows(r1, r2, r3,
                           tibble(Moderator_variables = "", Estimate=NA_real_, SE=NA_real_, z=NA_real_,
                                  p_txt=NA_character_, LCI=NA_real_, UCI=NA_real_, sig_CI=NA))
  }
  bind_rows(out) %>% rename(p = p_txt)
}

Table_S1 <- build_TableS1(models, model_labels)
out_dir <- "Export_rm3_rm20_TableS1"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
write.csv(Table_S1, file.path(out_dir, "Table_S1_rm3_rm20.csv"), row.names = FALSE)

if (requireNamespace("knitr", quietly = TRUE) && requireNamespace("kableExtra", quietly = TRUE)) {
  library(knitr)
  library(kableExtra)
  Table_S1_html <- Table_S1 %>%
    mutate(
      Moderator_variables = ifelse(sig_CI %in% TRUE,
                                   kableExtra::cell_spec(Moderator_variables, bold = TRUE),
                                   Moderator_variables),
      Estimate = ifelse(sig_CI %in% TRUE,
                        kableExtra::cell_spec(as.character(Estimate), bold = TRUE),
                        as.character(Estimate)),
      SE = ifelse(sig_CI %in% TRUE,
                  kableExtra::cell_spec(as.character(SE), bold = TRUE),
                  as.character(SE)),
      z = ifelse(sig_CI %in% TRUE,
                 kableExtra::cell_spec(as.character(z), bold = TRUE),
                 as.character(z)),
      p = ifelse(sig_CI %in% TRUE,
                 kableExtra::cell_spec(as.character(p), bold = TRUE),
                 as.character(p)),
      LCI = ifelse(sig_CI %in% TRUE,
                   kableExtra::cell_spec(as.character(LCI), bold = TRUE),
                   as.character(LCI)),
      UCI = ifelse(sig_CI %in% TRUE,
                   kableExtra::cell_spec(as.character(UCI), bold = TRUE),
                   as.character(UCI))) %>%
    select(Moderator_variables, Estimate, SE, z, p, LCI, UCI)
  html_file <- file.path(out_dir, "Table_S1_rm3_rm20.html")
  kbl(Table_S1_html, escape = FALSE,
      col.names = c("Moderator variables","Estimate","SE","z","p","LCI","UCI")) %>%
    kableExtra::kable_styling(full_width = FALSE) %>%
    kableExtra::save_kable(file = html_file)
}
cat("Done. Outputs saved in: ", normalizePath(out_dir), "\n")

# =========================================
# Pagel's lambda check for tannin induction
# =========================================

options(stringsAsFactors = FALSE)
library(readxl)
library(dplyr)
library(stringr)
library(ape)
library(phytools)
library(V.PhyloMaker2)

N_ITER <- 1000
SEED   <- 20251221
WEIGHTED_PICK <- FALSE

cat(">>> Step 1. Choose your file...\n")
data_path <- file.choose()
out_dir   <- dirname(data_path)
cat(">>> Reading:", data_path, "\n")
dat0 <- readxl::read_excel(data_path)
dat0 <- as.data.frame(dat0)

# ---- Check required columns ----
needed <- c("lnRR", "rvar", "Family", "Plant_species", "Plant_type")
missing <- setdiff(needed, names(dat0))
if (length(missing) > 0) {
  stop("Missing required columns in Excel: ", paste(missing, collapse = ", "), call. = FALSE)
}

# Clean & standardize
# - Species uniqueness: keep only the first two words (Genus + species epithet)
dat <- dat0 %>%
  mutate(
    Plant_species = as.character(Plant_species),
    Plant_species = gsub("\ufeff", "", Plant_species),   # remove BOM
    Plant_species = str_squish(Plant_species),
    
    Family        = str_squish(as.character(Family)),
    Plant_type    = str_squish(as.character(Plant_type)),
    
    lnRR          = suppressWarnings(as.numeric(lnRR)),
    rvar          = suppressWarnings(as.numeric(rvar)),
    
    genus      = str_to_title(word(Plant_species, 1)),
    species_ep = str_to_lower(word(Plant_species, 2)),
    
    Plant_species_binom = paste(genus, species_ep),           # "Genus species"
    Plant_species_us    = paste(genus, species_ep, sep = "_") # "Genus_species"
  ) %>%
  filter(!is.na(genus), !is.na(species_ep)) %>%
  filter(!(tolower(species_ep) %in% c("sp", "sp.", "spp", "spp."))) %>%
  filter(is.finite(lnRR))

if (nrow(dat) == 0) {
  stop("Data are empty after cleaning: please check that Plant_species uses standard binomial names (Genus species).", call. = FALSE)
}

# ---- Standardize Plant_type labels (extend synonyms as needed for your dataset) ----
dat <- dat %>%
  mutate(
    Plant_type = case_when(
      str_to_lower(Plant_type) %in% c("woody") ~ "Woody",
      str_to_lower(Plant_type) %in% c("non-woody", "nonwoody", "non woody", "herbaceous") ~ "Non-woody",
      TRUE ~ Plant_type
    )
  )

# Species counts (Total / Woody / Non-woody)
cat("\n--- BASIC COUNTS ---\n")
cat("Rows after cleaning:", nrow(dat), "\n")
cat("Unique species (Plant_species_us):", n_distinct(dat$Plant_species_us), "\n")

species_counts <- dat %>%
  distinct(Plant_species_us, Plant_type) %>%
  count(Plant_type, name = "n_species") %>%
  arrange(desc(n_species))

cat("\nUnique species by Plant_type:\n")
print(species_counts)

# check whether the same species is labeled with multiple Plant_type values
type_conflict <- dat %>%
  distinct(Plant_species_us, Plant_type) %>%
  count(Plant_species_us, name = "n_types") %>%
  filter(n_types > 1) %>%
  left_join(
    dat %>%
      distinct(Plant_species_us, Plant_type) %>%
      group_by(Plant_species_us) %>%
      summarise(types = paste(sort(unique(Plant_type)), collapse = " | "), .groups = "drop"),
    by = "Plant_species_us"
  ) %>%
  arrange(desc(n_types), Plant_species_us)

if (nrow(type_conflict) > 0) {
  cat("\n!!! WARNING: Some species have MULTIPLE Plant_type labels:\n")
  print(type_conflict, row.names = FALSE)
  write.csv(type_conflict, file.path(out_dir, "species_type_conflicts.csv"), row.names = FALSE)
  cat(">>> Saved: species_type_conflicts.csv\n")
} else {
  cat("\nOK: No Plant_type conflicts across the same species.\n")
}

# check whether multiple raw names collapse into the same binomial
collapsed_map <- dat %>%
  distinct(Plant_species, Plant_species_us) %>%
  group_by(Plant_species_us) %>%
  summarise(
    n_raw = n(),
    raw_names = paste(Plant_species, collapse = " | "),
    .groups = "drop"
  ) %>%
  filter(n_raw > 1) %>%
  arrange(desc(n_raw))

if (nrow(collapsed_map) > 0) {
  cat("\n!!! Potential collapse detected (multiple raw names -> same Plant_species_us):\n")
  print(collapsed_map, row.names = FALSE)
  write.csv(collapsed_map, file.path(out_dir, "species_name_collapse_check.csv"), row.names = FALSE)
  cat(">>> Saved: species_name_collapse_check.csv\n")
} else {
  cat("\nOK: No collapse detected from Plant_species -> Plant_species_us.\n")
}

# Export species list (for use in the main text)
species_list_by_type <- dat %>%
  distinct(Plant_species_us, Plant_species_binom, Plant_type, Family) %>%
  arrange(Plant_type, Plant_species_us)

write.csv(species_list_by_type, file.path(out_dir, "species_list_by_type.csv"), row.names = FALSE)
cat(">>> Saved: species_list_by_type.csv\n")

# missing Family values (may cause failures in tree construction/matching)
fam_missing <- dat %>%
  filter(is.na(Family) | Family == "") %>%
  distinct(Plant_species, Plant_species_us, Family, Plant_type) %>%
  arrange(Plant_species_us)

if (nrow(fam_missing) > 0) {
  cat("\n!!! WARNING: Missing Family for these records (may prevent placement into the tree):\n")
  print(fam_missing, row.names = FALSE)
  write.csv(fam_missing, file.path(out_dir, "missing_family_records.csv"), row.names = FALSE)
  cat(">>> Saved: missing_family_records.csv\n")
}

# Build phylogenetic tree with V.PhyloMaker2
sp_list <- dat %>%
  distinct(genus, species_ep, Family) %>%
  filter(!is.na(Family), Family != "") %>%
  mutate(species = paste(genus, species_ep)) %>%   # Key fix
  transmute(
    species = species,  # "Genus species"
    genus   = genus,
    family  = Family
  )

# Safety check: detect duplicate binomials (often indicates the same species assigned to different Family values or inconsistent spelling)
dup_binom <- sp_list %>% count(species, name = "n") %>% filter(n > 1)
if (nrow(dup_binom) > 0) {
  cat("\n!!! ERROR: Duplicate binomials in sp_list (same 'Genus species' appears >1):\n")
  print(dup_binom)
  print(sp_list %>% semi_join(dup_binom, by = "species") %>% arrange(species))
  stop("Please fix Family values or species spelling (ensure each binomial appears only once).", call. = FALSE)
}

cat("\n>>> Step 2. Building phylogenetic tree with V.PhyloMaker2...\n")
phylo_res <- V.PhyloMaker2::phylo.maker(sp.list = sp_list)

tree    <- phylo_res$scenario.3
matched <- phylo_res$species.list

# ---- Make tip labels consistent with Plant_species_us ("Genus_species") ----
tip_src <- matched$species
if (!any(str_detect(tip_src, "\\s"))) {
  if (all(c("genus", "species") %in% names(matched))) {
    tip_src <- paste(matched$genus, matched$species)
  }
}
tip_src <- str_squish(tip_src)
tree$tip.label <- str_replace_all(tip_src, " ", "_")

# Standardize case: Genus capitalized, species lowercase
tip_parts <- str_split_fixed(tree$tip.label, "_", 2)
tree$tip.label <- paste0(str_to_title(tip_parts[,1]), "_", str_to_lower(tip_parts[,2]))

# Branch lengths safety
if (is.null(tree$edge.length)) {
  tree$edge.length <- rep(1, nrow(tree$edge))
  cat(">>> Tree has no branch lengths; set all edge lengths to 1.\n")
}

# duplicated tip labels (will reduce the number of tips)
dup_tip <- unique(tree$tip.label[duplicated(tree$tip.label)])
if (length(dup_tip) > 0) {
  cat("\n!!! ERROR: Duplicated tip labels in tree:\n")
  print(dup_tip)
  stop("Duplicated tip labels in the tree indicate a naming conflict (check genus/species spelling or duplicate binomials).", call. = FALSE)
}

cat(">>> Tips in tree:", length(tree$tip.label), "\n")
cat(">>> Unique tips in tree:", length(unique(tree$tip.label)), "\n")

# V.PhyloMaker2 match status
if ("status" %in% names(matched)) {
  cat("\n--- V.PhyloMaker2 species match status ---\n")
  print(table(matched$status))
  bad <- matched[matched$status != 1, , drop = FALSE]
  if (nrow(bad) > 0) {
    cat("\n!!! Species not fully matched (status != 1):\n")
    print(bad)
    write.csv(bad, file.path(out_dir, "vphylomaker_unmatched_species.csv"), row.names = FALSE)
    cat(">>> Saved: vphylomaker_unmatched_species.csv\n")
  } else {
    cat("OK: All species status == 1\n")
  }
} else {
  cat("\n(Note) matched$status not found; skip status QA.\n")
}

# DEBUG: which species are lost (DATA vs TREE)
sp_clean <- sort(unique(dat$Plant_species_us))
sp_tree  <- sort(tree$tip.label)

cat("\n--- Species counts by step ---\n")
cat("Cleaned species (data):", length(sp_clean), "\n")
cat("Tree tips:", length(sp_tree), "\n")

miss_in_tree <- setdiff(sp_clean, sp_tree)
if (length(miss_in_tree) > 0) {
  cat("\n!!! Species in DATA but NOT in TREE (this is the most common direct reason for a species-count drop, e.g., 78 -> 77):\n")
  print(miss_in_tree)
  
  miss_tbl <- dat %>%
    filter(Plant_species_us %in% miss_in_tree) %>%
    distinct(Plant_species, Plant_species_us, Family, Plant_type) %>%
    arrange(Plant_species_us)
  
  cat("\nTheir records (check Family / spelling):\n")
  print(miss_tbl, row.names = FALSE)
  
  write.csv(miss_tbl, file.path(out_dir, "species_missing_in_tree.csv"), row.names = FALSE)
  cat(">>> Saved: species_missing_in_tree.csv\n")
} else {
  cat("\nOK: All cleaned species are present in the tree.\n")
}

# Keep only species present in tree
common_sp <- intersect(tree$tip.label, dat$Plant_species_us)
dat_use <- dat %>% filter(Plant_species_us %in% common_sp)

cat("\n>>> Common species (tree ∩ data):", n_distinct(dat_use$Plant_species_us), "\n")

# Species-level aggregation (mean & inverse-variance weighted mean)
species_eff <- dat_use %>%
  group_by(Plant_species_us, Plant_type) %>%
  summarise(k = n(),
            lnRR_mean = mean(lnRR, na.rm = TRUE),
            lnRR_wmean = {
              w <- 1 / rvar
              w[!is.finite(w) | w <= 0] <- NA_real_
              if (all(is.na(w))) mean(lnRR, na.rm = TRUE) else weighted.mean(lnRR, w = w, na.rm = TRUE)
            },
            .groups = "drop")

write.csv(species_eff, file.path(out_dir, "species_level_effect_sizes.csv"), row.names = FALSE)
cat(">>> Saved: species_level_effect_sizes.csv\n")

calc_lambda <- function(tr, x_named) {
  common <- intersect(tr$tip.label, names(x_named))
  if (length(common) < 3) {
    return(data.frame(lambda = NA_real_, logL = NA_real_, P = NA_real_, n_tips = length(common)))
  }
  
  tr2 <- drop.tip(tr, setdiff(tr$tip.label, common))
  x2  <- as.numeric(x_named[tr2$tip.label])
  names(x2) <- tr2$tip.label
  
  out <- tryCatch(phytools::phylosig(tr2, x2, method = "lambda", test = TRUE),
                  error = function(e) NULL)

  if (is.null(out)) {
    return(data.frame(lambda = NA_real_, logL = NA_real_, P = NA_real_, n_tips = length(tr2$tip.label)))
  }
  data.frame(lambda = out$lambda, logL = out$logL, P = out$P, n_tips = length(tr2$tip.label))
}

# ---- lambda for species-level mean vs weighted mean ----
x_mean  <- setNames(species_eff$lnRR_mean,  species_eff$Plant_species_us)
x_wmean <- setNames(species_eff$lnRR_wmean, species_eff$Plant_species_us)

lambda_all_mean  <- calc_lambda(tree, x_mean)
lambda_all_wmean <- calc_lambda(tree, x_wmean)

cat("\n=== Pagel's lambda (All species) ===\n")
cat("Species-level MEAN:\n");  print(lambda_all_mean)
cat("Species-level INV-VAR weighted mean:\n"); print(lambda_all_wmean)

# ---- Group-specific (woody / non-woody) using weighted mean ----
calc_group_lambda <- function(group_label) {
  spg <- species_eff %>% filter(Plant_type == group_label)
  xg  <- setNames(spg$lnRR_wmean, spg$Plant_species_us)
  calc_lambda(tree, xg)
}

lambda_woody    <- calc_group_lambda("Woody")
lambda_nonwoody <- calc_group_lambda("Non-woody")

cat("\n=== Pagel's lambda (Group-specific; weighted mean) ===\n")
cat("Woody:\n");     print(lambda_woody)
cat("Non-woody:\n"); print(lambda_nonwoody)

# Stratified resampling (1 effect size per species)
pick_one_row <- function(df, weighted_pick = FALSE) {
  if (!weighted_pick) {
    df[sample.int(nrow(df), 1), , drop = FALSE]
  } else {
    w <- 1 / df$rvar
    w[!is.finite(w) | w <= 0] <- 0
    if (sum(w) <= 0) df[sample.int(nrow(df), 1), , drop = FALSE]
    else df[sample.int(nrow(df), 1, prob = w), , drop = FALSE]
  }
}

resample_lambda <- function(dat_long, tr, group = NULL, n_iter = 1000, seed = 1, weighted_pick = FALSE) {
  set.seed(seed)
  dd <- dat_long
  if (!is.null(group)) dd <- dd %>% filter(Plant_type == group)
  
  dd <- dd %>%
    filter(Plant_species_us %in% tr$tip.label) %>%
    filter(is.finite(lnRR))
  
  nsp <- n_distinct(dd$Plant_species_us)
  if (nsp < 3) {
    warning("Too few species to estimate lambda reliably. group = ", group)
    return(data.frame(iter = integer(), lambda = numeric(), P = numeric(), n_tips = integer(), group = character()))
  }
  
  out <- vector("list", n_iter)
  for (i in seq_len(n_iter)) {
    samp <- dd %>%
      group_by(Plant_species_us) %>%
      group_modify(~ pick_one_row(.x, weighted_pick = weighted_pick)) %>%
      ungroup()
    
    x <- setNames(samp$lnRR, samp$Plant_species_us)
    lam <- calc_lambda(tr, x)
    
    out[[i]] <- data.frame(
      iter   = i,
      lambda = lam$lambda,
      P      = lam$P,
      n_tips = lam$n_tips,
      group  = ifelse(is.null(group), "All", group)
    )
    
    if (i %% 100 == 0) cat("... resampling iter", i, "group =", ifelse(is.null(group), "All", group), "\n")
  }
  bind_rows(out)
}

cat("\n>>> Running stratified resampling...\n")
res_all      <- resample_lambda(dat_use, tree, group = NULL,        n_iter = N_ITER, seed = SEED, weighted_pick = WEIGHTED_PICK)
res_woody    <- resample_lambda(dat_use, tree, group = "Woody",     n_iter = N_ITER, seed = SEED, weighted_pick = WEIGHTED_PICK)
res_nonwoody <- resample_lambda(dat_use, tree, group = "Non-woody", n_iter = N_ITER, seed = SEED, weighted_pick = WEIGHTED_PICK)

# ---- Summaries ----
summarise_res <- function(res) {
  if (nrow(res) == 0) return(res)
  res %>%
    summarise(
      n = n(),
      n_valid = sum(!is.na(lambda)),
      lambda_median = median(lambda, na.rm = TRUE),
      lambda_q025   = quantile(lambda, 0.025, na.rm = TRUE),
      lambda_q975   = quantile(lambda, 0.975, na.rm = TRUE),
      prop_P_lt_0.05 = mean(P < 0.05, na.rm = TRUE),
      n_tips_median = median(n_tips, na.rm = TRUE)
    )
}

sum_all      <- summarise_res(res_all)
sum_woody    <- summarise_res(res_woody)
sum_nonwoody <- summarise_res(res_nonwoody)

cat("\n=== Resampling summary (All) ===\n");      print(sum_all)
cat("\n=== Resampling summary (Woody) ===\n");    print(sum_woody)
cat("\n=== Resampling summary (Non-woody) ===\n");print(sum_nonwoody)

lambda_summary <- bind_rows(
  sum_all      %>% mutate(group = "All"),
  sum_woody    %>% mutate(group = "Woody"),
  sum_nonwoody %>% mutate(group = "Non-woody")
) %>% select(group, everything())

write.csv(lambda_summary, file.path(out_dir, "lambda_resampling_summary.csv"), row.names = FALSE)
cat("\n>>> Saved: lambda_resampling_summary.csv\n")

# ---- Save full resampling outputs ----
write.csv(res_all,      file.path(out_dir, "lambda_resampling_all.csv"),      row.names = FALSE)
write.csv(res_woody,    file.path(out_dir, "lambda_resampling_woody.csv"),    row.names = FALSE)
write.csv(res_nonwoody, file.path(out_dir, "lambda_resampling_nonwoody.csv"), row.names = FALSE)

cat("\n>>> Saved:\n")
cat(" - lambda_resampling_all.csv\n")
cat(" - lambda_resampling_woody.csv\n")
cat(" - lambda_resampling_nonwoody.csv\n")

# ---- Quick histograms ----
pdf(file.path(out_dir, "lambda_resampling_histograms.pdf"), width = 7, height = 5)
hist(res_all$lambda, breaks = 30, main = "Resampled Pagel's lambda (All)", xlab = "lambda")
hist(res_woody$lambda, breaks = 30, main = "Resampled Pagel's lambda (Woody)", xlab = "lambda")
hist(res_nonwoody$lambda, breaks = 30, main = "Resampled Pagel's lambda (Non-woody)", xlab = "lambda")
dev.off()
cat(">>> Saved: lambda_resampling_histograms.pdf\n")

cat("\n✅ Done. All outputs saved to:\n", out_dir, "\n", sep = "")

n_sp_all <- dplyr::n_distinct(species_eff$Plant_species_us)

n_sp_woody <- species_eff %>%
  dplyr::filter(Plant_type == "Woody") %>%
  dplyr::pull(Plant_species_us) %>%
  dplyr::n_distinct()
n_sp_woody

n_sp_nonwoody <- species_eff %>%
  dplyr::filter(Plant_type == "Non-woody") %>%
  dplyr::pull(Plant_species_us) %>%
  dplyr::n_distinct()
n_sp_nonwoody

# Combine calc_lambda outputs (data.frame: lambda, logL, P, n_tips) into a single table
lambda_specieslevel_summary <- dplyr::bind_rows(
  lambda_all_mean  %>% dplyr::mutate(Group = "All",      Aggregation = "Species-level mean",                       n_species = n_sp_all),
  lambda_all_wmean %>% dplyr::mutate(Group = "All",      Aggregation = "Species-level inverse-variance weighted",  n_species = n_sp_all),
  lambda_woody     %>% dplyr::mutate(Group = "Woody",    Aggregation = "Species-level inverse-variance weighted",  n_species = n_sp_woody),
  lambda_nonwoody  %>% dplyr::mutate(Group = "Non-woody",Aggregation = "Species-level inverse-variance weighted",  n_species = n_sp_nonwoody)
) %>%
  dplyr::select(Group, Aggregation, n_species, n_tips, lambda, logL, P)

# Export
write.csv(lambda_specieslevel_summary,
          file.path(out_dir, "lambda_specieslevel_summary.csv"),
          row.names = FALSE)

cat(">>> Saved: lambda_specieslevel_summary.csv\n")
print(lambda_specieslevel_summary)

#=======================================================
# Subgroup dataset analysis (rm21 - rm38) & Export Table
#=======================================================

setwd("/Users/yaolin/Desktop/My papers/Manuscripts/2025 - Induced tannin/Submission - Plant Diversity - 0209/Data & Code")
data <- read_excel("Raw data - Induced tannin - 0211.xlsx")
write.csv(data, file = "data.csv", row.names = FALSE)
data <- read.csv("data.csv")

data1 <- data %>% filter(Plant_type == "Woody")
data2 <- data %>% filter(Plant_type == "Non-woody")

table(data$Plant_type, useNA = "ifany")

rm21 <- rma.mv(lnRR, rvar,
               mods = ~ factor(Tannin_type) - 1,
               random = list( ~ 1 | Publication/Effect_ID, ~ 1 | Plant_species),
               method = "REML", 
               data = data1)
rm21

rm22 <- rma.mv(lnRR, rvar,
               mods = ~ factor(Tannin_type) - 1,
               random = list( ~ 1 | Publication/Effect_ID, ~ 1 | Plant_species),
               method = "REML", 
               data = data2)
rm22

rm23 <- rma.mv(lnRR, rvar,
               mods = ~ factor(Leaf_type) - 1,
               random = list( ~ 1 | Publication/Effect_ID, ~ 1 | Plant_species),
               method = "REML", 
               data = data1)
rm23

rm24 <- rma.mv(lnRR, rvar,
               mods = ~ factor(Leaf_type) - 1,
               random = list( ~ 1 | Publication/Effect_ID, ~ 1 | Plant_species),
               method = "REML", 
               data = data2)
rm24

rm25 <- rma.mv(lnRR, rvar,
               mods = ~ paste(Leaf_maturity) - 1,
               random = list( ~ 1 | Publication/Effect_ID, ~ 1 | Plant_species),
               method = "REML", 
               data = data1)
rm25

rm26 <- rma.mv(lnRR, rvar,
               mods = ~ paste(Leaf_maturity) - 1,
               random = list( ~ 1 | Publication/Effect_ID, ~ 1 | Plant_species),
               method = "REML", 
               data = data2)
rm26

rm27 <- rma.mv(lnRR, rvar,
               mods = ~ 0 + factor(Life_span),
               random = list(~ 1 | Publication/Effect_ID, ~ 1 | Plant_species),
               method = "REML",
               data = data1)
rm27

rm28 <- rma.mv(lnRR, rvar,
               mods = ~ 0 + factor(Life_span),
               random = list(~ 1 | Publication/Effect_ID, ~ 1 | Plant_species),
               method = "REML",
               data = data2)
rm28

rm29 <- rma.mv(lnRR, rvar,
               mods = ~ 0 + paste(Treatment_type),
               random = list( ~ 1 | Publication/Effect_ID, ~ 1 | Plant_species),
               method = "REML", 
               data = data1)
rm29

rm30 <- rma.mv(lnRR, rvar,
               mods = ~ 0 + paste(Treatment_type),
               random = list( ~ 1 | Publication/Effect_ID, ~ 1 | Plant_species),
               method = "REML", 
               data = data2)
rm30

data1A <- data1 %>%
  filter(Treatment_type == "Artificial") %>%
  droplevels()                              

data2A <- data2 %>%
  filter(Treatment_type == "Artificial") %>%
  droplevels() 

rm31 <- rma.mv(lnRR, rvar,
               mods = ~ paste(Treatment_typeII) - 1,
               random = list( ~ 1 | Publication/Effect_ID, ~ 1 | Plant_species),
               method = "REML", 
               data = data1A)
rm31

rm32 <- rma.mv(lnRR, rvar,
               mods = ~ paste(Treatment_typeII) - 1,
               random = list( ~ 1 | Publication/Effect_ID, ~ 1 | Plant_species),
               method = "REML", 
               data = data2A)
rm32

data1H <- data1 %>%
  filter(Treatment_type == "Herbivore") %>%
  droplevels()                              

data2H <- data2 %>%
  filter(Treatment_type == "Herbivore") %>%
  droplevels()   

rm33 <- rma.mv(lnRR, rvar,
               mods = ~ paste(Treatment_typeII) - 1,
               random = list( ~ 1 | Publication/Effect_ID, ~ 1 | Plant_species),
               method = "REML", 
               data = data1H)
rm33

rm34 <- rma.mv(lnRR, rvar,
               mods = ~ paste(Treatment_typeII) - 1,
               random = list( ~ 1 | Publication/Effect_ID, ~ 1 | Plant_species),
               method = "REML", 
               data = data2H)
rm34

data1C <- data1 %>%
  filter(Treatment_type == "Chemical") %>%
  droplevels()                              

data2C <- data2 %>%
  filter(Treatment_type == "Chemical") %>%
  droplevels()                              

rm35 <- rma.mv(lnRR, rvar,
               mods = ~ paste(Treatment_typeII) - 1,
               random = list( ~ 1 | Publication/Effect_ID, ~ 1 | Plant_species),
               method = "REML", 
               data = data1C)
rm35

rm36 <- rma.mv(lnRR, rvar,
               mods = ~ paste(Treatment_typeII) - 1,
               random = list( ~ 1 | Publication/Effect_ID, ~ 1 | Plant_species),
               method = "REML", 
               data = data2C)
rm36

data1J <- data1 %>%
  filter(Treatment_type == "Joint") %>%
  droplevels()                              

data2J <- data2 %>%
  filter(Treatment_type == "Joint") %>%
  droplevels()    

rm37 <- rma.mv(lnRR, rvar,
               mods = ~ paste(Treatment_typeII) - 1,
               random = list( ~ 1 | Publication/Effect_ID, ~ 1 | Plant_species),
               method = "REML", 
               data = data1J)
rm37

rm38 <- rma.mv(lnRR, rvar,
               mods = ~ paste(Treatment_typeII) - 1,
               random = list( ~ 1 | Publication/Effect_ID, ~ 1 | Plant_species),
               method = "REML", 
               data = data2J)
rm38

# Export Table

library(dplyr)
library(tibble)

models <- list(rm21=rm21, rm22=rm22, rm23=rm23, rm24=rm24, rm25=rm25, rm26=rm26,
               rm28=rm28, rm29=rm29, rm30=rm30, rm31=rm31, rm32=rm32,
               rm33=rm33, rm34=rm34, rm35=rm35, rm36=rm36, rm37=rm37, rm38=rm38)

labels <- c(rm21="Model 21: Tannin type (Woody)",
            rm22="Model 22: Tannin type (Non-woody)",
            rm23="Model 23: Leaf type (Woody)",
            rm24="Model 24: Leaf type (Non-woody)",
            rm25="Model 25: Leaf maturity (Woody)",
            rm26="Model 26: Leaf maturity (Non-woody)",
            rm28="Model 28: Life span (Non-woody)",
            rm29="Model 29: Treatment type (Woody)",
            rm30="Model 30: Treatment type (Non-woody)",
            rm31="Model 31: Artificial triggers (Woody)",
            rm32="Model 32: Artificial triggers (Non-woody)",
            rm33="Model 33: Herbivore guild (Woody)",
            rm34="Model 34: Herbivore guild (Non-woody)",
            rm35="Model 35: Chemical elicitors (Woody)",
            rm36="Model 36: Chemical elicitors (Non-woody)",
            rm37="Model 37: Joint triggers (Woody)",
            rm38="Model 38: Joint triggers (Non-woody)")

fmt_p <- function(p) {
  ifelse(is.na(p), NA_character_,
         ifelse(p < 0.001, "<0.001",
                formatC(p, digits = 3, format = "f")))
}

clean_term <- function(x) {
  x2 <- gsub("^paste\\([^\\)]*\\)", "", x)
  x2 <- gsub("^factor\\([^\\)]*\\)", "", x)
  x2 <- gsub("^0 \\+ ", "", x2)
  x2 <- gsub("^\\s+|\\s+$", "", x2)
  x2
}

variance_line <- function(fit) {
  if (is.null(fit)) return("(model not fitted)")
  
  s2 <- tryCatch(as.numeric(fit$sigma2), error=function(e) numeric(0))
  s2 <- c(s2, rep(NA_real_, max(0, 3 - length(s2))))[1:3]
  
  QE  <- tryCatch(as.numeric(fit$QE),  error=function(e) NA_real_)
  QEp <- tryCatch(as.numeric(fit$QEp), error=function(e) NA_real_)
  dfQ <- tryCatch(as.numeric(fit$k - fit$p), error=function(e) NA_real_)
  
  paste0(
    "(τ1² = ", ifelse(is.na(s2[1]), "NA", sprintf("%.3f", s2[1])),
    ", τ2² = ", ifelse(is.na(s2[2]), "NA", sprintf("%.3f", s2[2])),
    ", τ3² = ", ifelse(is.na(s2[3]), "NA", sprintf("%.3f", s2[3])),
    "; Qᴇ = ",  ifelse(is.na(QE),  "NA", sprintf("%.3f", QE)),
    "; df = ",  ifelse(is.na(dfQ), "NA", sprintf("%.0f", dfQ)),
    "; p ",     ifelse(is.na(QEp), "= NA", paste0(ifelse(QEp < 0.001, "< ", "= "), fmt_p(QEp))),
    ")"
  )
}

coef_rows <- function(fit) {
  if (is.null(fit)) {
    return(tibble(
      Moderator_variables="(SKIP / model not fitted)",
      Estimate=NA_real_, SE=NA_real_, z=NA_real_, p=NA_character_,
      LCI=NA_real_, UCI=NA_real_, sig_CI=NA
    ))
  }
  
  term <- rownames(fit$b)
  est  <- as.numeric(fit$b)
  se   <- as.numeric(fit$se)
  z    <- as.numeric(fit$zval)
  pv   <- as.numeric(fit$pval)
  lci  <- as.numeric(fit$ci.lb)
  uci  <- as.numeric(fit$ci.ub)
  
  tibble(
    Moderator_variables = clean_term(term),
    Estimate = round(est, 3),
    SE       = round(se, 3),
    z        = round(z, 3),
    p        = fmt_p(pv),
    LCI      = round(lci, 3),
    UCI      = round(uci, 3),
    sig_CI   = ifelse(is.na(lci) | is.na(uci), NA, (lci > 0) | (uci < 0))
  )
}

block_one_model <- function(name, fit) {
  title <- labels[[name]]
  if (is.null(title) || is.na(title)) title <- name
  
  header <- tibble(
    Moderator_variables = title,
    Estimate=NA_real_, SE=NA_real_, z=NA_real_, p=NA_character_,
    LCI=NA_real_, UCI=NA_real_, sig_CI=NA
  )
  
  varline <- tibble(
    Moderator_variables = variance_line(fit),
    Estimate=NA_real_, SE=NA_real_, z=NA_real_, p=NA_character_,
    LCI=NA_real_, UCI=NA_real_, sig_CI=NA
  )
  
  blank <- tibble(
    Moderator_variables = "",
    Estimate=NA_real_, SE=NA_real_, z=NA_real_, p=NA_character_,
    LCI=NA_real_, UCI=NA_real_, sig_CI=NA
  )
  
  bind_rows(header, varline, coef_rows(fit), blank)
}

Table <- bind_rows(lapply(names(models), function(nm) block_one_model(nm, models[[nm]])))
write.csv(Table %>% select(-sig_CI), "Table_rm21_rm38_subgroup.csv", row.names = FALSE)
write.csv(Table, "Table_rm21_rm38_subgroup_with_sigCI.csv", row.names = FALSE)