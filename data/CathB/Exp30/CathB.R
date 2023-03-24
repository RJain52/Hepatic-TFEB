library(tidyverse)
library(ggpubr)
library(rstatix)
library(Cairo)
list.files()
mydat <- read_csv("2022_0201_CathB_Exp30_r.csv")

mydat$Group <- factor(mydat$Group, levels = c("RT", "Cold"))
mydat$Treatment <- factor(mydat$Treatment, 
                          levels = c("no_prot", "prot",
                                     "prot+E64D"))
mydat$Norm_final_ini <- mydat$Norm_final
mydat <- mydat |> group_by(Sample) |>
  mutate(Norm_final = Norm_final_ini/Norm_final_ini[Treatment == 'no_prot'])

mydat %>% ggplot(aes(x = Treatment, y = Norm_final, fill = Group)) 

p0 <- ggbarplot(mydat, x = "Treatment", y = "Norm_final", 
                add = c("mean_se", "dotplot"), color = "Group", fill = "Group", 
                add.params = list(width = 0.35, 
                                  binwidth = .2*min(mydat$Norm_final, na.rm = T)),
                alpha = 0.8, position = position_dodge(0.8), size = 0.5) +
      scale_y_continuous(limits = c(0.00, 1.2*max(mydat$Norm_final, na.rm = T)),
                     expand = c(0,0)) +
      scale_x_discrete(labels = c("- Protein",
                                  "+ Protein",
                                  "+ Protein + inhibitor"))

p1 <- ggpar(p0, palette = "npg", legend = "top", legend.title = "Temperature", 
            title = "Cathepsin B activity assay",  
            xlab = "Treatment", ylab = "Fluorescent Units (normalized to RT)") +
        theme_bw(base_size = 18, base_family = 'Arial') +
        theme(text = element_text(colour = 'black', 
                                  family = "Arial", 
                                  size = 18, 
                                  face = "bold"),
              axis.text.y = element_text(hjust = -0.5)
              )

stat.test <- mydat %>%
  group_by(Treatment) %>%
  na.omit() %>%
  t_test(Norm_final ~ Group) %>%
  adjust_pvalue(method = "none") %>%
  add_significance("p")

stat.test$p.adj.signif <- stat.test$p.signif
stat.test <- stat.test %>%
  na.omit() %>%
  add_xy_position(fun = "max", "Treatment", dodge = 0.8) 

p2 <- p1 + stat_pvalue_manual(
  stat.test, label = "p.signif", 
  tip.length = 0.02, size = 18,
  hide.ns = F)

p3 <- p2 # + theme_bw(base_size = 16) + theme(axis.text.x = element_text(angle = 45, hjust = 1.0))

CairoPDF(file = "2022_1107_CathepsinB_norm_noProt.pdf", height = 10, width = 10)
  print(p3)
dev.off()
