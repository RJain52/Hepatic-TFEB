library(tidyverse)
library(ggpubr)
library(rstatix)
library(Cairo)
list.files()
mydat <- read_csv("2022_0408_CathB_R.csv")

mydat$Group <- factor(mydat$Temp, levels = c("RT", "Cold"))

mydat %>% ggplot(aes(x = Group, y = Fluorescence_ugProt, fill = Group)) 

p0 <- ggbarplot(mydat, x = "Group", y = "Fluorescence_ugProt", 
                add = c("mean_se", "dotplot"), color = "Group", fill = "Group", 
                add.params = list(width = 0.35, 
                                  binwidth = .08*min(mydat$Fluorescence_ugProt, na.rm = T)),
                alpha = 0.8, position = position_dodge(0.8), size = 0.5) +
      scale_y_continuous(limits = c(0.00, 1.2*max(mydat$Fluorescence_ugProt, na.rm = T)),
                     expand = c(0,0))

p1 <- ggpar(p0, palette = "npg", legend = "right", legend.title = "Temperature", 
            title = "Cathepsin B activity assay (Exp 33)",  
            xlab = "Temperature", ylab = "Fluorescent Units (normalized to protein)") 

stat.test <- mydat %>%
  na.omit() %>%
  t_test(Fluorescence_ugProt ~ Group) %>%
  adjust_pvalue(method = "none") %>%
  add_significance("p")

stat.test$p.adj.signif <- stat.test$p.signif
stat.test <- stat.test %>%
  na.omit() %>%
  add_xy_position(fun = "max", "Group", dodge = 0.8) 

p2 <- p1 + stat_pvalue_manual(
  stat.test, label = "p", tip.length = 0.02, hide.ns = F)

p3 <- p2 + theme_bw(base_size = 16) + theme(axis.text.x = element_text(angle = 45, hjust = 1.0))

CairoPDF(file = "2022_0201_CathepsinB.pdf")
  print(p3)
dev.off()
