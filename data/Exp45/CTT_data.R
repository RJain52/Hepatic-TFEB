# TFEB KD CTT Data analysis
# Raghav Jain

library(tidyverse)
library(scales)
library(ggsci)
library(ggpubr)
library(rstatix)
library(Cairo)
show_col(pal_npg("nrc")(10))
?pal_npg

list.files()

mdat <- readxl::read_xlsx("Exp_45_data.xlsx", sheet = "Sheet2")

temp <- mdat |> select(`Con or KD`, `Exp ID #`, RT_Cold, `T-0h`:`T-6h`)
temp <- temp |> 
  group_by(`Con or KD`, RT_Cold) |>
  pivot_longer(cols = starts_with("T"),
               names_to = "Time", 
               values_to = "Temp")

temp_sum <- temp |> 
  group_by(`Con or KD`, RT_Cold, Time) |>
  summarise(
    mean = mean(Temp, na.rm = T),
    sd = sd(Temp,na.rm = T)
  )
  
temp_sum$group <- paste0(temp_sum$`Con or KD`, "_", temp_sum$RT_Cold)
temp_sum$group <- factor(temp_sum$group, 
                         levels = c("Con_RT", "Con_Cold",
                                    "KD_RT", "KD_Cold"))

p0 <- temp_sum |> ggplot(aes(x=Time, y=mean, group=group, color=group, shape =`Con or KD`)) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1, 
                position=position_dodge(0.00)) +
  geom_line() + geom_point(size = 3) +
  scale_color_manual(values = c("#F39B7FFF", "#8491B4FF",
                                "#DC0000FF", "#3C5488FF")) +
  scale_x_discrete(labels = c("0", "1", "2", "3", "4", "5", "6"), expand = c(0.05, 0.05)) +
  xlab("Time (h)") +
  ylab("Core temperature (C)") +
  ggtitle("Cold Tolerance Test", 
          subtitle = "Male B6 mice given AAV8 with eGFP-shRNA containing scramble or TFEB target sequence.") +
  theme_bw(base_size = 16, base_family = "Arial") + 
  theme(text = element_text(face = "bold"))
  

CairoPDF(file = "TFEB_KD_CTT.pdf", height = 8, width = 10)
  print(p0)
dev.off()

###################################################
  
liver <- mdat |> select(`Con or KD`, `Exp ID #`, RT_Cold, `Liver wet weight (g)`)
liver$`Con or KD` <- factor(liver$`Con or KD`, 
                           levels = c("Con", "KD"))
liver$RT_Cold <- factor(liver$RT_Cold, 
                            levels = c("RT", "Cold"))

liver$group <- paste0(liver$`Con or KD`, "_", liver$RT_Cold)
liver$group <- factor(liver$group, 
                         levels = c("Con_RT", "Con_Cold",
                                    "KD_RT", "KD_Cold"))
  
mypal <- c("#DC0000FF", "#3C5488FF")

p0 <- liver |> ggbarplot(x = "Con or KD", y = "Liver wet weight (g)", color = "RT_Cold",
                  fill = "RT_Cold", position = position_dodge(0.8), 
                  add.params = list(width = 0.35, 
                                    binwidth = 0.04),
                  alpha = 0.9, size = 0.5,
                  add = c("mean_sd", "dotplot"))


p1 <- ggpar(p0, palette = mypal, 
            legend.title = expression(bold("Temperature"))) +
  scale_y_continuous(limits = c(0.0, 1.6), 
                     expand=c(0,0)) +
  theme_bw(base_size = 16, base_family = "Arial") +
  theme(legend.position = "right", 
        text = element_text(face="bold", size = 14, 
                            family = "Arial", 
                            colour = "black")
  ) + 
  ylab("Wet weight (g)") + 
  xlab("AAV8 eGFP-shRNA target") +
  scale_x_discrete(label = c("scramble", "TFEB")) +
  ggtitle("Liver weight from mice at RT or cold for 6h", 
          subtitle = "Male B6 mice treated with AAV8 eGFP-scramble or TFEB target sequence.")
print(p1)

stat.test <- liver |>
  group_by(`Con or KD`) |>
  na.omit() %>%
  t_test(`Liver wet weight (g)` ~ RT_Cold) %>%
  adjust_pvalue(method = "none") %>%
  add_significance("p")

stat.test$p.adj.signif <- stat.test$p.signif
stat.test <- stat.test %>%
  na.omit() %>%
  add_xy_position(fun = "max", "Con or KD", dodge = 0.8) 

p2 <- p1 + stat_pvalue_manual(
  stat.test, label = "p", tip.length = 0.02, hide.ns = F)
p2

CairoPDF(file = "Liver_weight.pdf", height = 8, width = 8)
  print(p2)
dev.off()

##########################################################

liver$initial <- mdat$`initial weight (g)`

p0 <- liver |> ggbarplot(x = "Con or KD", y = "initial", color = "RT_Cold",
                         fill = "RT_Cold", position = position_dodge(0.8), 
                         add.params = list(width = 0.35, 
                                           binwidth = 0.80),
                         alpha = 0.9, size = 0.5,
                         add = c("mean_sd", "dotplot"))


p1 <- ggpar(p0, palette = mypal, 
            legend.title = expression(bold("Temperature"))) +
  scale_y_continuous(limits = c(0.0, 35), 
                     expand=c(0,0)) +
  theme_bw(base_size = 16, base_family = "Arial") +
  theme(legend.position = "right", 
        text = element_text(face="bold", size = 14, 
                            family = "Arial", 
                            colour = "black")
  ) + 
  ylab("Body weight (g)") + 
  xlab("AAV8 eGFP-shRNA target") +
  scale_x_discrete(label = c("scramble", "TFEB")) +
  ggtitle("Initial body weight of mice", 
          subtitle = "Male B6 mice treated with AAV8 eGFP-scramble or TFEB target sequence.")
print(p1)

stat.test <- liver |>
  group_by(`Con or KD`) |>
  na.omit() %>%
  t_test(initial ~ RT_Cold) %>%
  adjust_pvalue(method = "none") %>%
  add_significance("p")

stat.test$p.adj.signif <- stat.test$p.signif
stat.test <- stat.test %>%
  na.omit() %>%
  add_xy_position(fun = "max", "Con or KD", dodge = 0.8) 

p2 <- p1 + stat_pvalue_manual(bracket.nudge.y = 1.0,
  stat.test, label = "p", tip.length = 0.02, hide.ns = F)
p2

CairoPDF(file = "Initial_weight.pdf", height = 8, width = 8)
  print(p2)
dev.off()
##########################################################

liver$final <- mdat$`final weight (g)`
liver <- liver |> filter(final < 35)

p0 <- liver |> ggbarplot(x = "Con or KD", y = "final", color = "RT_Cold",
                         fill = "RT_Cold", position = position_dodge(0.8), 
                         add.params = list(width = 0.35, 
                                           binwidth = 0.80),
                         alpha = 0.9, size = 0.5,
                         add = c("mean_sd", "dotplot"))

p1 <- ggpar(p0, palette = mypal, 
            legend.title = expression(bold("Temperature"))) +
  scale_y_continuous(limits = c(0.0, 40), 
                     expand=c(0,0)) +
  theme_bw(base_size = 16, base_family = "Arial") +
  theme(legend.position = "right", 
        text = element_text(face="bold", size = 14, 
                            family = "Arial", 
                            colour = "black")
  ) + 
  ylab("Body weight (g)") + 
  xlab("AAV8 eGFP-shRNA target") +
  scale_x_discrete(label = c("scramble", "TFEB")) +
  ggtitle("Final body weight of mice at RT or cold for 6h", 
          subtitle = "Male B6 mice treated with AAV8 eGFP-scramble or TFEB target sequence.")
print(p1)

stat.test <- liver |>
  group_by(`Con or KD`) |>
  na.omit() %>%
  t_test(final ~ RT_Cold) %>%
  adjust_pvalue(method = "none") %>%
  add_significance("p")

stat.test$p.adj.signif <- stat.test$p.signif
stat.test <- stat.test %>%
  na.omit() %>%
  add_xy_position(fun = "max", "Con or KD", dodge = 0.8) 

p2 <- p1 + stat_pvalue_manual(bracket.nudge.y = 1.0,
                              stat.test, label = "p", 
                              tip.length = 0.02, hide.ns = F)
p2


CairoPDF(file = "final_weight.pdf", height = 8, width = 8)
  print(p2)
dev.off()
