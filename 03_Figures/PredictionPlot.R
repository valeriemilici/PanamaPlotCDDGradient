###Initialize Workspace --------------------------------------------------------
rm(list =ls())
#load packages
library(mgcv) #for fitting gam
library(marginaleffects) #for making predictions from gam
library(tidyverse) #tidy code
library(ggplot2) #plotting predictions
library(patchwork) #combining plots
library(mgcViz)
library(gratia) #for plotting the gams

load("02_Analyses/SpeciesSpecific/res_mod")

#plot the gam estimated splines for the main effects
spline_list <- list()

for(i in 1:length(res_mod)){
  pl <- draw(res_mod[[i]], select = 3, ) +
    labs(title = names(res_mod)[i]) +
    xlab("") +
    ylab("") +
    theme_classic() +
    theme(plot.caption = element_blank())
  
  spline_list[[i]] <- pl
}

for (i in 1:length(res_mod)) {
  
  # Get the vizmod for the current species
  vizmod <- getViz(res_mod[[i]], post = T, unconditional = T)
  pl <- plot(vizmod, nsim = 20, allTerms = T, select = 3) + 
    # Add confidence interval line/ fit line/ simulation line
    l_ciLine() + l_fitLine() + l_simLine() +
    #Add confidence interval bar/# Add fitted points
    l_ciBar() + l_fitPoints(size = 1) +  
    l_rug() +                             # Add rug plot
    # Add title with the name of the current species
    labs(title = names(res_mod)[i])       
  
  spline_list[[i]] <- pl
}

spline1 <-  wrap_plots(spline_list[1:10], ncol = 2) +
 plot_annotation(caption = "Conspecific Seedling Density",
                 theme = theme(plot.caption = element_text(size = 16,
                                                           hjust = 0.5)))
spline2 <-  wrap_plots(spline_list[11:20], ncol = 2) +
  plot_annotation(caption = "Conspecific Seedling Density",
                  theme = theme(plot.caption = element_text(size = 16,
                                                            hjust = 0.5)))
spline3 <-  wrap_plots(spline_list[21:30], ncol = 2) +
  plot_annotation(caption = "Conspecific Seedling Density",
                  theme = theme(plot.caption = element_text(size = 16,
                                                            hjust = 0.5)))
spline4 <-  wrap_plots(spline_list[31:40], ncol = 2) +
  plot_annotation(caption = "Conspecific Seedling Density",
                  theme = theme(plot.caption = element_text(size = 16,
                                                            hjust = 0.5)))
spline5 <-  wrap_plots(spline_list[41:50], ncol = 2) +
  plot_annotation(caption = "Conspecific Seedling Density",
                  theme = theme(plot.caption = element_text(size = 16,
                                                            hjust = 0.5)))
spline6 <-  wrap_plots(spline_list[51:60], ncol = 2) +
  plot_annotation(caption = "Conspecific Seedling Density",
                  theme = theme(plot.caption = element_text(size = 16,
                                                            hjust = 0.5)))
spline7 <-  wrap_plots(spline_list[61:70], ncol = 2) +
  plot_annotation(caption = "Conspecific Seedling Density",
                  theme = theme(plot.caption = element_text(size = 16,
                                                            hjust = 0.5)))

ggsave(plot = spline1, filename = "03_Figures/figures/sp.part1.png")
ggsave(plot = spline2, filename = "03_Figures/figures/sp.part2.png")
ggsave(plot = spline3, filename = "03_Figures/figures/sp.part3.png")
ggsave(plot = spline4, filename = "03_Figures/figures/sp.part4.png")
ggsave(plot = spline5, filename = "03_Figures/figures/sp.part5.png")
ggsave(plot = spline6, filename = "03_Figures/figures/sp.part6.png")
ggsave(plot = spline7, filename = "03_Figures/figures/sp.part7.png")
# Rapidly create prediction plots from res_mod output
plot_list <- list()
for (i in 1:length(res_mod)) {
  
inset <- names(res_mod)[i]
 
   p <- plot_predictions(res_mod[[i]], condition = 'con_dens.p', 
                   type = 'response' ) +
    labs(y = "Predicted Probability of Mortality",
         x = "Conspecific Seedling Density") +
    #ylim(c(-0.1, 1.5)) +
    ggtitle((inset))+
    theme_classic(12) +
    theme(plot.title=element_text(vjust = -3, hjust = 0.05, size = 10))
 
   plot_list[[i]] <- p
}

part1 <- wrap_plots(plot_list[1:10], ncol = 2) + plot_layout(axis_titles = "collect")
part2 <- wrap_plots(plot_list[11:20], ncol = 2) + plot_layout(axis_titles = "collect")
part3 <- wrap_plots(plot_list[21:30], ncol = 2) + plot_layout(axis_titles = "collect")
part4 <- wrap_plots(plot_list[31:40], ncol = 2) + plot_layout(axis_titles = "collect")
part5 <- wrap_plots(plot_list[41:50], ncol = 2) + plot_layout(axis_titles = "collect")
part6 <- wrap_plots(plot_list[51:60], ncol = 2) + plot_layout(axis_titles = "collect")
part7 <- wrap_plots(plot_list[61:70], ncol = 2) + plot_layout(axis_titles = "collect")

ggsave(plot = part1, filename = "03_Figures/figures/pp.part1.png")
ggsave(plot = part2, filename = "03_Figures/figures/pp.part2.png")
ggsave(plot = part3, filename = "03_Figures/figures/pp.part3.png")
ggsave(plot = part4, filename = "03_Figures/figures/pp.part4.png")
ggsave(plot = part5, filename = "03_Figures/figures/pp.part5.png")
ggsave(plot = part6, filename = "03_Figures/figures/pp.part6.png")
ggsave(plot = part7, filename = "03_Figures/figures/pp.part7.png")

ANAM_pred <- plot_predictions(res_mod[[1]], condition = 'con_dens.p', 
                 type = 'response' ) +
  labs(y = "Predicted Probability of Mortality",
       x = "Conspecific Seedling Density") +
  ylim(c(-0.1, 0.8)) +
  ggtitle(bquote("2007 mm"~yr^-1))+
  theme_classic(12) +
  theme(plot.title=element_text(vjust = -3, hjust = 0.05, size = 10))

BV_pred <- plot_predictions(res_mod[[2]], condition = 'con_dens.p', 
                         type = 'response') +
  labs(y = "Predicted Probability of Mortality",
       x = "Conspecific Seedling Density") +
  ylim(c(-0.1, 0.8)) +
  ggtitle(bquote("2595 mm"~yr^-1))+
  theme_classic(12) +
  theme(plot.title=element_text(vjust = -3, hjust = 0.05, size = 10))

Charco_pred <- plot_predictions(res_mod[[3]], condition = 'con_dens.p', 
                                type = 'response') +
  labs(y = "Predicted Probability of Mortality",
       x = "Conspecific Seedling Density") +
  ylim(c(-0.1, 0.8)) +
  ggtitle(bquote("2051 mm"~yr^-1))+
  theme_classic(12) +
  theme(plot.title=element_text(vjust = -3, hjust = 0.05, size = 10))

Met_pred <- plot_predictions(res_mod[[4]], condition = 'con_dens.p', 
                                type = 'response') +
  labs(y = "Predicted Probability of Mortality",
       x = "Conspecific Seedling Density") +
  ylim(c(-0.1, 0.8)) +
  ggtitle(bquote("1874 mm"~yr^-1))+
  theme_classic(12) +
  theme(plot.title=element_text(vjust = -3, hjust = 0.05, size = 10))

Oleo_pred <- plot_predictions(res_mod[[5]], condition = 'con_dens.p', 
                                type = 'response') +
  labs(y = "Predicted Probability of Mortality",
       x = "Conspecific Seedling Density") +
  ylim(c(-0.1, 0.8)) +
  ggtitle(bquote("2330 mm"~yr^-1))+
  theme_classic(12) +
  theme(plot.title=element_text(vjust = -3, hjust = 0.05, size = 10))

PP_pred <- plot_predictions(res_mod[[6]], condition = 'con_dens.p', 
                                type = 'response') +
  labs(y = "Predicted Probability of Mortality",
       x = "Conspecific Seedling Density") +
  ylim(c(-0.1, 0.8)) +
  ggtitle(bquote("1756 mm"~yr^-1))+
  theme_classic(12) +
  theme(plot.title=element_text(vjust = -3, hjust = 0.05, size = 10))

SR_pred <- plot_predictions(res_mod[[7]], condition = 'con_dens.p', 
                                type = 'response') +
  labs(y = "Predicted Probability of Mortality",
       x = "Conspecific Seedling Density") +
  ylim(c(-0.1, 0.8)) +
  ggtitle(bquote("3054 mm"~yr^-1))+
  theme_classic(12) +
  theme(plot.title=element_text(vjust = -3, hjust = 0.05, size = 10))

Sherman_pred <- plot_predictions(res_mod[[8]], condition = 'con_dens.p', 
                                type = 'response') +
  labs(y = "Predicted Probability of Mortality",
       x = "Conspecific Seedling Density") +
  ylim(c(-0.1, 0.8)) +
  ggtitle(bquote("3203 mm"~yr^-1))+
  theme_classic(12) +
  theme(plot.title=element_text(vjust = -3, hjust = 0.05, size = 10))

combo_pred <- PP_pred + Met_pred + ANAM_pred + Charco_pred +
Oleo_pred + BV_pred + SR_pred + Sherman_pred +
  plot_layout(axis_titles = "collect",ncol = 4)

ggsave(combo_pred, filename ="03_Figures/MS_Figs/Composite_Pred_Fig.png",
       height = 4, width = 8, units = "in")
