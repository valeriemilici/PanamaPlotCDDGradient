###Initialize Workspace --------------------------------------------------------
rm(list =ls())

library(tidyverse)
library(ggplot2)
library(patchwork)


load("02_Analyses/SiteSpecific/res_mod")

  vizmod_1 <- getViz(res_mod[[1]], post = T, unconditional = T)
  vizmod_2 <- getViz(res_mod[[2]], post = T, unconditional = T)
  vizmod_3 <- getViz(res_mod[[3]], post = T, unconditional = T)
  vizmod_4 <- getViz(res_mod[[4]], post = T, unconditional = T)
  vizmod_5 <- getViz(res_mod[[5]], post = T, unconditional = T)
  vizmod_6 <- getViz(res_mod[[6]], post = T, unconditional = T)
  vizmod_7 <- getViz(res_mod[[7]], post = T, unconditional = T)
  vizmod_8 <- getViz(res_mod[[8]], post = T, unconditional = T)
 

  pl1 <- plot(vizmod_1, ylim = c(-3,3), nsim = 1000, select = 3) + 
    # Add confidence interval line/ fit line/ simulation line
    l_ciLine() + l_fitLine() + l_simLine() +
    #Add confidence interval bar/# Add fitted points
    l_ciBar() + l_fitPoints(size = 1) +  
    l_rug() +
    ylim(c(-2,3)) +
    labs(title = names(res_mod)[[1]]) 
    
  
  pl2 <- plot(vizmod_2, ylim = c(-3,3), nsim = 1000, select = 3) + 
    # Add confidence interval line/ fit line/ simulation line
    l_ciLine() + l_fitLine() + l_simLine() +
    #Add confidence interval bar/# Add fitted points
    l_ciBar() + l_fitPoints(size = 1) +  
    l_rug() +
    ylim(c(-2,3)) +
    labs(title = names(res_mod)[[2]]) 

  pl3 <- plot(vizmod_3, ylim = c(-3,3), nsim = 1000, select = 3) + 
    # Add confidence interval line/ fit line/ simulation line
    l_ciLine() + l_fitLine() + l_simLine() +
    #Add confidence interval bar/# Add fitted points
    l_ciBar() + l_fitPoints(size = 1) +  
    l_rug() +
    ylim(c(-2,3)) +
    labs(title = names(res_mod)[[3]]) 
    
  
  pl4 <- plot(vizmod_4, ylim = c(-3,3), nsim = 1000, select = 3) + 
    # Add confidence interval line/ fit line/ simulation line
    l_ciLine() + l_fitLine() + l_simLine() +
    #Add confidence interval bar/# Add fitted points
    l_ciBar() + l_fitPoints(size = 1) +  
    l_rug() +
    ylim(c(-2,3)) +
    labs(title = names(res_mod)[[4]]) 
   
  pl5 <- plot(vizmod_5, ylim = c(-3,3), nsim = 1000, select = 3) + 
    # Add confidence interval line/ fit line/ simulation line
    l_ciLine() + l_fitLine() + l_simLine() +
    #Add confidence interval bar/# Add fitted points
    l_ciBar() + l_fitPoints(size = 1) +  
    l_rug() +
    ylim(c(-2,3)) +
    labs(title = names(res_mod)[[5]]) 
    
  
  pl6 <- plot(vizmod_6, ylim = c(-3,3), nsim = 1000, select = 3) + 
    # Add confidence interval line/ fit line/ simulation line
    l_ciLine() + l_fitLine() + l_simLine() +
    #Add confidence interval bar/# Add fitted points
    l_ciBar() + l_fitPoints(size = 1) +  
    l_rug() +
    ylim(c(-2,3)) +
    labs(title = names(res_mod)[[6]]) 
     
  
  pl7 <- plot(vizmod_7, ylim = c(-3,3), nsim = 1000, select = 3) + 
    # Add confidence interval line/ fit line/ simulation line
    l_ciLine() + l_fitLine() + l_simLine() +
    #Add confidence interval bar/# Add fitted points
    l_ciBar() + l_fitPoints(size = 1) +  
    l_rug() +
    ylim(c(-2,3)) +
    labs(title = names(res_mod)[[7]]) 
  
   pl8 <- plot(vizmod_8, ylim = c(-3,3), nsim = 1000, select = 3) + 
    # Add confidence interval line/ fit line/ simulation line
    l_ciLine() + l_fitLine() + l_simLine() +
    #Add confidence interval bar/# Add fitted points
    l_ciBar() + l_fitPoints(size = 1) +  
    l_rug() +
    ylim(c(-2,3)) +
    labs(title = names(res_mod)[[8]]) 
     

# I have had to manually save these images, neither ggsave nor save work. 
  #patchwork also won't combine the plots although the internet says it should
  #very frustrating. 
