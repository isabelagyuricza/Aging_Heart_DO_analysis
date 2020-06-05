################################################################################

# This function changes the panel's colors to CC the colors
# when using facet plots according to haplotypes. 

################################################################################

# It requires ggplot2 and qtl2 libraries.

# The order of haplotypes used for the facets should be in the correct order 
# (A, B, C, ... , H).

# The function takes as arguments the plot and the number of columns the user
# wants the plot to be split on (same number of columns as facet_wrap or facet_grid).
# The options are 8, 4, 1 or 1. 

# If the user didn't specify the number of columns in the facet_grid on ggplot, 
# the ncol should be 8. 


facet_CC_colour <- function(plot = plot, ncol = c(8, 4, 2, 1)){
  
  g <- ggplot_gtable(ggplot_build(plot))
  
  data(CCcolors)
  
  fills <- CCcolors
  
  if(ncol == 8) {
    
    stript <- which(grepl('strip-t', g$layout$name))
    
  } else if (ncol == 4) {
    
    stript <- which(grepl('strip-t', g$layout$name))
    
    stript <- c(stript[5:8],strip[1:4])
    
  } else if (ncol == 2) {
    
    stript <- which(grepl('strip-t', g$layout$name))
    
    stript <- c(stript[c(7,8)],stript[c(5,6)],stript[c(3,4)],stript[c(1,2)])
    
  } else { 
    
    stript <- which(grepl('strip-t', g$layout$name))
    
    stript <- stript[8:1]
    
    
  }
    
    k <- 1
    for (i in stript) {
      j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
      g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
      k <- k+1
    }

  return(g)
  
}


# How to use it:

# final_df %>% 
#   ggplot(aes(x=Age, y=estimate)) +
#   geom_point() +
#   geom_line(aes(group=Founder)) +
#   facet_wrap(~Founder, ncol=4) +
#   geom_errorbar(aes(x=Age, ymin=estimate-se, ymax=estimate+se), width=0.2) +
#   ggtitle("test") -> p
# 
# p <- facet_CC_colour(plot = p, ncol = 4)
# 
# grid.draw(p)
# 
# 
# final_df %>% 
#   ggplot(aes(x=Age, y=estimate)) +
#   geom_point() +
#   geom_line(aes(group=biotype)) +
#   facet_grid(biotype ~ Founder, scales = "free")+
#   geom_errorbar(aes(x=Age, ymin=estimate-se, ymax=estimate+se), width=0.2) +
#   theme_bw() + ggtitle("test") -> p
# 
# p <- facet_CC_colour(plot = p, ncol = 8)
# 
# grid.draw(p)
