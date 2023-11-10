simple_heatmap <- function(x){
  library(reshape2)
  melted_x <- melt(x)
  min_x = min(x)
  mid_x = mean(x)
  max_x = max(x)
  
  head(melted_x)
  library(ggplot2)
  p = ggplot(data = melted_x, aes(x=Var1, y=reorder(Var2, -value), fill=value)) + 
    geom_tile(color = "grey")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = mid_x, limit = c(min_x, max_x), space = "Lab", 
                         name="Scale") +
    theme_light()+ 
    theme(axis.text.y = element_text(size = 9))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                     size = 8, hjust = 1))+
    # theme(axis.text.x = element_blank())+
    coord_fixed()
   return(p)
}
