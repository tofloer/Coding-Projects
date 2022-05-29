cor_matrix <- round(cor(Boston),2)
mcor <- melt(cor_matrix)
corp <- ggplot(mcor, aes(x=Var1, y=Var2, z=value)) + 
          theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank()) + 
          geom_tile(aes(fill= value)) + 
          scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1), 
                               space = "Lab", name="Pearson\nCorrelation") + 
          geom_text(aes(Var2, Var1, label = value), color = "black", size = 3)
plot(corp)