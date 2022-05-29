PCRplot <- CV_PCR$results
PCRplot <- as.data.frame(PCRplot$RMSE^2)
names(PCRplot)[1] <- "MSE"
PCRp <- ggplot(PCRplot, aes(y=MSE, x=as.numeric(row.names(PCRplot)))) + theme_bw() + 
          theme(plot.title = element_text(hjust = 0.5)) + 
          labs(y="PRESS", x = "Number of Components") + 
          geom_vline(xintercept=N, linetype="dashed") +
          theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16),  text = element_text(size=18),
                axis.title.y = element_text(margin = margin(t = 0, r = 12.5, b = 0, l = 0)),
                axis.title.x = element_text(margin = margin(t = 12.5, r = 0, b = 0, l = 0)),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                panel.background = element_rect(colour = "black", size=0.5)) + 
          geom_line()+ geom_point(color='red3',size=1.5) + 
          scale_x_continuous(breaks = round(seq(0, 140, by = 10),1)) + 
          scale_y_continuous(breaks = round(seq(0, 40, by = 5),1)) + 
          coord_cartesian(xlim=c(10, 107.5),ylim=c(10,30))
plot(PCRp)