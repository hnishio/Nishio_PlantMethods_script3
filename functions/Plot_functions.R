

##### Correlation between visual and predicted #####
cor_vis_pred <- function(df, title){
  
  cor <- cor.test(as.numeric(df$predicted), as.numeric(df$visual), method="pearson")
  cor_estimate <- round(cor$estimate, digits=2)
  mylabs <- parse(text=sprintf('r=="%.2f"', cor_estimate))
  if(cor$p.value<0.05){
    plabs <- expression(paste(italic(P), " < 0.05", sep=""))
  }else{plabs <- "n.s."}
  
  min_axes <- 0
  max_axes <- max(df$visual)*1.1
  
  g <- ggplot(df, aes(x=visual, y=predicted)) +
    geom_abline(intercept=0, slope=1, linetype="dashed", col="gray") +
    geom_point(size=0.5) +
    coord_cartesian(xlim=c(min_axes, max_axes), ylim=c(min_axes, max_axes), clip='off') +
    annotate("text", x=10, y=max_axes*0.95, 
             label=mylabs, size=7/ggplot2::.pt) +
    annotate("text", x=10, y=max_axes*0.82, 
             label=plabs, size=7/ggplot2::.pt) +
    theme_bw() +
    theme(plot.title=element_text(size=7),
          axis.title=element_text(size=7), 
          axis.text=element_text(size=7),
          plot.tag = element_text(size = 12, face = "bold")) + 
    labs(title=title, x="Visual", y="Statistical")
  return(g)
}





##### Line plot pairwise (x: Total reaction time, y: distance) #####
pairwise_start_dis <- function(df, i){
  
  title <- paste0("Cell ", i)
  data <- df[df$cell==i,]
  align <- as.data.frame(t(combn(x = data$chloroplast, m = 2)))
  names(align) <- c("from", "to")
  
  min_axis_y <- min(data$distance) - diff(range(data$distance))*0.1
  max_axis_y <- max(data$distance) + diff(range(data$distance))*0.1
  min_axis_x <- min(data$start_time)
  max_axis_x <- max(data$start_time)
  
  g <- ggplot(data, aes(x=start_time, y=distance)) +
    geom_point(size=0.8, alpha=0.5)
  
  no_negative <- 0
  for(j in 1:nrow(align)){
    x_from <- data[data$chloroplast == align$from[j],"start_time"]
    x_to <- data[data$chloroplast == align$to[j],"start_time"]
    y_from <- data[data$chloroplast == align$from[j],"distance"]
    y_to <- data[data$chloroplast == align$to[j],"distance"]
    if((y_to - y_from) / (x_to - x_from) >= 0){
      g <- g +
        geom_segment(x=x_from, y=y_from, xend=x_to, yend=y_to, size=0.1, col="black")
    }else{
      g <- g +
        geom_segment(x=x_from, y=y_from, xend=x_to, yend=y_to, size=0.1, col="red")
      no_negative <- no_negative + 1
    }
  }
  
  g <- g +  
    coord_cartesian(xlim=c(min_axis_x, max_axis_x), ylim=c(min_axis_y, max_axis_y), clip='off') +
    annotate("text", x=max_axis_x - (max_axis_x - min_axis_x)*0.25, y=min_axis_y + (max_axis_y - min_axis_y)*0.27, 
             label="Negative", size=7/ggplot2::.pt) +
    annotate("text", x=max_axis_x - (max_axis_x - min_axis_x)*0.2, y=min_axis_y + (max_axis_y - min_axis_y)*0.15, 
             label="slopes:", size=7/ggplot2::.pt) +
    annotate("text", x=max_axis_x - (max_axis_x - min_axis_x)*0.17, y=min_axis_y + (max_axis_y - min_axis_y)*0.03, 
             label=paste(round(no_negative/nrow(align)*100, 1), "%", sep=""), size=7/ggplot2::.pt) +
    theme_bw() +
    theme(plot.title = element_text(size=7, face = "bold"),
          axis.title=element_text(size=7), 
          axis.text=element_text(size=7),
          plot.tag = element_text(size = 12, face = "bold")) + 
    labs(title=title, x="Start time (min)")
  if(i==1 | i==4 | i==7){
    g <- g +  
      labs(y=expression(paste("Distance from microbeam (", mu, "m)", sep = "")))
  }else{
    g <- g +  
      theme(axis.title.y = element_blank())
  }
  
  return(g)
}





##### Line plot (x: Total reaction time, y: distance) for aligned chloroplasts #####
aligned_start_dis <- function(df, df_align, i){
  
  title <- paste0("Cell ", i)
  data <- df[df$cell==i,]
  align <- df_align[df_align$cell==i,]
  
  rem_chl <- base::setdiff(unique(c(align$from, align$to)), data$chl)
  if(length(rem_chl) != 0){
    align <- align[!(align$from==rem_chl | align$to==rem_chl),]
  }
  
  min_axis_y <- min(data$distance) - diff(range(data$distance))*0.1
  max_axis_y <- max(data$distance) + diff(range(data$distance))*0.1
  min_axis_x <- min(data$start_time)
  max_axis_x <- max(data$start_time)
  
  g <- ggplot(data, aes(x=start_time, y=distance)) +
    geom_point(size=0.8, alpha=0.5)
  
  no_negative <- 0
  for(j in 1:nrow(align)){
    x_from <- data[data$chloroplast == align$from[j],"start_time"]
    x_to <- data[data$chloroplast == align$to[j],"start_time"]
    y_from <- data[data$chloroplast == align$from[j],"distance"]
    y_to <- data[data$chloroplast == align$to[j],"distance"]
    if((y_to - y_from) / (x_to - x_from) >= 0){
      g <- g +
        geom_segment(x=x_from, y=y_from, xend=x_to, yend=y_to, size=0.1, col="black")
    }else{
      g <- g +
        geom_segment(x=x_from, y=y_from, xend=x_to, yend=y_to, size=0.1, col="red")
      no_negative <- no_negative + 1
    }
  }#for(j in 1:nrow(align)){
  
  g <- g +  
    coord_cartesian(xlim=c(min_axis_x, max_axis_x), ylim=c(min_axis_y, max_axis_y), clip='off') +
    annotate("text", x=max_axis_x - (max_axis_x - min_axis_x)*0.25, y=min_axis_y + (max_axis_y - min_axis_y)*0.27, 
             label="Negative", size=7/ggplot2::.pt) +
    annotate("text", x=max_axis_x - (max_axis_x - min_axis_x)*0.2, y=min_axis_y + (max_axis_y - min_axis_y)*0.15, 
             label="slopes:", size=7/ggplot2::.pt) +
    annotate("text", x=max_axis_x - (max_axis_x - min_axis_x)*0.17, y=min_axis_y + (max_axis_y - min_axis_y)*0.03, 
             label=paste(round(no_negative/nrow(align)*100, 1), "%", sep=""), size=7/ggplot2::.pt) +
    theme_bw() +
    theme(plot.title = element_text(size=7, face = "bold"),
          axis.title=element_text(size=7), 
          axis.text=element_text(size=7),
          plot.tag = element_text(size = 12, face = "bold")) + 
    labs(title=title, x="Start time (min)")
  if(i==1 | i==4 | i==7){
    g <- g +  
      labs(y=expression(paste("Distance from microbeam (", mu, "m)", sep = "")))
  }else{
    g <- g +  
      theme(axis.title.y = element_blank())
  }
  
  return(g)
}

