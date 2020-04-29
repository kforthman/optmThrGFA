#' A function to plot the circle plots.
#'

Circbar <- function(mydata, ebar, graphtitle, textlabel, minx, maxx){

  data <- mydata
  labelsize <- textlabel

  # Assign min and max

  mymin <- ifelse(-1.5 + min(mydata$Lower) < -2,-1.5 + min(mydata$Lower),-2)
  mymax <- ifelse(1 + max(mydata$Upper) < 1, 1, 1 + max(mydata$Upper))

  # Set the level lines
  if (is.null(minx)){
    ifelse(min(mydata$Lower) < -.1, minx <- min(mydata$Lower), minx <- -.1)
    ifelse(min(mydata$Upper) < minx, minx <- min(mydata$Upper), minx <- minx)
  }
  if(is.null(maxx)){
    ifelse(max(mydata$Upper) > .1, maxx <- max(mydata$Upper), maxx <- .1)
    ifelse(max(mydata$Lower) > maxx, maxx <- max(mydata$Lower), maxx <- maxx)
  }

  # Set a number of 'empty bar' to add at the end of each group
  empty_bar <- ebar

  to_add = data.frame( matrix(NA, empty_bar*nlevels(data$GFAgroups), ncol(data)) )
  colnames(to_add) = colnames(data)
  to_add$GFAgroups=rep(levels(data$GFAgroups), each=empty_bar)
  data=rbind(data, to_add)
  data=data %>% arrange(GFAgroups)
  data$id=seq(1, nrow(data))

  # Get the name and the y position of each label
  label_data=data
  number_of_bar=nrow(label_data)
  angle= 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
  label_data$hjust<-ifelse( angle < -90, 1, 0)
  label_data$angle<-ifelse(angle < -90, angle+180, angle)

  # prepare a data frame for base lines

  base_data=data %>%
    group_by(GFAgroups) %>%
    dplyr::summarize(start=min(id), end=max(id) - empty_bar) %>%
    rowwise() %>%
    mutate(Title=mean(c(start, end)))

  # prepare a data frame for grid (scales)
  grid_data = base_data
  grid_data$end = grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
  grid_data$start = grid_data$start - 1
  #grid_data=grid_data[-1,]

  # Make the plot
  p = ggplot(data, aes(x=as.factor(id), y=Median, fill=GFAgroups)) +
    # Note that id is a factor. If x is numeric, there is some space between the first bar

    geom_bar(aes(x=as.factor(id), y=Median, fill=GFAgroups), stat="identity", alpha=0.5)  +
    geom_errorbar(aes(x=as.factor(id),ymin=Lower,ymax=Upper,color="Gray")) +

    # Add a level lines.
    geom_segment(data=grid_data, aes(x = end, y = minx, xend = start, yend = minx), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = maxx, xend = start, yend = maxx), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE )

  # Add text showing the value of each level lines
  p = p + annotate("text", x = rep(max(data$id),3), y = c(minx,0, maxx), label = format(c(minx,0, maxx),digits=2) , color="black", size=3 , angle=0, fontface="bold", hjust=1) +

    ylim(mymin,mymax) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(-1,4), "cm")
    ) +
    coord_polar() +
    geom_text(data=label_data, aes(x=id, y=mymax-.9, label=GFAvarlabs, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=labelsize, angle= label_data$angle, inherit.aes = FALSE )  +

    # Add base line information
    geom_segment(data=base_data, aes(x = start, y = -1.1 , xend = end, yend = -1.1, colour=GFAgroups), alpha=0.8, size=1 , inherit.aes = FALSE )  +

    geom_segment(data=base_data, aes(x = start, y = 0 , xend = end, yend = 0), colour = "Gray", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +

    geom_text(data=grid_data, aes(x = Title, y = -1.3, label=GFAgroups),color="black", fontface="bold") +

    annotate("text", x = 0, y = -1.8, label = c(graphtitle) , color="red", size=4 , fontface="bold")

  return(p)

}
