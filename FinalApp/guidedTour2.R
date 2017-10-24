library(ggplot2)
library(plotly)
library(tourr)
library(crosstalk) 
library(shiny)
#library(htmltools)
library(tidyr) # For PCs reshaping output for plot
library(GGally)

# Load Auckland schoos data
load("data/akl.RData")

# Function for reordering first two cols (df can be sphere'd or not)
# All possible combinations of Xs as first two variables
# Output is a list with p(p-1)/2 components that are (n by p) data frames (ready for touring)
col_reorder <- function (df) {
  p <- ncol(df)
  Xnames <- colnames(df)
  Xs <- list()
  tour_names <- vector(mode="character", length=p*(p-1)/2)
  k <- 1
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      Xs[[k]] <- df[c(i,j,(1:p)[-c(i,j)])]
      tour_names[k] <- paste(Xnames[i], Xnames[j], sep="&")
      k <- k + 1
    }
  }
  reordered <- list(Xs, tour_names)
  return(reordered)
}

# Central mass index and projected data, XA
# 'base' is a projection matrices from the 3rd slice of array components in the list tinterp
# 'Xdf' is the scaled X matrix
cmass_tour <- function(base, Xdf) {
  X_matrix <- as.matrix(Xdf)
  XA <- X_matrix%*%matrix(base, ncol=2) # (n by d)
  cmass_index <- (sum(exp(-0.5*diag(XA%*%t(XA))))/dim(XA)[1]-exp(-dim(XA)[2]/2))/(1-exp(-dim(XA)[2]/2))
  list(rescale(XA), cmass_index)
}

# Holes index and projected data, XA
holes_tour <- function(base, Xdf) {
  X_matrix <- as.matrix(Xdf)
  XA <- X_matrix%*%matrix(base, ncol=2) # (n by d)
  holes_index <- (1-sum(exp(-0.5*diag(XA%*%t(XA))))/dim(XA)[1])/(1-exp(-dim(XA)[2]/2))
  list(rescale(XA), holes_index)
}

# 't_tour' fn produces a list of p(p-1)/2 lists with varying lengths (depend on the length of each tour)
# p(p-1)/2 is the number of combos for starting projections
# 't_tour' collect together data on the projection pursuit index,
# x and y co-ords for the tour projections and tour axes.
# The bases are components in the list `tinterp`
t_tour <- function(i, bases, Xdfs, index) {
  if (index=="cmass") {
    lapply(bases[[i]], cmass_tour, Xdfs[[i]])
  } else if (index=="holes") {
    lapply(bases[[i]], holes_tour, Xdfs[[i]])
  }
}

# Function to collect variable required for the 'indexPlot'
# ie. initial orthogonal projection used (init_pair), iteration step and pursuit index 
# Sub-function: Unlists the proj pursuit index for each tour (stored within t_tour)
pp_index <- function(j, a_tour) {unlist(a_tour[[j]][2])}
# Apply pursuit_info to the first component of:
# base = tinterp[[i]], t_name = tour_names[i], a_tour = t_tour[[i]]
pursuit_info <- function(base, t_name, a_tour) {
  # steps is the number of iterations for that tour
  steps <- length(base)
  init_pair <- rep(t_name, steps)
  iteration <- 1:steps
  pursuit_index <- do.call(c, lapply(1:steps, pp_index, a_tour))
  pp_info <- list(data.frame(init_pair, iteration, pursuit_index))
  return(pp_info)
}

# Extracts xy coords for "tour" projection plots of a single tour (ie. proj coords of a tour)
XYsingle <- function (single_tour) {
  n <- length(single_tour)
  #Below is equivalent? If so change AXsingle too.
  #x <- sapply(1:n, function(i) {unlist(single_tour[[i]][[1]][, 1])})
  x <- do.call(c, lapply(1:n, function(i) {unlist(single_tour[[i]][[1]][, 1])}))
  y <- do.call(c, lapply(1:n, function(i) {unlist(single_tour[[i]][[1]][, 2])}))
  XY <- data.frame(x=x, y=y, iteration = rep(1:n, each=length(x)/n))
}

# Extracts xy coords for "axes" plots of a single tour 
AXsingle <- function (single_tinterp) {
  n <- length(single_tinterp)
  p <- dim(single_tinterp)[1] # Number of real valued Xs
  # Take first 'p' values to be x-coords for axes
  AX_x <- do.call(c, lapply(1:n, function(i) {unlist(single_tinterp[[i]][1:p])}))
  # Take remaining 'p' values to be y-coords for axes
  AX_y <- do.call(c, lapply(1:n, function(i) {unlist(single_tinterp[[i]][(p+1):(p+p)])}))
  AX <- data.frame(x = AX_x, y = AX_y, iteration = rep(1:n, each = p),
                   measure=rep(colnames(attr(single_tinterp, "data")), n))
}

# Functions for extracting last projection of a (single) tour (for scattermatrix)
# Extract the last projection for a tour
last_XY <- function(XY) {
  m <- max(XY$iteration)
  XY[which(XY$iteration==m), -3]
}

# Extract the last axes plot coords for a tour
last_AX <- function(AX) {
  m <- max(AX$iteration)
  AX[which(AX$iteration==m), -3]
}

# Draws "tour" (proj) plot 
tour_plot <- function(df) {
  tx <- list(
    title="", range = c(-0.1, 1.2), 
    zeroline=F, showticklabels=F, showgrid=F
  )
  #All tour plots linked by brushing
  df %>% SharedData$new(key=~ID, group="2Dtour") %>% 
    plot_ly(x=~x, y=~y, frame=~iteration, color=I("black")) %>%
    add_markers(text=~ID, hoverinfo="text") %>%
    layout(xaxis=tx, yaxis=tx, 
           margin=list(l=0, r=0, b=0, t=0, pad=0))
}

# Draws tour "axes" plot
# k is the tour number from the "plotly_click"
tour_axes <- function(df, tour_cols, k) {
  ax <- list(
    title="", range=c(-1.1, 1.2), 
    zeroline=F, showticklabels=F, showgrid=F
  )
  plot_ly(df, x=~x, y=~y, frame=~iteration, hoverinfo="none") %>%
    add_segments(xend = 0, yend = 0, color = I(tour_cols[k]), size = I(1)) %>%
    add_text(text = ~measure, color=I("black")) %>%
    layout(xaxis = ax,
           yaxis = list(title="", range=c(-2.1, 2.2), 
                        zeroline=F, showticklabels=F, showgrid=F), 
           margin=list(l=0, r=0, b=0, t=0, pad=0))
}

# Draws LAST "tour" (proj) plot for scattermatrix
last_plot <- function(proj_list) {
  names(proj_list) <- c("x", "y", "col")
  df <- data.frame(x=proj_list$x, y=proj_list$y)
  ggplot(df, aes(x=x, y=y)) +
    geom_point(col=proj_list$col) +
    scale_x_continuous(limits = c(-0.1, 1.1)) +
    scale_y_continuous(limits = c(-0.1, 1.1)) +
    theme_void() +
    theme(legend.position = "none") +
    coord_fixed()
}

# Draws LAST tour proj's axes for scattermatrix
last_axis_plot <- function(axis_list) {
  names(axis_list) <- c("x", "y", "m", "col")
  df <- data.frame(x=axis_list$x, y=axis_list$y, measure=axis_list$m)
  ggplot(df, aes(x=x, y=y)) +
    geom_segment(aes(xend=0, yend=0), colour=axis_list$col) +
    geom_text(aes(label=measure)) +
    scale_x_continuous(limits = c(-1.1, 1.1)) +
    scale_y_continuous(limits = c(-1.1, 1.1)) +
    theme_void() +
    coord_fixed()
}

# Transform data frame to plot PCP
ggpcp2 <- function(df, a=1, sd.group=NULL) {
  # Transform df to a long data frame
  df$ID <- rownames(df)
  long <- gather(df, "Qualification", "Achievement", c("L1", "L2", "L3", "UE"))
  # Static PCP
  pcp.static <- long %>% SharedData$new(key=~ID, group=sd.group) %>%
    ggplot(aes(x=Qualification, y=Achievement, group=ID, colour=Decile)) +
    geom_line(alpha=a) +
    geom_point(alpha=a, size=0.01) +
    labs(x="Qualification level", y="Achievement rate") +
    theme(legend.position="none")
  return(pcp.static)
}

guidedTour_app <- function(dataset=akl[-1], index="cmass", factors=1, PC=TRUE, ...) {
  # Subset real-valued vars (type="double") for touring
  realXs <- dataset[, which(sapply(dataset, typeof)=="double")]
  rownames(realXs) <- rownames(dataset)
  Xdataset <- as.data.frame(scale(realXs)) # No need if using grand_tour
  if (PC) {
    Xdataset <- as.data.frame(apply(predict(prcomp(realXs, scale.=T)), 2, scale))
  }
  # All possible combinations of pairs of Xs as first two variables
  reordered <- col_reorder(Xdataset) 
  Xdfs <- reordered[[1]]  
  tour_names <- reordered[[2]] 
  # Names for each tour using the initial pair's names (eg."PC1&PC2"). 
  
  # Save tour (rescale=F since Xdfs already rescaled to range [0, 1])
  if (index=="cmass") {
    t <- lapply(Xdfs, function (x) save_history(x, guided_tour(cmass, d=2, max.tries = 50), rescale=FALSE, max=50, ...))
  } else if (index=="holes") {
    t <- lapply(Xdfs, function (x) save_history(x, guided_tour(holes, d=2, max.tries = 50), rescale=FALSE, max=50, ...))
  } else {
    stop("Invalid 'index' argument. Choose either: cmass or holes")
  }
  
  # 'Inserting' initial base as orthogonal projection of first two measures
  p <- ncol(Xdataset) # Number of real-valued Xs 
  t0 <- matrix(c(1, rep(0, p), 1, rep(0, (p-2))), ncol = 2)
  class(t0) <- "history_array"
  t1 <- lapply(t, function (x) array(c(t0, x), dim = dim(x) + c(0,0,1)))
  # Assign class and data attributes to new basis 
  for (i in 1:length(t1)) {
    class(t1[[i]]) <- "history_array"
    attr(t1[[i]], "data") <- Xdfs[[i]]
  }
  tinterp <- lapply(t1, interpolate)
  # t, t1, tinterp, t_tour, Xdfs are all lists with p(p-1)/2 components 
  # The components contain info for the tour resulting from each starting position
  
  # t_tour contains data for the 'tour' plot and its subplots 'axes' and 'tour'. 
  # t_tour contains data for the pursuit_index: for 'indexPlot'
  t_tour <- lapply(1:length(tinterp), t_tour, tinterp, Xdfs, index)
  
  # pp_info contains data for the 'indexPlot'
  pp_info <- do.call(rbind.data.frame, 
                     mapply(FUN=pursuit_info, 
                            base=tinterp, t_name=tour_names, a_tour=t_tour))
  
  # "tour" proj plot xy coords for ALL tours
  XYall <- lapply(t_tour, XYsingle) 
  # Last tour projection for each tour.
  last_projs <- lapply(XYall, last_XY) # No need for ID since not in SharedData
  
  # Add "ID" variable, rownames of original dataset
  XYall <- lapply(XYall, function (x) {
    ID <- rownames(dataset)
    cbind(x, ID)
  })
  
  # "axes" plot coords for ALL tours
  AXall <- lapply(tinterp, AXsingle)
  
  # Pull out last tour axes plot for each tour.
  last_axes <- lapply(AXall, last_AX)
  
  # Colour scale used in the index ggplot2 plot
  tour_cols <- hcl(h=seq(15, 360, 360/(p*(p-1)/2)), c=100, l=65)
  # Colour each last tour projection to match index plot cols
  if (p>5) {
    index_cols <- tour_cols[1:10] # Maximum of 10 pairwise combos
  } else {
    index_cols <- tour_cols 
  }
  
  # x, y, and col for each last proj is in a list
  last_projs <- mapply(c, last_projs, as.list(index_cols), SIMPLIFY = FALSE)
  last_axes <- mapply(c, last_axes, as.list(index_cols), SIMPLIFY = FALSE)
  
  P <- ifelse(p > 5, 5, p) # Maximum of 10 pairwise combos
  # Reorder plots for last_axes so that the scattermatrix will be diagonally symmetrical
  if (P==4) {
    last_axes <- list(last_axes[[1]], last_axes[[2]], last_axes[[4]],
                      last_axes[[3]], last_axes[[5]], last_axes[[6]])
  } else if (P==5) {
    last_axes <- list(last_axes[[1]], last_axes[[2]], last_axes[[5]],
                      last_axes[[3]], last_axes[[6]], last_axes[[8]],
                      last_axes[[4]], last_axes[[7]], last_axes[[9]],
                      last_axes[[10]])
  }
  
  # Scattermatrix plot data (max of 10 pairwise combos)
  # Lower diagonal is final proj and upper diagonal its axes components
  # List of P^2 (max 5) plots to plot using ggmatrix. 
  # Plots on diagonals (i,i) are empty
  diag_labels <- colnames(Xdfs[[1]])
  scatt_list <- list(ggally_text(diag_labels[1])) 
  tplot_list <- lapply(last_projs, last_plot)
  aplot_list <- lapply(last_axes, last_axis_plot)
  for (i in 1:(P-1)) {
    scatt_list <- c(scatt_list, tplot_list[1:(P-i)],
                    aplot_list[1:i], list(ggally_text(diag_labels[i+1]))) 
    # (P-i) is the number of tour plots and i the # of axes plots for iteration i
    # Delete merged plots
    tplot_list <- tplot_list[-c(1:(P-i))]
    aplot_list <- aplot_list[-c(1:i)]
  }
  # scatt_list contains all last projections (for each tour, up to 10 tours)
  
  # Subset categorical vars (for linked brushing by groups)
  Fcols <- which(sapply(dataset, class)=="factor")
  Fdataset <- as.data.frame(dataset[, Fcols])
  if (factors <= length(Fdataset)) {
    fac_n <- factors
  } else if (length(Fdataset) >= 2) {
    fac_n = 2 # Use default if 'factors' arg not sensible
    warning("'factors' argument exceeds number of factors in data.")
  } else {
    fac_n = length(Fdataset)
    warning("'factors' argument exceeds number of factors in data.")
  }
  colnames(Fdataset) <- attr(Fcols, "names")
  Fdataset$ID <- rownames(dataset)
  Fdataset$All <- factor(rep("1", nrow(dataset)))
  
  # SharedData for the factors plot
  sdF <- SharedData$new(Fdataset, key=~ID, group="2Dtour") 
  
  ## Shiny app
  ui <- fluidPage(
    splitLayout(cellWidths = c("25%", "45%", "30%"),
                plotlyOutput("indexPlot"),
                plotlyOutput("tour"),
                plotOutput("scatter")
    ),
    splitLayout(cellWidths = c("25%", "48%", "27%"),
                plotlyOutput("factors"), 
                plotlyOutput("pcp"),
                radioButtons("scaling", "Select axes scaling", 
                             c("Standardise", "Univariate_minmax", "Global_minmax"), selected="Standardise")
    )
  )
  server <- function(input, output, session) {
  # Proj pursuit index plot 
  output$indexPlot <- renderPlotly({
    ggplot(pp_info, aes(x=iteration, y=pursuit_index, group=init_pair, col=init_pair)) +
      geom_point() +
      geom_line() +
      labs(title="Pursuit index", subtitle="Click to select tour") +
      theme(legend.position="none") # Selecting by legend hard to reset.
      ggplotly(source="index_plot", tooltip = c("x", "y", "group")) %>%
        layout(autosize=F, width=400, height=400)
  })
  
  # PCP plot
  output$pcp <- renderPlotly({
    if (input$scaling == "Standardise") {
      nzqa.scale <- cbind(dataset[5], scale(dataset[1:4]))
    } else if (input$scaling == "Univariate_minmax") {
      nzqa.scale <- cbind(dataset[5], rescale(dataset[1:4]))
    } else {
      nzqa.scale <- dataset
    }
    static <- ggpcp2(nzqa.scale, sd.group="2Dtour") +
      ggtitle("Achievement of Auckland schools in 2016") +
      labs(x="Qualification level", y="Achievement")
    ggplotly(static, tooltip=c("group", "colour", "x", "y")) %>%
      highlight(on="plotly_click", off="plotly_doubleclick", color="red", persistent=T) %>%
      layout(autosize=F, width=800, height=400, showlegend=FALSE)
  })
  
  output$tour <- renderPlotly({
    s <- event_data("plotly_click", source = "index_plot")
    if (length(s)) {k <- s$curveNumber + 1} else {k <- 1}
    proj <- XYall[[k]]
    basis <- AXall[[k]]
    # tour plot
    tour <- tour_plot(proj)
    axes <- tour_axes(basis, tour_cols, k)
    subplot(tour, axes, titleY=F, widths=c(0.7, 0.3), margin=0) %>%
      animation_opts(33) %>% #33 milliseconds between frames
      hide_legend() %>%
      layout(dragmode="select",
             margin=list(l=0, r=0, b=0, t=0, pad=0)) %>%
      highlight(on="plotly_select", off="plotly_deselect", color="blue", 
                persistent=T, dynamic=T)
  })
  
  output$factors <- renderPlotly({
    # facet_grid allows up to 4 factors to be crossed into groups
    if (fac_n==1) {
      gp <- ggplot(sdF, aes_string(x=names(Fdataset)[1])) +
        geom_jitter(aes(y=Fdataset$All, label=ID), 
                    width=0, height=0.25) +
        labs(title="The first factor is displayed") 
    } else if (fac_n==2) {
      gp <- ggplot(sdF, aes_string(x=names(Fdataset)[1], y=names(Fdataset)[2])) +
        geom_jitter(aes(label=ID), 
                    width = 0.25, height = 0.25) +
        labs(title="The first 2 factors are displayed") 
    } else if (fac_n==3) {
      gp <- ggplot(sdF, aes_string(x=names(Fdataset)[1], y=names(Fdataset)[2],
                                  label1=names(Fdataset)[3])) +
        geom_jitter(aes(label=ID), 
                    width=0.2, height=0.2) +
        labs(title="The first 3 factors are displayed") +
        facet_grid(.~Fdataset[,3])
    } else if (fac_n>3) {
      gp <- ggplot(sdF, aes_string(x=names(Fdataset)[1], y=names(Fdataset)[2],
                                  label1=names(Fdataset)[3],
                                  label2=names(Fdataset)[4])) +
        geom_jitter(aes(label=ID), 
                    width = 0.2, height = 0.2) +
        labs(title="The first 4 factors are displayed") +
        facet_grid(Fdataset[,4]~Fdataset[,3], scales="free", space="free")
    }
    if (fac_n==0) { # fac_n=0 No factors in dataset
      plotly_empty()
      warning("There are no variables that are factors in the data set.")
    } else {
      ggplotly(gp, tooltip=c("x", "y", "label", "label1", "label2"), height=400, width=400) %>%
        layout(dragmode="select") %>%
        hide_legend() %>%
        highlight(on="plotly_select", off="plotly_deselect", color="blue", 
                  persistent=T)
    }
  })
  
  output$scatter <- renderPlot({
    ggmatrix(scatt_list, P, P, byrow = F, showAxisPlotLabels = F,
             title="Final tour projection and axes plots")
  })
  
  }
  shinyApp(ui, server)
}

guidedTour_app()
