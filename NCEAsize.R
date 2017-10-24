library(plotly)
library(ggplot2)
library(crosstalk)

# Load merge school data
load("files/nzqa.sch.RData")
nzqa.sch$Decile <- as.factor(nzqa.sch$Decile)

# Function for PCA plot
ncea_PCA <- function(df, scale.factor=20, nudge=-0.2, decile.pch=T, gp="Decile", sd.group=NULL) {
  # Perform PCA
  pca <- prcomp(df[2:5], scale.=T)
  # Explained variance
  exp.var <- round(pca$sdev^2/sum(pca$sdev^2)*100, 1)
  # Scores (rotated principal components)
  scores <- data.frame(t(t(pca$x[, 1:2])*(pca$sdev[1:2]^-1)), df[-(2:5)], Group=df[, gp])
  scores$Decile <- as.factor(df$Decile)
  axis.limit <- c(-round(max(abs(scores[1:2])), 1), round(max(abs(scores[1:2])), 1))
  # Set up for linked brushing 
  sdPCA <- SharedData$new(scores, key=~School, group=sd.group)
  # Variable axes (loadings)
  lambda <- pca$sdev[1:2]/sqrt(nrow(df)-1)
  loadings <- data.frame(varnames = rownames(pca$rotation),
                         t(t(pca$rotation[, 1:2])*lambda*scale.factor))
  # Plot of first two PCs
  p <- ggplot(loadings) + 
    geom_segment(aes(x=0, y=0, xend=PC1, yend=PC2), colour="grey35") +
    geom_text(aes(x=PC1, y=PC2, label=varnames), colour="grey35", nudge_x=nudge) +
    labs(x=paste("PC1 (", exp.var[1], "% explained var.)", sep=""),
         y=paste("PC2 (", exp.var[2], "% explained var.)", sep="")) +
    scale_x_continuous(limits=axis.limit) +
    scale_y_continuous(limits=axis.limit*1.1) +
    theme(
      panel.border = element_rect(colour="grey", fill=NA),
      panel.background = element_blank(),
      legend.position = "none",
      axis.text=element_blank(),
      axis.ticks=element_blank())
  # Plot decile as character or colour by another variable
  if (decile.pch) {
    plot <- p + geom_text(data=sdPCA, aes(x=PC1, y=PC2, label=Decile, colour=Decile, label1=School))
  } else {
    plot <- p + geom_point(data=sdPCA, aes(x=PC1, y=PC2, colour=Group, label1=School))
  }
  return(plot)
}

# Factor plot to link and brush with
ncea_factor <- function(df, gp="Decile", sd.group="AKLncea") {
  # Set up for linked brushing 
  sdGroup <- SharedData$new(data.frame(df, Group=df[, gp]), key=~School, group=sd.group) 
  # Static plot
  factor.p <- ggplot(sdGroup) +
    geom_jitter(aes(y=Group, x=1, colour=Group, label=School), width=0.1, height=0) +
    labs(x=as.character(gp), y="") +
    scale_x_continuous(limits=c(0.7, 1.3)) +
    theme(
      panel.border=element_rect(colour="grey", fill=NA),
      panel.background=element_blank(),
      legend.position="none",
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank())
  # Interactive linked plot
  link.factor <- ggplotly(factor.p, tooltip=c("colour", "label"))
  return(link.factor)
}

ui <- fluidPage(
  inputPanel(
    selectInput("group", "Select variable", choices=as.list(colnames(nzqa.sch)[6:15]), selected="Decile"),
    div(style="width: 600px;", sliderInput("cohort2", label="Minimum size of Year 11, 12, 13 cohorts", min=0, max=100, step=10, value=0, animate=T))
  ),
  mainPanel(plotlyOutput("plots"))
)

server <- function(input, output) {
  output$plots <- renderPlotly({
    # Unique SharedData group depending on variable chosen
    sd_group <- paste(input$group, "ncea", sep="")
    # Static PCA plot
    PCA <- ncea_PCA(subset(nzqa.sch, Min.Cohort>=input$cohort2), decile.pch=F, gp=input$group, scale.factor=50, sd.group=sd_group) 
    
    # Link PCA plot
    link.PCA <- ggplotly(PCA, tooltip=c("label1"))
    
    # Link group plot
    link.group <- ncea_factor(subset(nzqa.sch, Min.Cohort>=input$cohort2), gp=input$group, sd.group=sd_group) %>%
      layout(yaxis=list(side="right"), margin=list(r=150))
    
    subplot(link.PCA, link.group, widths=c(0.5, 0.5), margin=0, titleX=T, titleY=T) %>%
      layout(dragmode="select", autosize = F, width=1000, height=500) %>%
      highlight(on="plotly_select", off="plotly_deselect", color="red", persistent=T)
  })
}

shinyApp(ui, server)