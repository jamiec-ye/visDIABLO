#' Visualize
#' @author Jamie C. Ye
#
#' @param model A DIABLO model object.
#' @param featureMapping A list of data frames containing 'Data.Names', 'Gene.Symbols', 'Display.Names' for each of the data blocks.
#' @import shiny
#' @import shinythemes
#' @import shinydashboard
#' @import shinyBS
#' @import DT
#' @import igraph
#' @import visNetwork
#' @import plotly
#' @import ggmixOmics
#' @import network
#' @import sna
#' @import ggnetwork
#' @import ggplot2
#' @import tidyverse
#' @export

visualize <- function(model, rename = F, featureMapping = NULL) {
  if (rename == T){
    M <- convertHGNC(model)
  } else {
    M <- model
  }
  model1 <- M
  model2 <- M

  if (class(featureMapping) != "list"){
    featureMapping = NULL
    warning(sprintf("featureMapping parameter is not in correct format"))
  }

  # Get component names ----
  dataNames <- names(M$X)
  nEntries <- length(dataNames)
  nComp <- unique(M$ncomp)

  # Params ----
  geneEnrichment <- TRUE
  PPIIntegration <- FALSE

  # UI ----

  quarterWidth <- 3
  halfWidth <- 6
  tQuarterWidth <- 9
  fullWidth <- 12

  header <- dashboardHeader(title = "Dashboard"
  )
  # temp
  library(shinyBS)

  sidebar <- dashboardSidebar(
    sidebarMenu(
      # menuItem("DIABLO", tabName = "diablo"),
      # menuItem("Components", tabName = "components"),
      # menuItem("Variables", tabName = "variables"),
      menuItem(
        div(
          div(
            # edit1
            style="width:75%; display:inline-block; vertical-align: middle;",
            "BiPlot"
          ),
          div(
            # edit2
            style="display:inline-block; vertical-align: middle;",
            shinyBS::bsButton("q1", label = "", icon = icon("question"),
                     style = "info", size = "small"),
            shinyBS::bsPopover(id = "q1", title = "BiPlot",
                      content = paste0("A biplot is plot which aims to represent both the observations and variables of a matrix of multivariate data on the same plot."),
                      placement = "right",
                      trigger = "click",
                      options = list(container = "body")
            )
          )
        )
        , tabName = "biplot"),
      # menuItem("Loading Vectors", tabName = "loadingVectors"),
      # menuItem("Heatmap", tabName = "heatmap"),
      menuItem(
        div(
          div(
            # edit1
            style="width:75%; display:inline-block; vertical-align: middle;",
            "Network"
          ),
          div(
            # edit2
            style="display:inline-block; vertical-align: middle;",
            shinyBS::bsButton("q2", label = "", icon = icon("question"),
                     style = "info", size = "small"),
            shinyBS::bsPopover(id = "q2", title = "Network",
                      content = paste0("Lasso a group of nodes to perform geneset enrichment analysis"),
                      placement = "right",
                      trigger = "click",
                      options = list(container = "body")
            )
          )
        )
        , tabName = "network")
      # menuItem("Circos", tabName = "circos")
    )
  )

  body <- dashboardBody(
    tabItems(
      tabItem("diablo",
              fixedRow(
                column(width = 9,
                       box(width = NULL, solidHeader = TRUE,
                           plotOutput("diablo1", height = 800)
                           ,

                           conditionalPanel(condition = "input.compareDiablo == true",
                                            plotOutput("diablo2", height = 800))

                       )
                ),
                column(width = 3,
                       box(width = NULL, status = "warning",
                           checkboxInput("compareDiablo", "Compare"),
                           p(class = "text-muted",
                             "Compare displays and contrasts another model."
                           )
                       )
                )
              )
      ),

      tabItem("components",
              fluidRow(
                column(width = 9,
                       box(width = NULL, solidHeader = TRUE,
                           plotlyOutput("compplot1", height= "auto")
                           ,
                           conditionalPanel(condition = "input.compareIndiv == true",
                                            plotOutput("indiv2")
                           )

                       ),
                       box(width = NULL,
                           verbatimTextOutput("hoverComp"),
                           verbatimTextOutput("clickComp"),
                           verbatimTextOutput("brushComp"))
                ),
                column(width = 3,
                       box(width = NULL, status = "warning",
                           checkboxInput("compareIndiv", "Compare"),
                           p(class = "text-muted",
                             "Compare displays and contrasts another model."
                           )),
                       box(width = NULL, status = "warning",
                           selectInput("selectDataComp", label = h3("Select data"),
                                       choices = dataNames,
                                       selected = 1),
                           br(),
                           radioButtons("compXComp", label = h3("X Component"),
                                        choices = as.list(1:nComp),
                                        selected = 1),
                           radioButtons("compYComp", label = h3("Y Component"),
                                        choices = as.list(1:nComp),
                                        selected = 2),
                           br(),
                           checkboxInput("showIndNames", label = "Show Ind. Names", value = FALSE),
                           p(class = "text-muted",
                             "Show individual names"
                           )
                       )
                )
              )
      ),

      tabItem("variables",
              fluidRow(
                column(width = 9,
                       box(width = NULL, solidHeader = TRUE,
                           plotlyOutput("var1", height = 800)
                           ,

                           conditionalPanel(condition = "input.compareVar == true",
                                            plotOutput("var2", height = 800))

                       ),
                       box(width = NULL,
                           DT::dataTableOutput("varTable")
                       )
                ),
                column(width = 3,
                       box(width = NULL, status = "warning",
                           checkboxInput("compareVar", "Compare"),
                           p(class = "text-muted",
                             "Compare displays and contrasts another model."
                           )),
                       box(width = NULL, status = "warning",
                           selectInput("selectDataVar", label = h3("Select data"),
                                       choices = dataNames,
                                       selected = 1),
                           br(),
                           selectInput("selectComp", label = h3("Select component"),
                                       choices = as.list(1:nComp),
                                       selected = 1),
                           p(
                             class = "text-muted",
                             paste("Select data to display in table."
                             )
                           ),
                           br(),
                           checkboxInput("showVarNames", label = "Show Var. Names", value = FALSE),
                           p(class = "text-muted",
                             "Show variable names"
                           )
                       )
                )
              )
      ),

      tabItem("biplot",
              fluidRow(
                column(width = 9,
                       box(width = NULL, solidHeader = TRUE,
                           plotlyOutput("biplot1", height = 800)
                           ,
                           conditionalPanel(condition = "input.compareIndiv == true",
                                            plotOutput("biplot2", height = 800))

                       )
                ),
                column(width = 3,
                       # box(width = NULL, status = "warning",
                       #     checkboxInput("compareIndiv", "Compare"),
                       #     p(class = "text-muted",
                       #       "Compare displays and contrasts another model."
                       #     )),
                       box(width = NULL, status = "warning",
                           selectInput("selectDataBi", label = h3("Select data"),
                                       choices = dataNames,
                                       selected = 1),
                           br(),
                           radioButtons("compXBi", label = h3("X Component"),
                                        choices = as.list(1:nComp),
                                        selected = 1),
                           radioButtons("compYBi", label = h3("Y Component"),
                                        choices = as.list(1:nComp),
                                        selected = 2)
                           # br(),
                           # checkboxInput("showIndNames", label = "Show Ind. Names", value = FALSE),
                           # p(class = "text-muted",
                           #   "Show individual names"
                           # )
                       )
                )
              )
      ),

      tabItem("loadingVectors",
              fluidRow(
                column(width = 9,
                       box(width = NULL, solidHeader = TRUE,
                           plotOutput("loadings1", height = 800)
                           ,

                           conditionalPanel(condition = "input.compareLoading == true",
                                            plotOutput("loadings2", height = 800))

                       )
                ),
                column(width = 3,
                       box(width = NULL, status = "warning",
                           checkboxInput("compareLoading", "Compare"),
                           p(class = "text-muted",
                             "Compare displays and contrasts another model."
                           )
                       )
                )
              )
      ),

      tabItem("heatmap",
              fluidRow(
                column(width = 9,
                       box(width = NULL, solidHeader = TRUE,
                           plotOutput("heatmap1", height = 800)
                           ,

                           conditionalPanel(condition = "input.compareHeat == true",
                                            plotOutput("heatmap2", height = 800))

                       )
                ),
                column(width = 3,
                       box(width = NULL, status = "warning",
                           checkboxInput("compareHeat", "Compare"),
                           p(class = "text-muted",
                             "Compare displays and contrasts another model."
                           )
                       )
                )
              )
      ),

      tabItem("network",
              fluidRow(
                column(width = 9,
                       box(width = NULL, solidHeader = TRUE,
                           plotlyOutput("network", height = 800),
                           dataTableOutput("nodes_data_from_shiny")
                       ),
                       box(width = NULL,
                           DT::dataTableOutput("brushNet"),
                           verbatimTextOutput("brushNodes")
                           # verbatimTextOutput("clickNet")
                       ),
                       box(width = NULL,
                           # valueBoxOutput("hoverNetM", width = 6),
                           # valueBoxOutput("clickNetM", width = 6),
                           valueBoxOutput("density"),
                           valueBoxOutput("transitivity"),
                           valueBoxOutput("modularity")
                       )
                ),
                column(width = 3,
                       # box(width = NULL, status = "warning",
                       #     checkboxInput("compare", "Compare"),
                       #     p(class = "text-muted",
                       #       "Compare displays and contrasts another model."
                       #     )
                       # ),
                       box(width = NULL, status = "warning",
                           sliderInput("threshold", label = h3("Cutoff"), min = 0.4,
                                       max = 1, value = 0.7),
                           # plotOutput("histoSlider", height = 100),
                           p(
                             class = "text-muted",
                             paste("Note: Cutoff removes any edges with weight less than the indicated value."
                             )
                           )
                       ),
                       box(width = NULL,
                           fluidRow(valueBoxOutput("hoverNetM", width = 12)),
                           fluidRow(valueBoxOutput("clickNetM", width = 12)))
                )
              )
      ),

      tabItem("circos",
              fluidRow(
                column(width = 9,
                       box(width = NULL, solidHeader = TRUE,
                           plotOutput("circos1", height = 800),

                           conditionalPanel(condition = "input.compareCircos == true",
                                            plotOutput("circos2", height = 800))
                       )
                ),
                column(width = 3,
                       box(width = NULL, status = "warning",
                           checkboxInput("compareCircos", "Compare"),
                           p(class = "text-muted",
                             "Compare displays and contrasts another model."
                           )
                       ),
                       box(width = NULL, status = "warning",
                           sliderInput("cutoff", label = h3("Cutoff"), min = 0.5,
                                       max = 1, value = 0.90),
                           p(
                             class = "text-muted",
                             paste("Note: Cutoff removes any edges with weight less than the indicated value."
                             )
                           )
                       )
                )
              )

      )
    )
  )

  ui <- dashboardPage(
    header,
    sidebar,
    body
  )


  # Server ----
  server <- function(input,output, session) {

    # Compplot ----
    output$compplot1 <- renderPlotly({
      p <- get(input$selectDataComp, ggmixOmics::ggcompplot(M, comps = c(as.numeric(input$compXComp),
                                                                         as.numeric(input$compYComp))))

      ggplotly(p) %>%
        layout(dragmode = "lasso", height = 800, width = 800)
    }
    #, height = function() {
    #  session$clientData$output_compplot1_width
    #}
    )

    output$hoverComp <- renderPrint({
      d <- event_data("plotly_hover")
      if (is.null(d)) "Hover events appear here (unhover to clear)" else d
    })

    output$clickComp <- renderPrint({
      d <- event_data("plotly_click")
      if (is.null(d)) "Click events appear here (double-click to clear)" else d
    })

    output$brushComp <- renderPrint({
      d <- event_data("plotly_selected")
      if (is.null(d)) "Click and drag events (i.e., select/lasso) appear here (double-click to clear)" else d
    })


    # Varplot ----
    plotVar <- mixOmics::plotVar(model1)

    # Get rownames for each Block
    # legacy
    # flowCytometry <- row.names(subset(plotVar, Block == "Flow Cytometry"))
    # luminexCytokine <- row.names(subset(plotVar, Block == "Luminex Cytokine"))
    # metabolomics <- row.names(subset(plotVar, Block == "Metabolomics"))
    # proteomics <- row.names(subset(plotVar, Block == "Proteomics"))
    # transcriptomics <- row.names(subset(plotVar, Block == "Transcriptomics"))

    output$var1 <- renderPlotly({
      p <- get(input$selectDataVar, ggmixOmics::ggvarplot(M))

      ggplotly(p) %>%
        layout(dragmode = "lasso")
    })

    output$varTable <- DT::renderDataTable({
      compN <- paste(c("comp ", input$selectComp), sep="", collapse="")
      table <- as.matrix(model1$loadings[[1]][,compN])
      for(i in 2:nEntries){
        table <- rbind(table, as.matrix(model1$loadings[[i]][,compN]))
      }
      if(input$compare == TRUE){
        tempTable <- as.matrix(model2$loadings$`Flow cytometry`[,compN])
        tempTable <- rbind(tempTable, as.matrix(model2$loadings$'Luminex cytokine'[,compN]))
        tempTable <- rbind(tempTable, as.matrix(model2$loadings$'Metabolomics'[,compN]))
        tempTable <- rbind(tempTable, as.matrix(model2$loadings$'Proteomics'[,compN]))
        tempTable <- rbind(tempTable, as.matrix(model2$loadings$'Transcriptomics'[,compN]))
        table <- cbind(as.matrix(table), as.matrix(tempTable))
      }
      table

      if(input$compare == FALSE)
        colnames(table) <- c('Eigenvector')
      if(input$compare == TRUE)
        colnames(table) <- c('Eigenvector 1', 'Eigenvector 2')

      # table <- as.matrix(table[apply(table[,-1], 1, function(x) !all(x==0)),])
      DT::datatable(data = table,
                    options = list(scrollX = TRUE, scrollY = "275px", autoWidth = FALSE)
      )
    })

    # Biplot ----
    # Make names (TODO)
    samples <- rownames(M$X[[1]])
    output$biplot1 <- renderPlotly({
      p <- get(input$selectDataBi, ggmixOmics::ggbiplot(M, comps = c(as.numeric(input$compXBi),
                                                                     as.numeric(input$compYBi))))
      p$data <- cbind(p$data, samples)

      ggplotly(p) %>%
        layout(dragmode = "lasso")
    })

    # Loadings ----
    output$loadings1 <- renderPlot({mixOmics::plotLoadings(model1, title = "Plot of Loading vectors 1")$graph},
                                   bg = "transparent")
    output$loadings2 <- renderPlot({mixOmics::plotLoadings(model2, title = "Plot of Loading vectors 2")$graph},
                                   bg = "transparent")

    # Circos ----
    output$circos1 <- renderPlot({mixOmics::circosPlot(model1, cutoff = input$cutoff)}, bg = "transparent")
    output$circos2 <- renderPlot({mixOmics::circosPlot(model2, cutoff = input$cutoff)}, bg = "transparent")


    # Heatmap ----
    ncomp <- 2
    corr <- getCorMat(M, method = 'pearson')
    d <- hclust(dist(corr))
    m <- corr[d$order, d$order]
    g1 <- GGally::ggcorr(data = NULL, cor_matrix = m, hjust = 0, size = 0.1, colour = 'white', layout.exp = 0)

    output$heatmap1 <- renderPlot({plot(g1)})

    # Diablo ----
    output$diablo1 <- renderPlot({mixOmics::plotDiablo(model1)$graph}, bg = "transparent")
    output$diablo2 <- renderPlot({mixOmics::plotDiablo(model2)$graph}, bg = "transparent")

    # Network ----
    output$network <- renderPlotly({
      #minimal example
      corThreshold <- input$threshold

      corMat <- getCorMat(M)
      corMat[abs(corMat) < corThreshold] <- 0
      diag(corMat) <- 0

      rownames(corMat) <- make.names(rownames(corMat), unique=TRUE)
      colnames(corMat) <- make.names(colnames(corMat), unique=TRUE)

      graph <- graph.adjacency(abs(corMat), weighted = TRUE, mode = "lower")
      graph <- simplify(graph)

      # Assume ordered
      keeps <- 1:unique(M$ncomp) %>%
        purrr::map(~ mixOmics::selectVar(M, comp = .)) %>%
        purrr::at_depth(2, ~ .x[[1]]) %>%
        purrr::transpose() %>%
        purrr::map(purrr::reduce, union) %>%
        purrr::map(length) %>%
        head(-1)

      V(graph)$group <- rep(names(keeps), keeps)

      # graph information
      output$density <- renderValueBox({
        valueBox(
          value = round(edge_density(graph, loops = FALSE), digits = 3),
          subtitle = "Edge Density",
          icon = icon("anchor")
        )
      })

      output$transitivity <- renderValueBox({
        valueBox(
          value = round(transitivity(graph, type = "global", vids = NULL,
                                     weights = NULL, isolates = c("NaN", "zero")), digits = 3),
          subtitle = "Transitivity",
          icon = icon("wifi")
        )
      })

      output$modularity <- renderValueBox({
        valueBox(
          value = round(modularity(graph, membership(cluster_edge_betweenness(graph))), digits = 3),
          subtitle = "Modularity",
          icon = icon("gavel")
        )
      })

      # convert plot to ggnetwork
      nodesNedges <- ggnetwork(graph)
      nodesNedges$xend <- as.numeric(nodesNedges$xend)
      nodesNedges$yend <- as.numeric(nodesNedges$yend)
      nodesNedges$y <- as.numeric(nodesNedges$y)
      nodesNedges$x <- as.numeric(nodesNedges$x)

      nodes <- nodesNedges[is.na(nodesNedges$weight), c("x", "y", "group" ,"vertex.names")]
      edges <- nodesNedges[!is.na(nodesNedges$weight), c("x", "y", "vertex.names", "xend", "yend")]

      linkNode <- merge(edges, nodes, by.x = c("xend", "yend"), by.y = c("x","y"))

      nodesNedges <- merge(nodesNedges, linkNode, all.x = TRUE)
      nodesNedges <- nodesNedges[order(nodesNedges[,"na.y"], nodesNedges[, "x"], na.last = FALSE),]

      # reorder
      nodesNedges <- nodesNedges[,c("x", "y", "group", "na.x", "vertex.names", "xend", "yend", "na.y", "weight", "vertex.names.x", "vertex.names.y")]


      # PPI Integration ----
      if(PPIIntegration == TRUE){
        # import PPI data
        data <- PPIList

        # match edges and PPI data
        matches <- matchPPI(nodesNedges[, c("vertex.names.x", "vertex.names.y")], data)
        nodesNedges <- cbind(nodesNedges, matches)
      }


      # plotly interactivity
      toPrint <- NULL
      for(i in 1:nEntries){
        toPrint[[length(toPrint)+1]] <- nodesNedges[nodesNedges$group == dataNames[i],]
      }

      output$hoverNetM <- renderValueBox({
        n <- event_data("plotly_hover")
        if (is.null(n)){
          value <- "N/A"
        }
        else
          value <- as.character(toPrint[[n$curveNumber]][n$pointNumber + 1,]$vertex.names)
        valueBox(
          value = value,
          subtitle = "Hovered Node"
        )
      })

      output$clickNetM <- renderValueBox({
        n <- event_data("plotly_click")
        if (is.null(n)){
          value <- "N/A"
        }
        else
          value <- as.character(toPrint[[n$curveNumber]][n$pointNumber + 1,]$vertex.names)
        valueBox(
          value = value,
          subtitle = "Clicked Node"
        )
      })



      output$hoverNet <- renderPrint({
        n <- event_data("plotly_hover")
        if (is.null(n)) "Hover events appear here (unhover to clear)"
        else
          as.character(toPrint[[n$curveNumber]][n$pointNumber + 1,]$vertex.names)
      })

      output$clickNet <- renderPrint({
        n <- event_data("plotly_click")
        if (is.null(n)) "Click events appear here (double-click to clear)"
        else
          as.character(toPrint[[n$curveNumber]][n$pointNumber + 1,]$vertex.names)
      })

      output$brushNet <- DT::renderDataTable({
        n <- event_data("plotly_selected")
        if (is.null(n)){
          # "Click and drag events (i.e., select/lasso) appear here (double-click to clear)"
        }
        else {
          n <- n[c("curveNumber", "pointNumber")]
          n[,2] <- n[,2] + 1

          print <- nodes %>%
            dplyr::tbl_df() %>%
            dplyr::mutate(group = factor(group, levels = dataNames),
                          curveNumber = as.numeric(group)) %>%
            dplyr::group_by(group) %>%
            dplyr::mutate(pointNumber = 1:n()) %>%
            as.data.frame() %>%
            dplyr::inner_join(., n)

          # Geneset Enrichment ----
          if(geneEnrichment == TRUE){
            p <- genesetEnrichment(as.vector(print$vertex.names), featureMapping)
          }
          else{
            p <- as.vector(print$vertex.names)
          }
          DT::datatable(data = p, class = 'compact stripe', colnames = c('Rank' = 1, 'Collection' = 2, 'Geneset' = 3, 'FDR' = 4),
                        caption = htmltools::tags$caption(
                          style = 'caption-side: top; text-align: center;',
                          'Table 1: ', htmltools::em('Geneset Enrichment Results')),
                        options = list(searching = FALSE, autoWidth = TRUE, scrollY = 400, paging = FALSE))
        }
      })

      output$brushNodes <- renderPrint({
        n <- event_data("plotly_selected")
        if (is.null(n)) "Click and drag events (i.e., select/lasso) appear here (double-click to clear)"
        else {
          n <- n[c("curveNumber", "pointNumber")]
          n[,2] <- n[,2] + 1

          print <- nodes %>%
            dplyr::tbl_df() %>%
            dplyr::mutate(group = factor(group, levels = dataNames),
                          curveNumber = as.numeric(group)) %>%
            dplyr::group_by(group) %>%
            dplyr::mutate(pointNumber = 1:n()) %>%
            as.data.frame() %>%
            dplyr::inner_join(., n)

          # Geneset Enrichment ----
          p <- as.vector(print$vertex.names)

          p
        }
      })

      # plot graph
      # library dependency bug
      library(ggnetwork)

      # for (i in 1:nrow(nodesNedges)){
      #   if(!is.na(nodesNedges[i, "weight"])){
      #     nodesNedges[i, "vertex.names"] <- NA
      #   }
      # }

      plot <- ggplot2::ggplot(nodesNedges, ggplot2::aes(x = x, y= y, xend = xend, yend = yend, text = vertex.names)) +
        ggnetwork::geom_edges(size = 0.1, color = "grey50") +
        ggnetwork::geom_nodes(ggplot2::aes(fill = group), size = 6, shape = 21, color = 'white') +
        viridis::scale_fill_viridis('', discrete = TRUE) +
        ggplot2::theme_void()

      ggplot <- ggplotly(plot, tooltip = "text") %>%
        layout(dragmode = "lasso")

      ggplot$x$data[[1]]$hoverinfo <- "none"

      ggplot
    })
  }

  shinyApp(ui = ui, server = server)
}
