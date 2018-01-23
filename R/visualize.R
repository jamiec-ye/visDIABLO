#' Visualize DIABLO models in an interactive environment
#' @author Jamie C. Ye <jamiec.ye@@gmail.com>
#
#' @param model A DIABLO model object.
#' @param featureMapping A list of data frames containing 'Data.Names', 'Gene.Symbols', 'Display.Names' for each of the data blocks.
#' @import shiny
#' @import shinythemes
#' @import shinydashboard
#' @import shinyBS
#' @importFrom DT dataTableOutput renderDataTable datatable
#' @importFrom igraph graph.adjacency simplify edge_density transitivity cluster_edge_betweenness
#' @import visNetwork
#' @importFrom plotly renderPlotly ggplotly plotlyOutput
#' @import ggmixOmics
#' @import network
#' @import sna
#' @import ggnetwork
#' @import ggplot2
#' @import tidyverse
#' @export

visualize <- function(model, featureMapping = NULL) {
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
    )
  )

  body <- dashboardBody(
    tabItems(
      tabItem("biplot",
              fluidRow(
                column(width = 9,
                       box(width = NULL, solidHeader = TRUE,
                           plotly::plotlyOutput("biplot1", height = 800)
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

      tabItem("network",
              fluidRow(
                column(width = 9,
                       box(width = NULL, solidHeader = TRUE,
                           plotly::plotlyOutput("network", height = 800),
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

    # Biplot ----
    # Make names (TODO)
    samples <- rownames(M$X[[1]])
    output$biplot1 <- plotly::renderPlotly({
      p <- get(input$selectDataBi, ggmixOmics::ggbiplot(M, comps = c(as.numeric(input$compXBi),
                                                                     as.numeric(input$compYBi))))
      p$data <- cbind(p$data, samples)

      plotly::ggplotly(p) %>%
        layout(dragmode = "lasso")
    })

    # Network ----
    output$network <- plotly::renderPlotly({
      #minimal example
      corThreshold <- input$threshold

      corMat <- getCorMat(M)
      corMat[abs(corMat) < corThreshold] <- 0
      diag(corMat) <- 0

      rownames(corMat) <- make.names(rownames(corMat), unique=TRUE)
      colnames(corMat) <- make.names(colnames(corMat), unique=TRUE)

      graph <- igraph::graph.adjacency(abs(corMat), weighted = TRUE, mode = "lower")
      graph <- igraph::simplify(graph)

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
          value = round(igraph::edge_density(graph, loops = FALSE), digits = 3),
          subtitle = "Edge Density",
          icon = icon("anchor")
        )
      })

      output$transitivity <- renderValueBox({
        valueBox(
          value = round(igraph::transitivity(graph, type = "global", vids = NULL,
                                     weights = NULL, isolates = c("NaN", "zero")), digits = 3),
          subtitle = "Transitivity",
          icon = icon("wifi")
        )
      })

      output$modularity <- renderValueBox({
        valueBox(
          value = round(modularity(graph, membership(igraph::cluster_edge_betweenness(graph))), digits = 3),
          subtitle = "Modularity",
          icon = icon("gavel")
        )
      })

      # convert plot to ggnetwork
      nodesNedges <- ggnetwork::ggnetwork(graph)
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

      ggplot <- plotly::ggplotly(plot, tooltip = "text") %>%
        layout(dragmode = "lasso")

      ggplot$x$data[[1]]$hoverinfo <- "none"

      ggplot
    })
  }

  shinyApp(ui = ui, server = server)
}
