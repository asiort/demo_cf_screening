## app.R ##
library(shiny)
library(ComplexUpset)
library(ggplot2)

## Load DATA
load("data/mat_genes.RData")
load("data/sets_MDS.RData")
options = names(all_sets)


ui = fluidPage(
  titlePanel("Demo CF screening"),
  
  ## Panel 1
  fluidRow(
    sidebarPanel(width=2, 
       helpText("Plots the intersecion."),
       
       checkboxGroupInput(inputId = "selection_plot", 
                          label = "Plot options",
                          # fill = TRUE,
                          choices = options,
                          selected = options,
                          inline = FALSE),
       
       actionButton("reset_button_plot_all", "Select all", class = "btn-lg btn-success"),
       actionButton("reset_button_plot_null", "Clear all", class = "btn-danger")
    ),
    
    mainPanel(
      plotOutput("result_plot", height = 500)),
    br()
  ),
  
  ## Panel 2
  sidebarLayout(
    sidebarPanel(width=2,
      helpText("CF in MDS."),
      
      checkboxGroupInput(inputId = "selection", 
                  label = "Choose your option(s).",
                  choices = options),
      
      actionButton("reset_button", "Clear all", class = "btn-danger")
    ),
    
    mainPanel(
      br(), br(),
      textOutput("result")
    )
  )
)
## END OF UI

server <- function(input, output) {
  
  ########
  ##TEXT##
  ########
  len_sel = reactive(length(input$selection))
  
  result_sel = reactive(
    if (len_sel() == 0) {
      "Please select at least an option."
      
    } else if (len_sel() == 1) {
        if (length(all_sets[input$selection][[1]]) == 0) {
          paste(input$selection, " do not have significant genes.")
        } else {
          paste("The gene(s) of ", input$selection, " are: ", paste(all_sets[input$selection][[1]], collapse=" ,"), ".")
        }
        
    } else if (len_sel() == 2) {
        if (length(intersect(all_sets[input$selection[1]][[1]], all_sets[input$selection[2]][[1]])) == 0) {
          paste(input$selection [1], " and ", input$selection[2]," do not share significant genes.")
          
        } else {
          paste("The gene(s) are: ", paste(intersect(all_sets[input$selection[1]][[1]], all_sets[input$selection[2]][[1]]), collapse = " ,"), ".")
          }
      
    } else {
        custom_list = list()
        
        for (i in c(1:len_sel())) {
          custom_list = append(custom_list, list(all_sets[input$selection[i]][[1]])) ## Generate a list with all the options
        }
        
        if (length(Reduce(intersect, custom_list))==0) {
          paste(paste(input$selection[c(1:len_sel()-1)], collapse=", "), "and ", input$selection[len_sel()], " do not share significant genes.")
          
        } else {
        paste("The gene(s) are:",paste(Reduce(intersect, custom_list)), collapse=" ,", ".")
          
        }
      }
  )
  
  ## Reset button
  observeEvent(
    input$reset_button, {updateCheckboxGroupInput(getDefaultReactiveDomain(), "selection", choices = options, selected = NULL)}
  )
  
  output$result = renderText({
    result_sel()
    })
  
  ########
  ##PLOT##
  ########
  len_sel_plot = reactive(length(input$selection_plot))

  ## Reset button plot
  observeEvent(
    input$reset_button_plot_all, {updateCheckboxGroupInput(getDefaultReactiveDomain(), "selection_plot", choices = options, selected = options)}
  )

  observeEvent(
    input$reset_button_plot_null, {updateCheckboxGroupInput(getDefaultReactiveDomain(), "selection_plot", choices = options, selected = NULL)}
  )
  
  output$result_plot = renderPlot({
    tryCatch({
      ComplexUpset::upset(
        mat_genes, input$selection_plot,
        min_degree = 2, n_intersections = 50, mode = "intersect", keep_empty_groups=T, 
        height_ratio = 1.2, themes=upset_default_themes(text=element_text(size=18)))
      
      }, error = function(e) {
        ComplexUpset::upset(
          mat_genes, input$selection_plot,
          min_degree = 1, n_intersections = 50, mode = "intersect", keep_empty_groups=T,
          height_ratio = 1.2, themes=upset_default_themes(text=element_text(size=18)))
        
      })
  })
  
}
## END OF SERVER

shinyApp(ui, server)
