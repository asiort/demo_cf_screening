## app.R ##
library(shinydashboard)
library(shiny)
library(ComplexUpset)
library(ggplot2)
library(vroom)
library(dplyr)
library(DT)
library(patchwork)
library(shinyBS)
library(grDevices)

## Load DATA
load("data/matrix_and_sets.RData") ## All required list of genes and sets
load("data/inputTables.RData") ## All required tables
load("data/customImagesPart1.RData") ## All images of CF
load("data/customImagesPart2.RData") ## All images of CF
load("data/customImagesComplex.RData") ## All images of complexes

all_matrix_set = list(mat_all, mat_all_filtered, mat_all_screening,
                      mat_all_filtered_screening, filter_list, screening_list,
                      mds_score, mat_dgea)

options_mat = c("Complete", "Filtered", "Screening",
                "Filetered_and_Screening", "Filter_Lara_list", "Screening_list",
                "MDS_score", "DEG_complete_table")

options_only_mat = c("Complete", "Filtered", "Screening",
                "Filetered_and_Screening", "mds_score", "mat_dgea")

names(all_matrix_set) = options_mat
options = names(all_together_set)

## Prepare matrix
rownames(mds_score) = tolower(mds_score[,"Gene"])
mds_score = subset(mds_score, select=c(HSC, CMP, MEP, GMP, AVERAGE))
rownames(mat_dgea) = mat_dgea[,"Gene"]
mat_dgea = mat_dgea %>% mutate_if(is.numeric, format, digits=3) ## Change digits number as desired output
mds_score = mds_score %>% mutate_if(is.numeric, format, digits=3) 

## Declarate variables in order to create the customes images
plot_contain_list = tolower(filter_list[1:122,1])
rownames(cf_complex_df) = tolower(cf_complex_df$cf)
list_png_files = c(list_png_files_part1, list_png_files_part2) ## combine two lists

## START OF UI
ui <- dashboardPage(
  dashboardHeader(title = "Demo CF screening"),
  
  ## Sidebar content
  dashboardSidebar(
    sidebarMenu(
      menuItem("Intersection", tabName = "intersection", icon = icon("list-alt")),
      menuItem("Download content", tabName = "download", icon = icon("table"))
    )
  ),
  
  ## Body content
  dashboardBody(
    tabItems(
      #####################
      ##First tab content##
      #####################
      tabItem(tabName = "intersection",
              ## Intersection plot
              fluidRow(
                box(plotOutput("result_plot", height = 500,),
                    status = "primary",
                    width = 7
                ),
                
                column(width = 4,
                  box(
                    width = NULL,
                    radioButtons("filter", "Lara et al., 2022 TF & CF ", 
                                 c("TRUE", "FALSE"), selected = "TRUE"),
                    br(),
                    radioButtons("filter_screening", "CRISPR screening genes", 
                                 c("TRUE", "FALSE"), selected = "FALSE"),
                    status = "success"
                  ),
                  
                  box(
                    width = NULL,
                    checkboxGroupInput(inputId = "selection_plot", 
                                       label = "Intersection plot options",
                                       # fill = TRUE,
                                       choices = options,
                                       selected = options[-(length(options))], # Except last option
                                       inline = FALSE),
                    
                    actionButton("reset_button_plot_all", "Select all", class = "btn-lg btn-success"),
                    actionButton("reset_button_plot_null", "Clear all", class = "btn-danger"),
                    status = "primary"
                  )
                )
              ),
              
              ## The list of genes
              fluidRow(
                box(
                  title = "Gene names for intersection",
                  textOutput("result_text"),
                  status = "warning",
                  width = 7
                ),
                box(
                  checkboxGroupInput(inputId = "selection_text", 
                                     label = "Gene list.",
                                     # fill = TRUE,
                                     choices = options,
                                     selected = NULL,
                                     inline = FALSE),
                  
                  actionButton("reset_button_text", "Clear all", class = "btn-danger"),
                  downloadButton("download_butt_geneList", "Download .tsv"),
                  actionButton("generate_image_butt" ,"CF Image", icon("paper-plane"),
                               style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                  actionButton("generate_image_complex_butt" ,"Complex Image", icon("paper-plane"),
                               style="color: #fff; background-color: #538e62; border-color: #36543d"),
                  status = "warning",
                  width = 4
                )
              ),
              ## DGEA table
              fluidRow(
                box(
                   tabsetPanel(
                     id = "dgea_table_out",
                     tabPanel("AML", dataTableOutput("dgea_table_outAML")),
                     tabPanel("SCORE_MDS", dataTableOutput("dgea_table_outSCORE")),
                     tabPanel("HSC_MDS", dataTableOutput("dgea_table_outHSC")),
                     tabPanel("CMP_MDS", dataTableOutput("dgea_table_outCMP")),
                     tabPanel("MEP_MDS", dataTableOutput("dgea_table_outMEP")),
                     tabPanel("GMP_MDS", dataTableOutput("dgea_table_outGMP"))
                   ),
                   status  ="warning",
                   width = 11
                )
              ),
              mainPanel( ## Pop up cf
                bsModal(id="custom_image_modal", title="Lara et al., 2022 Supp Fig 1f (Adapted for selected genes)", 
                        trigger="generate_image_butt", size = "large", 
                        plotOutput("image_genes")
                  ),
                
              ),
              mainPanel( ## Pop up complex
                bsModal(id="custom_image_complex_modal", title="Lara et al., 2022 Fig 1d (Adapted for selected complexes)", 
                        trigger="generate_image_complex_butt", size = "xl", 
                        plotOutput("image_complex")
                ),
                
              )
      ),
      
      # Second tab content
      tabItem(tabName = "download",
              fluidRow(
                box(
                  selectInput("download_mat", "Select a matrix", options_mat),
                  dataTableOutput("preview_mat"),
                  downloadButton("download_butt", "Download .tsv"),
                  width = 11
                )
              )
      )
    )
  )
)
## END OF UI

## START OF SERVER
server <- function(input, output) {
#-------------------------------------------------------------------------------
# Tab intersection
  ########################
  ##FILTER SELECTION SET##
  ########################
  
  selected_set = reactive(
    if (input$filter == "TRUE" && input$filter_screening == "TRUE") {
      all_together_filtered_screening_set 
    } else if (input$filter == "TRUE" && input$filter_screening == "FALSE") {
      all_together_filtered_set
    } else if (input$filter == "FALSE" && input$filter_screening == "TRUE") {
      all_together_screening_set
    } else {
      all_together_set
    }
  )
  
  selected_mat = reactive(
    if (input$filter == "TRUE" && input$filter_screening == "TRUE") {
      mat_all_filtered_screening 
    } else if (input$filter == "TRUE" && input$filter_screening == "FALSE") {
      mat_all_filtered
    } else if (input$filter == "FALSE" && input$filter_screening == "TRUE") {
      mat_all_screening
    } else {
      mat_all
    }
  )  
  
  ###############
  ##TEXT OUTPUT##
  ###############
  len_sel = reactive(length(input$selection_text)) ## Number of selected options
  
  result_sel = reactive(
    if (len_sel() == 0) {
      "Please select at least an option."
      
    } else if (len_sel() == 1) {
      if (length(selected_set()[input$selection_text][[1]]) == 0) {
        paste(input$selection_text, " do not have significant genes.")
      } else {
        paste("The gene(s) of ", input$selection_text, " are: ", paste(selected_set()[input$selection_text][[1]], collapse=", "), ".")
      }
      
    } else if (len_sel() == 2) {
      if (length(intersect(selected_set()[input$selection_text[1]][[1]], selected_set()[input$selection_text[2]][[1]])) == 0) {
        paste(input$selection_text[1], " and ", input$selection_text[2]," do not share significant genes.")
        
      } else {
        paste("The gene(s) are: ", paste(intersect(selected_set()[input$selection_text[1]][[1]], selected_set()[input$selection_text[2]][[1]]), collapse = ", "), ".")
      }
      
    } else {
      custom_list = list()
      
      for (i in c(1:len_sel())) {
        custom_list = append(custom_list, list(selected_set()[input$selection_text[i]][[1]])) ## Generate a list with all the options
      }
      
      if (length(Reduce(intersect, custom_list))==0) {
        paste(paste(input$selection_text[c(1:len_sel()-1)], collapse=", "), "and ", input$selection_text[len_sel()], " do not share significant genes.")
        
      } else {
        paste("The gene(s) are:", paste(Reduce(intersect, custom_list), collapse = ", "), collapse=", ", ".")
        
      }
    }
  )  

  
  result_sel_gene = reactive(
    if (len_sel() == 0) {
      NULL
    } else if (len_sel() == 1) {
      selected_set()[input$selection_text][[1]]
    } else {
      custom_list_intersect = list()
      for (i in c(1:len_sel())) { ## START for loop
        custom_list_intersect = append(custom_list_intersect, list(selected_set()[input$selection_text[i]][[1]])) 
      } ## END for loop
      Reduce(intersect, custom_list_intersect)
    }
  )
  
  
  result_plot = reactive(
    tryCatch({
      ComplexUpset::upset(
        selected_mat(), input$selection_plot,
        min_degree = 2, n_intersections = 50, mode = "intersect", keep_empty_groups=T, 
        height_ratio = 1.2, themes=upset_default_themes(text=element_text(size=18)))
      
    }, error = function(e) { ## In case there is no intersection
      if (conditionMessage(e) == "Needs at least two indicator variables" || 
          conditionMessage(e) == "'x' must be an array of at least two dimensions") {
        
      } else {
      ComplexUpset::upset(
        selected_mat(), input$selection_plot,
        min_degree = 1, n_intersections = 50, mode = "intersect", keep_empty_groups=T,
        height_ratio = 1.2, themes=upset_default_themes(text=element_text(size=18)))
      }
    })
  )
  
  
  ##################
  ##GENERATE IMAGE##
  ##################
  cf_paper_png = reactive(
    if (length(result_sel_gene()[which(result_sel_gene() %in% plot_contain_list)]) == 0) {
      NULL
      } else {
      result_sel_gene()[which(result_sel_gene() %in% plot_contain_list)]
      }
  )

  ccomplex_paper_png = reactive(
    if (length(result_sel_gene()[which(result_sel_gene() %in% rownames(cf_complex_df))]) == 0) {
      NULL
    } else {
      result_sel_gene()[which(result_sel_gene() %in% rownames(cf_complex_df))]
    }
  )

  obtain_complex = function(cf_paper_names) {
    if (is.null(cf_paper_names)) {
      return(NULL)
    } else {
      return(unique(cf_complex_df[cf_paper_names,"complex"]))
    }
  }
  
  
  create_image = function(cf_paper_names) {
    if (is.null(cf_paper_names)) {
      return(NULL)
    } else {
    return(patchwork::wrap_plots(list_png_files[c("header", cf_paper_names)], ncol = 1))
    }
  }

  create_image_complex = function(complex_names) {
    if (is.null(complex_names)) {
      return(NULL)
    } else {
      
      return(patchwork::wrap_plots(list_png_files_complex[c("header_complex", complex_names)], ncol = 1))
    }
  }
  

  ###########
  ##BUTTONS##
  ###########
  
  ## Reset button text
  observeEvent(
    input$reset_button_text, {updateCheckboxGroupInput(getDefaultReactiveDomain(), "selection_text", choices = options, selected = NULL)}
  )
  ## Select all button plot
  observeEvent(
    input$reset_button_plot_all, {updateCheckboxGroupInput(getDefaultReactiveDomain(), "selection_plot", choices = options, selected = options)}
  )
  ## Reset button plot
  observeEvent(
    input$reset_button_plot_null, {updateCheckboxGroupInput(getDefaultReactiveDomain(), "selection_plot", choices = options, selected = NULL)}
  )

  ###########################
  ##DOWNLOAD INTERSECT LIST##
  ###########################
  
  download_geneList = reactive(
    mat_dgea[result_sel_gene(),]
    )
  
  output$download_butt_geneList = downloadHandler(
    filename = function() {
      name_list = paste(input$selection_text[c(1:len_sel())], collapse="__")
      paste0(name_list, "_intersection_geneList.tsv")
    },
    content = function(file) {
      vroom::vroom_write(download_geneList(), file)
    }
  )
  
  
  #################
  ##RENDER OUTPUT##
  #################
  output$result_text = renderText({
    result_sel()
  })  
  
  output$result_plot = renderPlot({
    result_plot()
  })
  
  output$dgea_table_outAML = DT::renderDataTable({
    tryCatch({
    DT::datatable(subset(mat_dgea[result_sel_gene(),], 
                         select = c(CD34vsBlasto_AML_logFC, CD34vsBlasto_AML_padj,
                                    LCSvsBLASTO_AML_logFC, LCSvsBLASTO_AML_padj))) %>% formatStyle(
        columns = c("CD34vsBlasto_AML_logFC", "LCSvsBLASTO_AML_logFC"), 
        valueColumns = c("CD34vsBlasto_AML_logFC", "LCSvsBLASTO_AML_logFC"), 
        backgroundColor = 
          styleInterval(seq(-2, 2, by=0.25), 
                        unique(c(colorRampPalette(c("#1164ff", "white"))(length(seq(-2, 0, by=0.25))),
                                 colorRampPalette(c("white", "#f96244"))(length(seq(0, 2, by=0.25))+1))
                        ))
      )
    }, error=function(cond) { ## In case no subset is selected not display table
     return(NULL) 
    })
  })
  
  output$dgea_table_outHSC = DT::renderDataTable({
    tryCatch({
      DT::datatable(subset(mat_dgea[result_sel_gene(),], 
                           select = c(HSC_EB_MDS_logFC, HSC_EB_MDS_padj,
                                      HSC_MLD_MDS_logFC, HSC_MLD_MDS_padj,
                                      HSC_RS_MDS_logFC, HSC_RS_MDS_padj,
                                      HSC_SLD_MDS_logFC, HSC_SLD_MDS_padj))) %>% formatStyle(
        columns = c("HSC_EB_MDS_logFC", "HSC_MLD_MDS_logFC", 
                    "HSC_RS_MDS_logFC", "HSC_SLD_MDS_logFC"), 
        valueColumns = c("HSC_EB_MDS_logFC", "HSC_MLD_MDS_logFC", 
                         "HSC_RS_MDS_logFC", "HSC_SLD_MDS_logFC"), 
        backgroundColor = 
          styleInterval(seq(-2, 2, by=0.25), 
                        unique(c(colorRampPalette(c("#1164ff", "white"))(length(seq(-2, 0, by=0.25))),
                                 colorRampPalette(c("white", "#f96244"))(length(seq(0, 2, by=0.25))+1))
                        ))
      )
    }, error=function(cond) {
      return(NULL) 
    })
  })

  output$dgea_table_outCMP = DT::renderDataTable({
    tryCatch({
      DT::datatable(subset(mat_dgea[result_sel_gene(),], 
                           select = c(CMP_EB_MDS_logFC, CMP_EB_MDS_padj,
                                      CMP_MLD_MDS_logFC, CMP_MLD_MDS_padj,
                                      CMP_RS_MDS_logFC, CMP_RS_MDS_padj,
                                      CMP_SLD_MDS_logFC, CMP_SLD_MDS_padj))) %>% formatStyle(
        columns = c("CMP_EB_MDS_logFC", "CMP_MLD_MDS_logFC", 
                    "CMP_RS_MDS_logFC", "CMP_SLD_MDS_logFC"), 
        valueColumns = c("CMP_EB_MDS_logFC", "CMP_MLD_MDS_logFC", 
                         "CMP_RS_MDS_logFC", "CMP_SLD_MDS_logFC"), 
        backgroundColor = 
          styleInterval(seq(-2, 2, by=0.25), 
                        unique(c(colorRampPalette(c("#1164ff", "white"))(length(seq(-2, 0, by=0.25))),
                                 colorRampPalette(c("white", "#f96244"))(length(seq(0, 2, by=0.25))+1))
                        ))
      )
    }, error=function(cond) { 
      return(NULL) 
    })
  })
  
  output$dgea_table_outMEP = DT::renderDataTable({
    tryCatch({
      DT::datatable(subset(mat_dgea[result_sel_gene(),], 
                           select = c(MEP_EB_MDS_logFC, MEP_EB_MDS_padj,
                                      MEP_MLD_MDS_logFC, MEP_MLD_MDS_padj,
                                      MEP_RS_MDS_logFC, MEP_RS_MDS_padj,
                                      MEP_SLD_MDS_logFC, MEP_SLD_MDS_padj))) %>% formatStyle(
        columns = c("MEP_EB_MDS_logFC", "MEP_MLD_MDS_logFC", 
                    "MEP_RS_MDS_logFC", "MEP_SLD_MDS_logFC"), 
        valueColumns = c("MEP_EB_MDS_logFC", "MEP_MLD_MDS_logFC", 
                         "MEP_RS_MDS_logFC", "MEP_SLD_MDS_logFC"), 
        backgroundColor = 
          styleInterval(seq(-2, 2, by=0.25), 
                        unique(c(colorRampPalette(c("#1164ff", "white"))(length(seq(-2, 0, by=0.25))),
                                 colorRampPalette(c("white", "#f96244"))(length(seq(0, 2, by=0.25))+1))
                        ))
      )
    }, error=function(cond) {
      return(NULL) 
    })
  })
  
  output$dgea_table_outGMP = DT::renderDataTable({
    tryCatch({
      DT::datatable(subset(mat_dgea[result_sel_gene(),], 
                           select = c(GMP_EB_MDS_logFC, GMP_EB_MDS_padj,
                                      GMP_MLD_MDS_logFC, GMP_MLD_MDS_padj,
                                      GMP_RS_MDS_logFC, GMP_RS_MDS_padj,
                                      GMP_SLD_MDS_logFC, GMP_SLD_MDS_padj))) %>% formatStyle(
        columns = c("GMP_EB_MDS_logFC", "GMP_MLD_MDS_logFC", 
                    "GMP_RS_MDS_logFC", "GMP_SLD_MDS_logFC"), 
        valueColumns = c("GMP_EB_MDS_logFC", "GMP_MLD_MDS_logFC", 
                         "GMP_RS_MDS_logFC", "GMP_SLD_MDS_logFC"), 
        backgroundColor = 
          styleInterval(seq(-2, 2, by=0.25), 
                        unique(c(colorRampPalette(c("#1164ff", "white"))(length(seq(-2, 0, by=0.25))),
                          colorRampPalette(c("white", "#f96244"))(length(seq(0, 2, by=0.25))+1))
                          ))
      )
    }, error=function(cond) {
      return(NULL)
    })
  })
  
  output$dgea_table_outSCORE = DT::renderDataTable({
    tryCatch({
      DT::datatable(mds_score[result_sel_gene(),]) %>% formatStyle( 
        columns = colnames(mds_score), 
        valueColumns = colnames(mds_score), 
        backgroundColor = 
          styleInterval(seq(0, 1, by=0.1), 
            colorRampPalette(c("white", "#f96244"))(length(seq(0, 1, by=0.1))+1)) 
        )
    }, error=function(cond) {
      return(NULL)
    })
  })
  
  output$image_genes = renderPlot({
    create_image(cf_paper_png())
  })
    
  output$image_complex = renderPlot({
    create_image_complex(obtain_complex(ccomplex_paper_png()))
  })
  

#-------------------------------------------------------------------------------
# Tab download
  
  ######################
  ##DOWNLOAD REFERENCE##
  ######################
  
  download_mat_sel = reactive(input$download_mat)
  download_df = reactive(
    if (input$download_mat %in% options_only_mat) {
    cbind(Genes=rownames(all_matrix_set[[input$download_mat]]), 
          all_matrix_set[[input$download_mat]])
    } else {
      all_matrix_set[[input$download_mat]]
    }
  )
  
  output$preview_mat = renderDataTable({
    download_df()
  })
  
  output$download_butt = downloadHandler(
    filename = function() {
      paste0(input$download_mat, ".tsv")
    },
    content = function(file) {
      vroom::vroom_write(download_df(), file)
    }
  )
}
## END OF SERVER

shinyApp(ui, server)
## END APP
