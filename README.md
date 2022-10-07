# Demo_CF_screening
Demo-APP_CF_screening

### In order to run the app in you local R. 

´´´
list.of.packages = c("shiny", "ComplexUpset", "ggplot2")
new.packages ?list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

shiny::runGitHub(repo = "demo_CF_screening", username = "asiort", ref = "main")
´´´
