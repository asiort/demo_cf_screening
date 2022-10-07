# Demo_CF_screening

### In order to run the app in you local R. 

```
## Install required packages
list.of.packages = c("shiny", "ComplexUpset", "ggplot2")
new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## Run the app
shiny::runGitHub(repo = "demo_CF_screening", username = "asiort", ref = "main")
```


<img src="https://i.pinimg.com/originals/ee/4a/45/ee4a45e1886c42549fcdcc67f4372651.gif" width="350" height="350">
