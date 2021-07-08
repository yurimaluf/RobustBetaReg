# .onLoad <- function(...){
#   #msg <- paste("Loading", pkgname)
#   cat("Em caso de algum bug por favor reportar \n") 
# }

.onAttach<-function(...){
  cat("Em caso de bug por favor reportar para: yurimaluf@gmail.com \n")
  #rstudioapi::showDialog( "BugReports","Report bugs at:", url="https://github.com/yurimaluf/RobustBetaReg/issues")
}