
#' @useDynLib JANE, .registration = TRUE
NULL

JANE_startup_message <- function(){

  msg <- c(paste0(path <- r"(
+-------------------------------------------+
|    ___  ________  ________   _______      |
|   |\  \|\   __  \|\   ___  \|\  ___ \     |
|   \ \  \ \  \|\  \ \  \\ \  \ \   __/|    |
| __ \ \  \ \   __  \ \  \\ \  \ \  \_|/__  |
||\  \\_\  \ \  \ \  \ \  \\ \  \ \  \_|\ \ |
|\ \________\ \__\ \__\ \__\\ \__\ \_______\|
| \|________|\|__|\|__|\|__| \|__|\|_______||
+-------------------------------------------+  version )", 
  utils::packageVersion("JANE")),
  "\nType 'citation(\"JANE\")' for citing this R package in publications.")
  return(msg)
}

.onAttach <- function(lib, pkg)
{
  # startup message
  msg <- JANE_startup_message()
  packageStartupMessage(msg)      
  invisible()
}