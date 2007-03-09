.First.lib <- function(libname, pkgname){
  # Load C++ code:
  library.dynam("Ringo", pkgname, libname)
  
  if(.Platform$OS.type=="windows" && require(Biobase) && interactive() && .Platform$GUI=="Rgui"){
    addVigs2WinMenu("Ringo")
  }
  
} # .First.lib
