.First.lib <- function(libname, pkgname){
  # now use NAMESPACE, without NAMESPACE this loads C++ code:
  ## library.dynam("Ringo", pkgname, libname)

  ## show vignette in windows menu
  if(.Platform$OS.type=="windows" && require(Biobase) && interactive() && .Platform$GUI=="Rgui"){
    addVigs2WinMenu("Ringo")
  }
} # .First.lib
