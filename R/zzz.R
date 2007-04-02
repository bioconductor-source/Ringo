# now use NAMESPACE, without NAMESPACE this loads C++ code:
#  .First.lib <- function(libname, pkgname)
##    library.dynam("Ringo", pkgname, libname)

.onLoad <- function(libname, pkgname) {
  ## nothing to do here for the moment
}

.onAttach <- function(libname, pkgname) {
  ## show vignette in windows menu
  if(.Platform$OS.type=="windows" && require(Biobase) && interactive() && .Platform$GUI=="Rgui"){
    addVigs2WinMenu("Ringo")
  }
} # .onAttach

.onUnload <- function( libpath ) {
  library.dynam.unload( "Ringo", libpath )
} # .onUnload
