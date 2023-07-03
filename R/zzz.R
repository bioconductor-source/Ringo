# now use NAMESPACE, without NAMESPACE this loads C++ code:
#  .First.lib <- function(libname, pkgname)
##    library.dynam("Ringo", pkgname, libname)

.onLoad <- function(libname, pkgname) {
  ## nothing to do here for the moment
}

.onAttach <- function(libname, pkgname) {
    msg <- sprintf(
        "Package '%s' is deprecated and will be removed from Bioconductor
         version %s", pkgname, "3.19")
    .Deprecated(msg=paste(strwrap(msg, exdent=2), collapse="\n"))
  ## show vignette in windows menu
  if(.Platform$OS.type=="windows" && interactive() && .Platform$GUI=="Rgui"){
    addVigs2WinMenu("Ringo")
  }
} # .onAttach

.onUnload <- function( libpath ) {
  library.dynam.unload( "Ringo", libpath )
} # .onUnload
