# .First.lib <- function(lib, pkgname, where){
#   library.dynam(pkgname, pkgname, lib)
#   initialize()   
# }

.onLoad <- function(lib, pkgname){
  library.dynam(pkgname, pkgname, lib)
  initialize()   
}
