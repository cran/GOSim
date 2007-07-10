# .First.lib <- function(lib, pkgname, where){
#   library.dynam(pkgname, pkgname, lib)
#   initialize()   
# }

.onLoad <- function(lib, pkgname, where){
  library.dynam(pkgname, pkgname, lib)
  initialize()   
}
