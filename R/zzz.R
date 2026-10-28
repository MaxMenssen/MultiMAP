#-------------------------------------------------------------------------------
# Startup message shown when the package is attached (library(MultiMAP)).
# Flags that MultiMAP has been renamed to dbbMultiArm and is no longer maintained.
#-------------------------------------------------------------------------------
.onAttach <- function(libname, pkgname){
        packageStartupMessage(
                "Deprecated — renamed. MultiMAP has been renamed to dbbMultiArm ",
                "and is no longer maintained under this name. Please install and ",
                "use dbbMultiArm instead.")
}
