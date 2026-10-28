#-------------------------------------------------------------------------------
# Startup message shown when the package is attached (library(MultiMAP)).
# Flags that MultiMAP has been renamed to dbbMultiArm and is no longer maintained.
#-------------------------------------------------------------------------------
.onAttach <- function(libname, pkgname){
        packageStartupMessage(
                "MultiMAP has been renamed to 'dbbMultiArm' and is no longer ",
                "maintained. Please install and use dbbMultiArm instead:\n",
                "  remotes::install_github(\"MaxMenssen/dbbMultiArm\")")
}
