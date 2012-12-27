#*** zzz.R ***/
##
##     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
##              arnost.komarek[AT]mff.cuni.cz
##
#* ********************************************************************************* */

.onAttach <- function(libname, pkgname)
#.First.lib <- function(libname, pkgname)
{
   ###library.dynam("glmmAK", pkgname, libname)   ## no more needed, load is provided by useDynLib in NAMESPACE

   packageStartupMessage(paste(
       "\n",
       "### Generalized Linear Mixed Models\n",
       "### Arnost Komarek\n\n",
       "### See citation(\"glmmAK\") or toBibtex(citation(\"glmmAK\")) for the best way to cite\n",
       "### the package if you find it useful.\n\n", sep=""))
   #cat("\n")
   #cat("### Generalized Linear Mixed Models\n")
   #cat("### Arnost Komarek\n\n")
   #cat("### See citation(\"glmmAK\") or toBibtex(citation(\"glmmAK\")) for the best way to cite\n")
   #cat("### the package if you find it useful.\n\n")
   
   invisible()
}


