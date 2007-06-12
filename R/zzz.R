#*** zzz.R ***/
##
##     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
##              akom@email.cz
##
#* ********************************************************************************* */

.First.lib <- function(lib, pkg)
{
#   require(coda)
   library.dynam("glmmAK", pkg, lib)

   invisible()
}

