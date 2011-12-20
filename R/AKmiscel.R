#*** AKmiscel.R ***/
##
##     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
##              arnost.komarek[AT]mff.cuni.cz
##
##           CREATED:  02/11/2006
##
## PURPOSE: Miscellaneous smaller functions
##
#* ********************************************************************************* */

### "The best" division of the plot region according to the number of figures
### version for landscape plots
mfLand <- function(np)
{
  RET <-
    if (np == 1) c(1, 1)
    else if (np <= 2) c(1, 2)
         else if (np <= 4) c(2, 2)
              else if (np <= 6) c(2, 3)
                   else if (np <= 8) c(2, 4)
                        else if (np <= 9) c(3, 3)
                             else if (np <= 12) c(3, 4)
                                  else if (np <= 16) c(4, 4)
                                       else if (np <= 20) c(4, 5)
                                            else stop("np argument is too high")
  return(RET)  
}  

### "The best" division of the plot region according to the number of figures
### version for portrait plots
mfPort <- function(np)
{
  RET <-
    if (np == 1) c(1, 1)
    else if (np <= 2) c(2, 1)
         else if (np <= 4) c(2, 2)
              else if (np <= 6) c(3, 2)
                   else if (np <= 8) c(4, 2)
                        else if (np <= 9) c(3, 3)
                             else if (np <= 12) c(4, 3)
                                  else if (np <= 16) c(4, 4)
                                       else if (np <= 20) c(5, 4)
                                            else stop("np argument is too high")
  return(RET)  
}  

