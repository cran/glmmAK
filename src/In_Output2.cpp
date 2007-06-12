/*** In_Output2.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//    CREATED:  14/09/2006
// 
//              openFile: 14/09/2006
//      open_File_toRead: 18/01/2007
//
//
// PURPOSE: Utilities for input and output
//          * partially taken from openFile.cpp, in_output.cpp and in_output_GS.cpp of the 'bayesSurv' package
//
/* ********************************************************************************* */

#include "In_Output2.h"

namespace In_Output{

// ====================================================================================
// ***** openFile: open file for output                                          *****/
// ====================================================================================
void
openFile(std::ofstream &ofile,  const std::string &path,  const char &write_flag)
{
  try{
    bool err = false;

    std::string errmess;
    if (write_flag == 'n') {
      std::fstream temp(path.c_str(), std::ios::in);
      if (!temp) {
	ofile.open(path.c_str(), std::ios::out);
      } 
      else {
	temp.close();
	err = true;
      }
    }
    else if (write_flag == 'o')
           ofile.open(path.c_str(), std::ios::out | std::ios::trunc);
         else if (write_flag == 'a')
                ofile.open(path.c_str(), std::ios::out | std::ios::app);
              else {
                errmess = std::string("C++ Error: Incorrect flag for writing to ") + path + ". ";
		returnR error(errmess, 99);
                throw error;
              }

    if (!ofile || err) {
      errmess = std::string("C++ Error: Could not open ") + path + " for writing. ";
      returnR error(errmess, 99);
      throw error; 
    }
    return;
  }
  catch(returnR){
    throw;
  }  
}


/***** ======================================================================== *****/
/***** open_File_toRead: Open file for reading and skip first 'skip' rows       *****/
/***** ======================================================================== *****/
void
open_File_toRead(std::ifstream &file,  const std::string &path,  const int &skip)
{
  std::string errmes;

  file.open(path.c_str(), std::ios::in);
  if (!file){
    errmes = std::string("Error: Could not open ") + path;
    throw returnR(errmes, 99);
  }

  char ch;
  for (int i = 0; i < skip; i++){
    if (file.eof()){
      int ihelp = i + 1;
      errmes = std::string("Error: Reached end of file ") + path + " before "
               + char(ihelp) + std::string(" rows were skipped.");
      throw returnR(errmes, 99);
    }
    file.get(ch);        
    while (ch != '\n') file.get(ch);
  }

  return;
}

}  /** end of the namespace In_Output ***/
