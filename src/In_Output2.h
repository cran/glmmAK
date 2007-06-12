/*** In_Output2.h ***/

#ifndef _IN_OUTPUT_TWEE_H_
#define _IN_OUTPUT_TWEE_H_

#include <R.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include "AK_Error.h"

namespace In_Output{

void
openFile(std::ofstream &ofile,  const std::string &path,  const char &write_flag);

void
open_File_toRead(std::ifstream &file,  const std::string &path,  const int &skip);

}

#endif
