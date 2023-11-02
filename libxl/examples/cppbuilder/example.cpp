#pragma hdrstop
#pragma argsused

#ifdef _WIN32
#include <tchar.h>
#else
  typedef char _TCHAR;
  #define _tmain main
#endif

#include <stdio.h>
#include "libxl.h"

using namespace libxl;

int _tmain(int argc, _TCHAR* argv[]) 
{
	Book* book = xlCreateXMLBook();
	Sheet* sheet = book->addSheet("my");
	sheet->writeStr(5, 5, "Hello, World!");
	book->save("out.xlsx");
	book->release();
	return 0;
}
