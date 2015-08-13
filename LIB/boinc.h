#ifndef BOINC_H
#define BOINC_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//BOINC APIs used by FlexAID
#ifdef ENABLE_BOINC
# ifdef _WIN32
#  include "boinc_win.h"
#  include "util.h"
# endif

# define BOINC_OUTFILE "out.txt"
# include "config.h"
# include "boinc_api.h"
# include "diagnostics.h"     // boinc_init_diagnostics()
# include "filesys.h"         // boinc_fopen(), etc...
# include "str_util.h"	     // for parse_command_line()
#endif

#define MAX_PATH__ 255
#define MAX_NUM_FILES 15

//File Handling Function Declarations
int    OpenFile_B(char* filename, const char* mode, FILE **f);    // Opens a file handle
void   CloseFile_B(FILE **f, const char* mode);                                 // Closes a file handle

//int    GetIndexFromEnum(char *);
void   Initialize();                   // BOINC initialize function
void   Terminate(int status);          // BOINC finish function

#endif // include guard
