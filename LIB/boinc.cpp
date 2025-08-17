#include "boinc.h"

void Initialize(){
#ifdef ENABLE_BOINC
  int rc;
  BOINC_OPTIONS options;


  boinc_init_diagnostics(BOINC_DIAG_REDIRECTSTDERR|
                         BOINC_DIAG_MEMORYLEAKCHECKENABLED|
                         BOINC_DIAG_DUMPCALLSTACKENABLED| 
                         BOINC_DIAG_TRACETOSTDERR);

    
  boinc_options_defaults(options);
  options.main_program = true;
  options.check_heartbeat = true;
  options.direct_process_action = true;

  //rc = boinc_init();

  rc = boinc_init_options(&options);
  
  if(rc)
	  Terminate(5);

#endif
}

void Terminate(int status){
#ifdef ENABLE_BOINC
  boinc_finish(status);
#endif
  printf("Exit code: %d\n",status);
  exit(status);
}

//Input files are uploaded on the client machine
//Thus, filename has to be resolved for each file

int OpenFile_B(char* filename, const char* mode, FILE **f){
  char  mode_[3];
  int   read_only;
  char  filename_old[MAX_PATH__];
  char  resolved_name[MAX_PATH__];

  resolved_name[0] = '\0';
  
#ifdef _WIN32
  strcpy(mode_,mode);
#else
  strcpy(mode_,mode);
#endif

  read_only=0;
  if (!strcmp(mode_,"r")) 
    read_only=1;
  
#ifdef _WIN32
  strcpy(filename_old,filename); 
#else
  strcpy(filename_old,filename); 
#endif

  //Resolve Filename
#ifdef ENABLE_BOINC
  if (!read_only) { //file is write-able
    boinc_resolve_filename(BOINC_OUTFILE,resolved_name,sizeof(resolved_name));  
    //always append BOINC output
    strcpy(mode_,"a");
  }else{
    boinc_resolve_filename(filename,resolved_name,sizeof(resolved_name));
  }

# ifdef _WIN32
  strcpy(filename,resolved_name);
# else
  strcpy(filename,resolved_name);
# endif

#endif
  
  // Open File
#ifdef ENABLE_BOINC
  *f = boinc_fopen(filename,mode_);
#else

# ifdef _WIN32
  *f = fopen(filename,mode_);
# else
  *f = fopen(filename,mode_);
# endif

#endif
  if (*f == NULL){
    fprintf(stderr, "ERROR: Could not open file %s.\n", filename);
    return 0;
  }
  
  // Write File Header
#ifdef ENABLE_BOINC
  if (!read_only)
    fprintf(*f, "-----TEXT FOR FILE %s-----\n",filename_old);
#endif

  return 1;
}

void CloseFile_B(FILE **f, const char* mode){
#ifdef ENABLE_BOINC
  if (strcmp(mode,"r"))
    fprintf(*f, "-----END OF TEXT-----\n");
#endif
  
  fclose(*f);
}

