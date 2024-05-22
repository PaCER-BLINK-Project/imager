#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../src/apply_calibration.h"

void usage()
{
   printf("test_calsol CAL_SOLUTIONS_FILE.txt\n");
}

int main(int argc,char* argv[])
{
   if( argc<=1 || strncmp(argv[1],"-h",2)==0 ){
      usage();
   }   

   const char* filename = argv[1];
//   parse_cmdline(argc-1,argv+1);
//   print_parameters();


   CCalSols calsols;
   if( calsols.read_calsolutions_text( filename ) <= 0 ){
      printf("ERROR : could not read calibration solutions from file %s\n",filename);
      return -1;
   }
   
   printf("OK : read %d calibration solutions from file %s\n",int(calsols.size()),filename);
   calsols.show();   
}

