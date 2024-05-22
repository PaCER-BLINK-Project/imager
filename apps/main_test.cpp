#include <stdio.h>
#include <stdlib.h>

#include <string>
using namespace std;

int main()
{
   string in_basename = "chan_204_20230205T034618_simul_XX";
   char szDateTime[64],szRest[128];
   int channel;
   
   int ret = sscanf(in_basename.c_str(),"chan_%d_%s_%s",&channel,szDateTime,szRest);
   printf("WARNING : could not parse file name %s using format chan_%%s_%%s_%%s -> only %d strings parsed out (%d,%s,%s)\n",in_basename.c_str(),ret,channel,szDateTime,szRest);
}