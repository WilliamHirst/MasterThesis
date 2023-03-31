
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * Configure.h *                                 galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#ifndef Configure_h
#define Configure_h

#include<iostream.h>
#include<math.h>
#include<string.h> //IMOS20020112

class Configure
{
 public:

   char *galdef_directory;
   char *fits_directory;
   char *adjunct_directory;

   int directory_length;

//interface function prototype
   int init();
};

#endif
