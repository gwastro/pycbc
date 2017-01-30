/*
  Helper tool for running Cygwin Apps under BOINC 

  This tool is meant to be called with an absolute (Windows) path to the
  project directory as the first/only argument

  It will create a file etc/fstab with a single entry pointing to the project directory
  such that "boinc-resolved" paths will point to the correct locations

  example:

    fstab C:\Program Data\BOINC Data\projects\albert.phys.uwm.edu

    => etc/fstab:
    C:/Program\040Data/BOINC\040Data/projects/albert.phys.uwm.edu /projects/albert.phys.uwm.edu dummy binary 0 0

    boinc-softlinks will resolve to
    ../../projects/albert.phys.uwm.edu/file

  (C) Bernd Machenschalk 2016,2017
*/

#include <stdio.h>
#include <string.h>
#include <limits.h>
#ifdef _WIN32
#include <windows.h>
#elif defined(__APPLE__)
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#else // Linux
#include <sys/stat.h>
#endif

#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

char path[PATH_MAX];

int main(int argc, char*argv[]) {

#if defined(_WIN32) || defined(TEST_WIN32)

  char *c = argv[1];
  char *n, *proj;
  unsigned int i=0;
  FILE* fp;

  // replace all "\" with "/" in argv[1]
  while (*c) {
    if (*c == '\\') {
      *c = '/';
    }
    c++;
  }

  // point proj to second-last "/"
  while (i<2) {
    c--;
    if (*c == '/')
      i++;
  }
  proj = c;

  // copy argv[1] to path, replacing " " with "\040"
  c = argv[1];
  while ((n = strchr(c,' '))) {
    *n = '\0'; // end the string to copy here
    strcat(path, c);
    strcat(path, "\\040");
    *n = ' '; // restore the original value
    c = n + 1;
  }
  strcat(path, c);

#ifdef _WIN32
  CreateDirectory ("etc", NULL);
  fp = fopen("etc\\fstab", "w");
#else
  mkdir("etc",0755);;
  fp = fopen("etc/fstab", "w");
#endif
  if (!fp)
    return(1);
  fprintf(fp, "%s %s dummy binary 0 0\n", path, proj);
  fclose(fp);

#endif // (TEST)_WIN32

  return(0);
}
