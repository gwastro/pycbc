#include <stdio.h>
#include <string.h>
#include <limits.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <sys/stat.h>
#endif

#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

char path[PATH_MAX];

int main(int argc, char*argv[]) {
  char* c = argv[1];
  char *n, *proj;
  unsigned int i=0;

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
  FILE* fp = fopen("etc\\fstab", "w");
#else
  mkdir("etc",0755);;
  FILE* fp = fopen("etc/fstab", "w");
#endif
  if (!fp)
    return(1);
  fprintf(fp, "%s /project dummy binary 0 0\n", path);
  fprintf(fp, "%s %s dummy binary 0 0\n", path, proj);
  fclose(fp);
  return(0);
}
