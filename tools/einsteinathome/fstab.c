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
  char* n = argv[1];
  for(int i=0; i<strlen(argv[1]); i++)
    if (argv[1][i] == '\\')
      argv[1][i] = '/';
  while ((n = strchr(c,' '))) {
    *n = '\0';
    strcat(path, c);
    strcat(path, "\\040");
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
  fclose(fp);
  return(0);
}
