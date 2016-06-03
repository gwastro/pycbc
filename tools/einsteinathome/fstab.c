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
  char *c = argv[1];
  char *n, *proj;
  unsigned int i=0;
  FILE* fp;

#if defined(_WIN32) || defined(TEST_WIN32)

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

#else // (TEST)_WIN32

  int ret;
  char buf[1024];

  n = strrchr(argv[1],'/') + 1;
#define LINK "/tmp/BOINC.project."
  i = strlen(LINK) + strlen(n) + 1;
  proj = (char*)malloc(i * sizeof(char));
  strcpy(proj, LINK);
  strcat(proj, n);
  ret = symlink(argv[1], proj);
  if ((ret != 0) && (errno != EEXIST)) {
    fprintf(stderr, "ERROR: Can't link %s %s: %d\n", argv[1], proj, errno);
    exit(1);
  }
  // replace "../../projects<proj>" by "<LINK>" in all files named on the
  // command-line after the project directory (argv[1])
  // such that boinc_resolves to these paths without blanks
  for (i=2; i<argc; i++) {
    fp = fopen(argv[i], "r");
    if (!fp) {
      fprintf(stderr, "ERROR: Can't open file %s for reading\n", argv[i]);
      exit(1);
    }
    if (!fgets(buf, sizeof(buf), fp) && ferror(fp)) {
      fprintf(stderr, "ERROR: Can't read from file %s\n", argv[i]);
      exit(1);
    }
    fclose(fp);
#define XML "<soft_link>../../projects/"
    if (strncmp(buf, XML, strlen(XML))) {
      fprintf(stderr, "WARNING: file %s is not a XML soft-link to project dir\n", argv[i]);
      fclose(fp);
      continue;
    }
    fp = fopen(argv[i], "w");
    if (!fp) {
      fprintf(stderr, "ERROR: Can't open file %s for writing\n", argv[i]);
      exit(1);
    }
    if (fprintf(fp, "<soft_link>" LINK "%s", buf+strlen(XML)) <= 0) {
      fprintf(stderr, "ERROR: Can't write to file %s\n", argv[i]);
      exit(1);
    }
    fclose(fp);
  }

#endif // (TEST)_WIN32

  return(0);
}
