#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#ifndef _WIN32
#include <signal.h>
#endif

// continue waiting, set to 0 by signal handler
int run = 1;

void sighandler(int sig) {
  run = 0;
}

int main(int argc, char*argv[]) {
  int debug = 0;
  int filter = 0;
  int arg = 1;
  char*fname = "stderr.txt";
  FILE*fp=NULL;
  FILE*fw=NULL;
  long new=0, old=0;
  double progress = 0.0001;
  int delay = 1;

  // parse command-line options
  while (arg < argc) {
    // -v: verbose
    if (!strcmp("-v",argv[arg])) {
      debug++;
    }
    // -f: 'filter' input file, write non-matching lines to stdout
    else if (!strcmp("-f",argv[arg])) {
      filter = 1;
    }
    // -c: 'copy' input file, write all lines to stdout
    else if (!strcmp("-c",argv[arg])) {
      filter = 2;
    }
    // -s <seconds>: set sleep period
    else if(!strcmp("-s",argv[arg])) {
      arg++;
      delay = atoi(argv[arg]);
    }
    // -p <fraction>: initial progress value
    else if(!strcmp("-p",argv[arg])) {
      arg++;
      progress = atof(argv[arg]);
    }
    else if(!strcmp("-h",argv[arg]) || !strcmp("--help",argv[arg])) {
      printf("Usage: progress -[vfsph]* [file]\n"
             "  scans 'file' for lines indicating progress indication\n"
             "  and writes 'fraction done' to a file progress.txt\n"
             " - enable debug output (to stderr) with '-v'\n"
             " - use '-f' ('filter') to write scanned lines not indicatiing progress to stdout\n"
             " - use '-c' ('copy') to write all scanned lines to stdout\n"
             " - adjust sleep period with '-s <seconds>' (defaults to 1s)\n"
             " - specify initial 'fraction done' at progream start with '-p <fraction>' (defaults to 0.0001)\n"
             " - pass name of the file to watch on command-line (last argument, defaults to 'stderr.txt')\n"
             " - terminates when progress reaches 1.0\n"
             );
      exit(0);
    }
    else if(argv[arg][0] != '-') {
      fname = argv[arg];
    }
    arg++;
  }

#ifndef _WIN32
  signal(SIGTERM, sighandler);
#endif

  if (debug) fprintf(stderr, "writing initial progress file\n");
  if((fw = fopen("progress.txt", "w"))) {
    fprintf(fw,"%f\n", progress);
    fclose(fw);
  }

  // wait for the file to appear
  while(!fp) {
    if (debug) fprintf(stderr, "waiting for '%s' to appear\n", fname);
    sleep(delay);
    fp=fopen(fname,"r");
  }

  while(run && progress < 1.0) {
    float a, b, c, d;
    int found=0;
    char buf[1024];

    // wait for the file to be extended
    while(run && new <= old) {
      if (debug) fprintf(stderr, "waiting for '%s' to be extended\n", fname);
      sleep(delay);
      fseek(fp, 0, SEEK_END);
      new = ftell(fp);
    }
    fseek(fp, old, SEEK_SET);

    // scan the last line that matches the format
    // "2016-02-04 22:29:03,745 Filtering template 2045/12454 segment 2/15"
    // - scan matching lines
    // - skip non-matching lines
    // - end scanning on eof
    found = 0;
    while (run) {
      if(fgets(buf, sizeof(buf), fp)) {
	if (4 == sscanf(buf, "%*s %*s Filtering template %f/%f segment %f/%f", &a, &b, &c, &d)) {
	  if (debug) fprintf(stderr, "parsed values from line: %f %f %f %f\n", a, b, c, d);
	  found = 1;
	  if (filter>1) {
	    char* c = &buf[strlen(buf)-1];
	    if (*c == '\n') *c = '\0';
	    printf("%s\n",buf);
	  }
	} else if (filter) {
	  char* c = &buf[strlen(buf)-1];
	  if (*c == '\n') *c = '\0';
	  printf("%s\n",buf);
	} else if (debug) {
	  fprintf(stderr, "non matching line: '%s'\n", buf);
	}
	old = ftell(fp);
      } else {
	break;
      }
    }

    // write progress file
    if (found) {
      if (debug) fprintf(stderr, "writing progress file\n");
      progress = (a-1) / b + (c-1) / (b * d);
      if((fw = fopen("progress.txt", "w"))) {
	fprintf(fw,"%f\n", progress);
	fclose(fw);
      }
    }

  } // while(progress < 1.0)

  fclose(fp);

  return 0;
}
