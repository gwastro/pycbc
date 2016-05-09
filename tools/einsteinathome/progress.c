#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

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
  while (argc >= arg) {
    // -v: verbose
    if (!strcmp("-v",argv[arg])) {
      debug++;
    }
    // -f: 'filter' input file, write non-matching lines to stdout
    if (!strcmp("-f",argv[arg])) {
	filter = 1;
    }
    // -s <seconds>: set sleep period
    if (!strcmp("-s",argv[arg])) {
      arg++;
      delay = atoi(argv[arg]);
    }
    // -p <fraction>: initial progress value
    if (!strcmp("-p",argv[arg])) {
      arg++;
      progress = atof(argv[arg]);
    }
    arg++;
  }
  if (argc >= arg) {
    fname = argv[arg];
  }

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

  while(progress < 1.0) {
    float a, b, c, d;
    int found=0;
    char buf[1024];

    // wait for the file to be extended
    while(new <= old) {
      if (debug) fprintf(stderr, "waiting for stderr.txt to be extended\n");
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
    while (1) {
      if(fgets(buf, sizeof(buf), fp)) {
	if (4 == sscanf(buf, "%*s %*s Filtering template %f/%f segment %f/%f", &a, &b, &c, &d)) {
	  if (debug) fprintf(stderr, "parsed values from line: %f %f %f %f\n", a, b, c, d);
	  found = 1;
	} else if (filter) {
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
  }
  fclose(fp);

  return 0;
}
