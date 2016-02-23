#include <stdio.h>
#include <string.h>
#include <unistd.h>

#define debug 0

main() {
  FILE*fp=NULL;
  long new=0, old=0;
  // wait for the file to appear
  while(!fp) {
    if (debug) fprintf(stderr, "waiting for stderr.txt to appear\n");
    sleep(1);
    fp=fopen("stderr.txt","r");
  }
  while(1) {
    float a, b, c, d;
    int found=0;
    FILE*fw;
    char buf[1024];
    // wait for the file to be extended
    while(new <= old) {
      if (debug) fprintf(stderr, "waiting for stderr.txt to be extended\n");
      sleep(1);
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
      if(fw = fopen("progress.txt", "w")) {
	fprintf(fw,"%f\n", (a-1) / b + (c-1) / (b * d));
	fclose(fw);
      }
    }
  }
  fclose(fp);
}
