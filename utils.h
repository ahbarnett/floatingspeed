#ifndef UTILS_H
#define UTILS_H

#include <sys/time.h>

class CNTime {
 public:
  void start();
  int restart();
  int elapsed();
  double elapsedsec();
 private:
  struct timeval initial;
};

// crappy random doubles in [0,1)
#define rand01() (((double)rand())/RAND_MAX)

#endif
