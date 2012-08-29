#pragma once

#include <sys/time.h>

inline double get_wtime()
{
  struct timeval tv;
  gettimeofday(&tv, 0);
  return tv.tv_sec + 1.e-6 * tv.tv_usec;
}

