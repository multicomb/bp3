#pragma once

#include <sys/time.h>

#ifndef __TIMER__H__
#define __TIMER__H__

inline double get_wtime(){
  struct timeval tv;
  gettimeofday(&tv, 0);
  return tv.tv_sec + 1.e-6 * tv.tv_usec;
}

typedef unsigned long long TIMERout;
inline double TIMERcvt(const TIMERout val)
{
  return (double)val*1.0e-9;
}
inline static unsigned long long TIMER()
{
  unsigned int low,high;
  __asm__ __volatile__("rdtsc" : "=a" (low), "=d" (high));
  unsigned long long count = (unsigned long long)low + (((unsigned long long) high)<<32);
  return count;
}

#include <map>
#include <string>
#include <cstdio>
#include <cassert>



struct Timer
{
  private:

    struct TimeVar
    {
      unsigned long long val;

      unsigned long long operator()() const {return val; }
      TimeVar(const unsigned long long _val = 0) : val(_val) {}
      unsigned long long operator+=(const unsigned long long x) {val += x;}
    };

    typedef std::map<std::string, TimeVar> Map;

    unsigned long long initial;
    unsigned long long elapsed;
    std::string name;
    Map stamps_beg;
    Map stamps_end;
    Map stamps_lapse;

    inline static unsigned long long TIMER()
    {
      unsigned int low,high;
      __asm__ __volatile__("rdtsc" : "=a" (low), "=d" (high));
      unsigned long long count = (unsigned long long)low + (((unsigned long long) high)<<32);
      return count;
    }
    inline static double cvt(const unsigned long long val) { return (double)val * 1.0e-9; }

  public:
    Timer(const std::string &_name) : name(_name), elapsed(0) { initial = TIMER(); }
    void start() { initial = TIMER(); elapsed = 0;}
    void clear() { stamps_beg.clear(); stamps_end.clear(); stamps_lapse.clear();}
    unsigned long long stop () 
    { 
      if (elapsed == 0)
        elapsed = TIMER() - initial; return elapsed; 
    }
    const std::string& start(const std::string &tag) 
    {
      const unsigned long long t = TIMER(); 
      stamps_beg[tag] = t;
      return tag;
    }
    void stop(const std::string &tag) 
    { 
      const unsigned long long t = TIMER(); 
      const Map::const_iterator &it = stamps_beg.find(tag);
      assert(it != stamps_beg.end());
      stamps_end[tag] = t;
      stamps_lapse[tag] += t - it->second();
    }

    void dump() 
    {
      fprintf(stderr, "\n");
      unsigned long long sum = 0;
      fprintf(stderr, " Timer: %s  : %selapsed= %g\n", name.c_str(), elapsed == 0 ? "*" :" ", cvt(stop()));
      fprintf(stderr, "-------------------------\n");
      for (Map::const_iterator it = stamps_lapse.begin(); it != stamps_lapse.end(); it++)
      {
        fprintf(stderr, "  > %s: %g \n", it->first.c_str(), cvt(it->second()));
        sum += it->second();
      }
      fprintf(stderr, " > Sum= %g \n", cvt(sum));
      fprintf(stderr, "=========================\n");
    }
};

#endif // __TIMER__H__
