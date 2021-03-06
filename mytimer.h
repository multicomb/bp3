#pragma once

#include <sys/time.h>

#ifndef __TIMER__H__
#define __TIMER__H__

inline static double get_wtime(){
  struct timeval tv;
  gettimeofday(&tv, 0);
  return tv.tv_sec + 1.e-6 * tv.tv_usec;
}

inline static double TIMERcvt(const unsigned long long val)
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

// #include <unordered_map>
#include <map>
#include <string>
#include <cstdio>
#include <cassert>

template <unsigned int N, unsigned int I>
struct FnvHash
{
  inline static unsigned int Hash(const char (&str)[N])
  {
    return (FnvHash<N, I-1>::Hash(str) ^ str[I-1])*16777619u;
  }
};

template <unsigned int N>
struct FnvHash<N, 1>
{
  inline static unsigned int Hash(const char (&str)[N])
  {
    return (2166136261u ^ str[0])*16777619u;
  }
};

class StringHash
{
  unsigned int m_hash;
  public:
    template <unsigned int N>
      inline StringHash(const char (&str)[N]) : m_hash(FnvHash<N, N>::Hash(str)) {}
    unsigned int operator()() const { return m_hash; }
    operator unsigned int() const {return m_hash;}
};

struct TimerT
{
  private:

    enum {NMAX = 100};   /* maximal number of counters */
    typedef unsigned long long Ttype;

    struct Int0
    {
      int val;
      operator int() const {return val; }
      Int0(const int _val = -1) : val(_val) {}
    };

    Ttype initial;
    Ttype elapsed;
    const char* name;
    const char* nameCnt[NMAX];

//    typedef std::unordered_map<unsigned int, Int0> Map;
    typedef std::map<unsigned int, Int0> Map;

    int nCnt;
    Map hash2idx;
    Ttype stamps_beg[NMAX], stamps_end[NMAX], stamps_lapse[NMAX];

    template<const bool FAILCHECK, const int N>
      int string2idx(const char (&tag)[N])
      {
        const unsigned int hash = StringHash(tag);
        Int0 &idx = hash2idx[hash];
        if (idx == -1)
        {
          assert(!FAILCHECK);
          assert(nCnt < NMAX);
          nameCnt[nCnt] = tag;
          idx = nCnt++;
        }
        return idx;
      }

  public:
    TimerT(const char* _name) : name(_name), elapsed(0), nCnt(0)
  {
    initial = TIMER();
//    hash2idx.reserve(10*NMAX);
  }

    template<const int N>
      int start(const char (&tag)[N])
      {
        const int idx = string2idx<false>(tag);
        return start(idx);
      }
    int start(const int idx)
    {
      stamps_beg[idx] = TIMER();
      return idx;
    }
    int stop(const int idx)
    {
      assert(idx >= 0 && idx < nCnt);
      stamps_end  [idx]  = TIMER();
      stamps_lapse[idx] += stamps_end[idx] - stamps_beg[idx];
    }

    void dump() 
    {
      elapsed = TIMER() - initial; 
      fprintf(stderr, "\n");
      fprintf(stderr, " Timer: %s  : elapsed= %g\n", name, TIMERcvt(elapsed));
      fprintf(stderr, "-------------------------\n");
      for (int i = 0; i < nCnt; i++)
        fprintf(stderr, "  > %s: %g \n", nameCnt[i], TIMERcvt(stamps_lapse[i]));
      fprintf(stderr, "=========================\n");
    }
};

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
    const std::string& record(const std::string &tag_old, const std::string &tag_new)
    {
      stop(tag_old);
      return start(tag_new);
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
      //      fprintf(stderr, " > Sum= %g \n", cvt(sum));
      fprintf(stderr, "=========================\n");
    }
};

#endif // __TIMER__H__
