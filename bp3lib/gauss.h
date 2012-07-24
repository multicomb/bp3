#ifndef GAUSS_H
#define GAUSS_H	1

#include <string>
#include <vector>

using std::string;
using std::vector;

class Gauss
{
  public:
    Gauss(const string, unsigned = 15);
    ~Gauss(){;}
    
    vector<double> Getnode()   const { return node  ; }
    vector<double> Getweight() const { return weight; }
    void print() const;

 private:
    string type;
    vector<double> node, weight;
    void gaussq();
    void gaussq2(vector<double> &);
};

#endif
