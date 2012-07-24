#include <math.h>
#include "misc.h"
#include "gauss.h"

Gauss::Gauss(const string typein, unsigned nodes)
{
  type      = typein;

  if (nodes > 1)
  {
    node.resize(nodes);
    weight.resize(nodes, ZERO);
    gaussq();
  }
  else
    Bp3Error("Gauss::Gauss", "more than 1 node is needed");
}

void Gauss::gaussq()
{
  // Computes weights and nodes for gaussian-type quadratures
  // routine modified version from one obtained from GAMS.

  double abi;
  vector<double> scratch(node.size());
 
  double muzero     = ZERO;
  unsigned nm1      = node.size() - 1;
  if (type         == "Legendre")
  {
    muzero          = TWO;
    for (unsigned i = 0; i < nm1; i++)
    {
      node[i]       = ZERO;
      abi           = (double) i + ONE;
      scratch[i]    = abi/sqrt(FOUR*abi*abi - ONE);
    }
    node[nm1]       = ZERO;
  }
  else if (type    == "Hermite")
  {
    muzero          = sqrt(PI);
    for (unsigned i = 0; i < nm1; i++)
    {
      node[i]       = ZERO;
      scratch[i]    = sqrt(((double) i + ONE)/TWO);
    }
    node[nm1]       = ZERO;
  }
  else if (type    == "LaGuerre")
  {
    muzero          = ONE;
    for (unsigned i = 0; i < nm1; i++)
    {
      abi           = (double) i + ONE;
      node[i]       = TWO*abi - ONE;
      scratch[i]    = abi;
    }
    node[nm1]       = TWO*((double) node.size()) - ONE;
  }
  else
    Bp3Error("Gaussq", "Gaussian Quadrature type not defined");

  weight[0]         = ONE;
  
  gaussq2(scratch);

  for (unsigned i   = 0; i < weight.size(); i++)
        weight[i]   = muzero * weight[i] * weight[i];

}

void Gauss::gaussq2(vector<double> &scratch)
{

  /* This subroutine is a translation of an algol procedure,
     num. math. 12, 377-383(1968) by martin and wilkinson,
     as modified in num. math. 15, 450(1970) by dubrulle.
     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).
     this is a modified version of the 'eispack' routine imtql2.

     This subroutine finds the eigenvalues and first components of the
     eigenvectors of a symmetric tridiagonal matrix by the implicit ql
     method.                                                               
  */

  int i, k, m, ii;
  double b, c, f, g, p, r, s;

  // BEGIN                               

  // catch n = 1 case                   
  if (node.size()      == 1)
    return;

  int nm1              = node.size() - 1;
  scratch[nm1]         = ZERO;
  int l                = 0;
  int j                = 0;

  while(l             <= nm1)
  {
    for (m             = l; (m < nm1) && (fabs(scratch[m]) > (EPSILON * (fabs(node[m]) + fabs(node[m+1])))) ; m++);
    p                  = node[l];
    if (m             != l)
    {
      if (j           == 30)
	Bp3Error("Gaussq2","No convergence in 30 Interations of gaussq2");
      j++;
      g                = (node[l+1] - p)/(TWO*scratch[l]);
      r                = sqrt(g*g+ONE);
      if (g           >= ZERO)
	g              = node[m] - p + scratch[l]/(g + r);
      else
	g              = node[m] - p + scratch[l]/(g - r);
      s                = ONE;
      c                = ONE;
      p                = ZERO;
      for(i            = m - 1; i >= l; i--)
      {
	f              = s*scratch[i];
	b              = c*scratch[i];
	if (fabs(f)   >= fabs(g))
	{
	  c            = g/f;
	  r            = sqrt(c*c + ONE);
	  scratch[i+1] = f*r;
	  s            = ONE/r;
	  c           *= s;
	}
	else
	{
	  s            = f/g;
	  r            = sqrt(s*s + ONE);
	  scratch[i+1] = g*r;
	  c            = ONE/r;
	  s           *= c;
	}
	g              = node[i+1] - p;
	r              = (node[i] - g)*s + TWO*c*b;
	p              = s*r;
	node[i+1]      = g + p;
        g              = c*r - b;
        f              = weight[i+1];
        weight[i+1]    = s * weight[i] + c*f;
        weight[i]      = c * weight[i] - s*f;
      }

      node[l]         -= p;
      scratch[l]       = g;
      scratch[m]       = ZERO;
    }
    else
    {
      j                = 0;
      l++;
    }
  }
  // Order Eigenvalues and Eigenvectors 
  for (ii              = 1; ii <= nm1; ii++)
  {
    i                  = ii - 1;
    k                  = i;
    p                  = node[i];

    for (j             = ii; j <= nm1; j++)
    {
      if (node[j]      < p)
      {
	k              = j;
	p              = node[j];
      }
    }
    if (k             != i)
    {
      node[k]          = node[i];
      node[i]          = p;
      p                = weight[i];
      weight[i]        = weight[k];
      weight[k]        = p;
    }
  }
}

void Gauss::print() const
{
  // print out nodes and weights

  printf("Gaussian quadrature type: %s\n", type.c_str());
  printf("number   node    weight\n");
  for (unsigned i   = 0; i < weight.size(); i++)
    printf("%3i  %10.5f  %10.5f\n",i+1, node[i], weight[i]);
  printf("\n");
  fflush(stdout);
}
