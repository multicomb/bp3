#include <complex>
#include <vector>
#include <string>
#include "tabfunc.h"
#include <iomanip>
#include <iostream>


using std::string;
using namespace std;

namespace multivar_llhood {

template <typename realnum>
class Gauss;

template <typename realnum>
class covar_matrix;

#if defined _MSC_VER && !defined atanh
#define atanh(x) 0.5*log((1+x)/(1-x));
#endif

/*
template <typename realnum>
struct realnum_traits { };
struct realnum_traits<float> 
{
  static inline char Type() { return realnum_type; }
  private:
  const static char realnum_type = 'f';
};
struct realnum_traits<double> 
{
  static inline char Type() { return realnum_type; }
  private:
  const static char realnum_type = 'd';
};
*/

template <typename realnum>
struct lapack_workspace
{
  int info;
  char ch1, ch2;
  int dim_w, dim_rw, dim_iw;
  int num1, num2, cmplx;
  realnum* rwork;
  int* iwork;
  complex<realnum>* inout1_c;
  complex<realnum>* inout2_c;
  complex<realnum>* work_c;
  realnum* inout1_r;
  realnum* inout2_r;
  realnum* work_r;

  lapack_workspace() 	{ ch1='V'; ch2='U'; dim_w=0; dim_rw=0; dim_iw=0; num1=0; num2=0; cmplx=-1; inout1_r=0; inout1_c=0; }
  lapack_workspace(lapack_workspace& l_in) 											{ 
	  ch1=l_in.ch1; ch2=l_in.ch2; dim_w=l_in.dim_w; dim_rw=l_in.dim_rw; dim_iw=l_in.dim_iw; 
	  cmplx=l_in.cmplx; num1=l_in.num1; num2=l_in.num2; Alloc(num1,num2,cmplx); 	}
  lapack_workspace<realnum>& operator=(lapack_workspace& l_in) 	{ if (inout1_r||inout1_c) DeAlloc();
	  ch1=l_in.ch1; ch2=l_in.ch2; dim_w=l_in.dim_w; dim_rw=l_in.dim_rw; dim_iw=l_in.dim_iw; 
	  cmplx=l_in.cmplx; num1=l_in.num1; num2=l_in.num2; Alloc(num1,num2,cmplx); return *this;	}
  ~lapack_workspace<realnum>()  {  DeAlloc(); }
  void Alloc(int n1, int n2, int cmplx );
  void DeAlloc();
  inline void alloc_test(void* name, char* func_name="") { if (!name) /* likelihood<realnum>::Error( 101, func_name )*/; }
};



template <typename realnum>
struct workspace
{
  int dim1, dim2, dim3, dim4, dim5, dim6, dim7, dim8, dim9, dim10;	// workspaces dimensions
  realnum *space1;			// Num
  realnum **space2;			// Num x Num
  realnum *space3;			// Gauss_number
  realnum ****space4;		// cov.N_part x Num x Num x Num
  realnum ****space5;		// cov.N_part x Gauss number x Num x Num
  realnum **space6;			// Gauss_number_2 x Gauss_number (for N_meas == 3)
  realnum **space7;			// Gauss number x Num
  realnum ***space8;		// cov.N_part x Num x Num
  realnum *space9;			// Gauss number_2
  realnum **space10;		// Gauss number_2 x Num
};


template <typename realnum>
struct Gauss_structure
{
	int num_points;
	std::vector<realnum> node;
    realnum *weight;
    realnum *sin, *cos;
    realnum *sin2, *cos2;
};


template <typename realnum>
class likelihood
{
  friend class covar_matrix<realnum>;
  friend struct workspace<realnum>;
  friend struct lapack_workspace<realnum>;
  friend struct Gauss_structure<realnum>;
  friend class Gauss<realnum>;
  public:
    int Num;						// total number of measurements + models
    int N_meas;						// number of measurements
  // hide this!!! and write what required...
    int Num_p, N_meas_p;			// stored values for prior phase info (ppi) function (for no ppi Rice 1+1 is always used)  
    int Gauss_number;				// number of Gaussian nodes for integration over the last phase
    int Gauss_number_2;				// number of Gaussian nodes for integration over the last-1 phase (in case of N_meas>=3)
	int cent;						// centricity of reflection
    int Rice;						// flag to force using Rice distribution ( default -1=no force, 0=no Rice(ppi), 1=Rice)
  // hide this and/or solve the duplicity with covmat.h
    int no_imag;					// if no_imag then only real parts will be taken into account (default 0)
	int checkOK;					// flag to check if requested values (of derivatives) are available (default 0=no check)
    covar_matrix<realnum>* cov;		// all for building covariance matrix and the matrix itself
    realnum* (eigenvalues1[2]); 
	realnum* (eigenvalues2[2]);
	realnum* inv_eigenval;
	TabFunc<realnum>* tab;

    likelihood( int N_measur, int Number, int do_der_Fph = 2, int do_der_D = 2, string Dtype = "single",
	   int N_part = 1, int Gauss_Number_Of_Nodes = 20, int Gauss_Number_Of_Nodes_2 = 1, 
	   int no_imagin = 1, int Tab_func = 1, Gauss_structure<realnum>* G_s_copy = 0 );
	likelihood() : cov_ppi(this), cov_rice(this) { a = 0; };
	likelihood( likelihood<realnum>& );
	likelihood& operator=( likelihood<realnum>& );
	~likelihood()   { DeAlloc(); }

	realnum GetFOM()  							{ int v[2] = {0,0};			if ( !checkOK || GetOK(v) ) ;  return FOM;  }
	realnum GetHLA()  							{ int v[2] = {4,0};			if ( !checkOK || GetOK(v) ) ;  return HLA;  }
	realnum GetHLB()  							{ int v[2] = {4,0};			if ( !checkOK || GetOK(v) ) ;  return HLB;  }
	realnum GetHLC()  							{ int v[2] = {4,0};			if ( !checkOK || GetOK(v) ) ;  return HLC;  }
	realnum GetHLD()  							{ int v[2] = {4,0};			if ( !checkOK || GetOK(v) ) ;  return HLD;  }
	realnum GetPHIB()  							{ int v[2] = {0,0};			if ( !checkOK || GetOK(v) ) ;  return PHIB;  }
	realnum GetDer_F( int i )					{ int v[3] = {1,1,i};		if ( !checkOK || GetOK(v) ) ;  return der_F[i];  }
	realnum GetDer_ph( int i )					{ int v[3] = {1,1,i};		if ( !checkOK || GetOK(v) ) ;  return der_ph[i];  }
	realnum GetDer2_FF( int i, int j )			{ int v[4] = {1,2,i,j};		if ( !checkOK || GetOK(v) ) ;  return der2_F2[i][j];  }
	realnum GetDer2_phF( int i, int j )			{ int v[4] = {1,2,i,j};		if ( !checkOK || GetOK(v) ) ;  return der2_phF[i][j];  }
	realnum GetDer2_phph( int i, int j )		{ int v[4] = {1,2,i,j};		if ( !checkOK || GetOK(v) ) ;  return der2_ph2[i][j];  }
	realnum GetDer_D( int p, int i, int j )		{ int v[5] = {3,1,p,i,j};	if ( !checkOK || GetOK(v) ) ;  return der_D[p][i][j];  }
	realnum GetDer2_DD( int p, int k, int l, int p2, int m, int n )	
	  { int v[8] = {3,2,p,k,l,p2,m,n};	if ( !checkOK || GetOK(v) ) ; return der2_D[p][k][l][p2][m][n]; }
	realnum GetMinEig() { return min_acc_eigen; }
	void SetMinEig( realnum new_min_acc ) { min_acc_eigen = new_min_acc; }
    int GetDerD() { return do_der_D; }
    int GetDerF() { return do_der_Fph; }
    int GetCalcHL() { return do_ABCD; }
    string GetDType() { return cov_ppi.Dtype; }		// returns type of prior phase information use
	int GetTabFunc() { return Tab_func; } 
	int GetNoIntegLast() { return no_integ_last; } 
	int GetImproveInteg() { return improve_integ; } 
	int GetMaxPoints() { return max_num_points_allowed; };
	int GetSaferNumPoints() { return safer_num_points; } 
    void InverseAndEigen();
    realnum EvaluateSR( realnum *F, realnum *phase ); 
    realnum EvaluateSR_MLD( realnum *F1, realnum *phase1, realnum *F2, realnum *phase2 ); // temporary, might want to change structure later...
//    void SetDer( int do_Fph, int do_D ) { do_der_Fph = do_Fph; do_der_D = do_D; }
    void SetCalcHL(int y) { do_ABCD=y; }
	void SetCentRice( int centricity, int rice = -1 );
	void SetNoIntegLast ( int nointlast );
	void SetImproveInteg ( int impint ) { if (N_meas<=2) improve_integ = 0; else improve_integ = impint;  }
	void SetMaxPoints ( int max_number ) { max_num_points_allowed = max_number; }
	void SetSaferNumPoints( int num ) { safer_num_points = num; } 
	void SetNum4PHIB ( int index ) { num4phib = index; }
	realnum TransfToProb();		// transforms the llhood values of function and given deriv. to probability function and it's deriv.
	void PrintInvMat( int modelpart_only ); // prints real part of actual inverse matrices (array a or c)
	void PrintEigenv( int modelpart_only ); // prints eigenvalues of actual cov matrices (array a or c)
    
//  private:

	int Tab_func;		// flag to use tabulated versions of exp,I0,I1 (quicker but require bit more memory) (default 1=tab.)
	int no_integ_last;		// if true then the last observed phase won't be integrated out (default is 0); FOM is not defined if true
	int improve_integ;		// if >= 1, the values of (w) coefs in integrand are checked and based on their values, 
							// the integration is optimised, giving higher precision a/o less points required
							// if 1, then the order of integrated variables is optimised and the sampling (of all variables) is automatically optimised (based on the w values)
							// if 2, then the order of integrated variables is optimised and the transformation is done so that there are only a few important values in a certain range of second variable for integration
							// if 3, then 1+2 
							// if 4, then 3 and only the important integ values for the third variable are taken in integration
							// only if N_meas > 2 (SIRAS function, maybe MAD function in the future)
							// default 1 for SIRAS, 0 otherwise (change this later?)
	int improve_integ_now, Gauss_number_now, Gauss_number_now_2;				// used internally inside EvaluateSR
	double max_num_points_allowed;	// maximum allowed number of points - if exceeded, then evalution is aborted and negative number of points is returned
								// only taken into account when improve_integ >= 1 (default 1e100 ie always evaluate)
	int safer_num_points;	// if > 0 then the number of points for each integration variable is increased (ie higher precision, slower evaluation)
							// the (positive) value gives the rate of increase (number of increase steps)
							// default 0, only effective if improve_integ>=1
	realnum fval;								// function value
    realnum *der_F, *der_ph;					// derivatives wrt amplitudes and phases
    realnum **der2_F2, **der2_phF, **der2_ph2;	// 2. derivatives wrt amplitudes and phases
    realnum ***der_D;							// derivatives wrt D
    realnum ******der2_D;						// 2. derivatives wrt D
    realnum *auxder_F, *auxder_ph, **auxder2_F2, **auxder2_phF, **auxder2_ph2, ***auxder_D, ******auxder2_D;	// auxiliary arrays for no_integ_last
    realnum FOM, PHIB;							// figure of merite and best phase; computed if do_der_Fph flag is on
    realnum HLA,HLB,HLC,HLD;					// Hendrickon-Lattman's ABCD coefficients
    realnum **a, **b, **c, **d;					// real and imaginary parts of inverses to covariance matrices cov1, cov2
    int do_der_Fph, do_der_D;					// internal flags - which derivatives are to be computed  
	int do_ABCD;								// calculate HL ABCD coefficients or not
    realnum *Gauss_weight;
    realnum *Gauss_weight_2;				// used if integrating over 2 phases
    realnum *cos_last, *sin_last, *sin_last2, *cos_last2, **cos_last12_del;
	Gauss_structure<realnum> *Gauss_struct;		// stores one or more Gaussian samplings and other precomputed Gauss related values
	int num_Gauss_struct;						// number of Gauss_struct = 1 if improve_integ < 0, more than 1 otherwise
    covar_matrix<realnum> cov_ppi;			// using prior phase information (f.i. for acentrics in SAD)
    covar_matrix<realnum> cov_rice;			// no prior phase information - Rice 1x1 (f.i. for centrics in SAD)
	realnum **w, *phw;	// stores the w coefs
	int maxw_ind[2]; 	// stores the indices of max. w 
	int num_changing;	// number of permutations in the order of variables (only used if improve_integ>0)
	int num4phib;		// the index (<N_meas) of variable that should be used for phib calculation
						// if num4phib<0 then do nothing otherwise force num4phib index to be used for phib calc when switching variables
						// this can slow the evalution so should only be set to nonnegative values when necessary
						// default -1
	realnum sum_w,sum_wm; 		// parameter estimating the height of the peak, used for automatic sampling (only used if improve_integ > 2)
	realnum w_crit, w_crit_2;  
	workspace<realnum> work;
	lapack_workspace<realnum> lpck;
	realnum min_acc_eigen, min_nonzero;			// treating numbers close to zero; set up at likelihood(...)
	realnum max_hyparg;							// max. argument for hyperbolic functions
	int Bess_prec1, Bess_prec2;					// precision of (not tabulated) Bessel functions (in number of skipped terms in i0e&i1e); best available = 0
	int debug_info;
	int nil_m;
    void InverseMatrix( unsigned N, realnum** Inp_Re, realnum** Inp_Im, 
		realnum **Inver_Re, realnum **Inver_Im, realnum *Eigenvalues, realnum* Inv_Eigenval ); 
	void InitIntegImpr( realnum *F, realnum *phase );
	void ChangeIntegOrder( realnum *F, realnum *phase, int back );
	int GetOK( int* );
	void Alloc( Gauss_structure<realnum>* );
	void DeAlloc();
	static void Error( int NumErr, const char* Func, void* Info = 0 );
    inline void alloc_test(void* name, char* func_name="") { if (!name) Error( 101, func_name ); }

};





template<typename realnum>
void likelihood<realnum>::Error( int NumErr, const char* Func, void* Info )
{
  char ErrMess[200];
  switch (NumErr) {
    case 101: strcpy( ErrMess, " Allocation of memory not succesful.(Not enough memory?) " ) ;	 	break; 
    case 102: strcpy( ErrMess, " At least 1 model is needed. " ) ;  								break; 
    case 103: strcpy( ErrMess, " More than 1 node is needed for Gauss. " ) ;  						break;
    case 104: strcpy( ErrMess, " Gaussian Quadrature type not defined " ) ;							break;
    case 105: strcpy( ErrMess, " No convergence in 30 Iterations of gaussq2 " ) ;					break;
    case 106: strcpy( ErrMess, " Don't know how to construct cov matrix of dimension " ) ; 			break;	
	case 107: strcpy( ErrMess, " The derivative flag must one of {0,1,2} . " ) ; 					break;	
	case 108: strcpy( ErrMess, " Unknown type of D parameters arrangement " );						break;
	case 109: strcpy( ErrMess, " Attempt to access the not calculated value of " );					break;
	case 110: strcpy( ErrMess, " Number of partial models must be 1, 2, 3 or 4" );					break;
	case 111: strcpy( ErrMess, " Too many part. models for SIR(AS) or MAD " );						break;
	case 112: strcpy( ErrMess, " Can't open file with tabulated values " );							break;
	case 113: strcpy( ErrMess, " The second derivatives not implemented for given target " );		break;
	case 114: strcpy( ErrMess, " Internal error when trying to set cov terms to zero " );			break;
	case 201: strcpy( ErrMess, " Negative cov. matrix term " );										break;
	case 202: strcpy( ErrMess, " The number of partials was reset - given target requires " );		break;
  }
  if ( NumErr < 200 ) 	// errors
  {
	cerr << "Error " << NumErr << " in function " << Func << " : " << ErrMess ;
	cout << "Error " << NumErr << " in function " << Func << " : " << ErrMess ;
  }
  else 					// warnings
	cout << "Warning " << NumErr << " from function " << Func << " : " << ErrMess ;  
  if ( NumErr == 106 ) { cerr << *(int*) Info ; cout << *(int*) Info ; }
  if ( NumErr == 109 ) 
  { 
	int info = *(int*) Info; 
	if ( info == 0 ) { cerr << "FOM or PHIB";  cout << "FOM or PHIB";} else
	if ( info == 1 ) { cerr << "derivative wrt F or phase"; cout << "derivative wrt F or phase";  } else
	if ( info == 3 ) { cerr << "derivative wrt D";  cout << "derivative wrt D"; }
	if ( info == 4 ) { cerr << "HL coefficients";  cout << "HL coefficients"; }
  }
  // if ( NumErr == 110 ) { cerr << "(" << *(int*) Info << " is not supported). "; cout << "(" << *(int*) Info << " is not supported). "; }
  if ( NumErr == 112 ) { cerr << (char*) Info; cout << (char*) Info; }
  if ( NumErr == 113 ) { cerr << "\"" << *(string*) Info << "\""; cout << "\"" << *(string*) Info << "\""; }
  if ( NumErr == 201 ) { cout << "(" << ((int*)Info)[0] << "," << ((int*)Info)[1] << ") . Set to 0. "; }
  if ( NumErr == 202 ) { cout << *(int*) Info << " partials! "; }
  if ( NumErr < 200 ) 	cerr << endl;
  cout << endl; 
  if ( NumErr < 200 )	exit(1);
}


template<typename realnum>
int likelihood<realnum>::GetOK( int* v )
{
  int OK = 1;
  int dimension = 0;
  if ( v[0] == 0 && !do_der_Fph && !do_der_D ) 	OK = 0;
  else if ( v[0] == 1 ) 
  {	
	dimension = Num;
	if ( do_der_Fph < v[1] ) 	OK = 0;
  }
  else if ( v[0] == 3 ) 
  {	
	dimension = cov->DNum;
	if ( do_der_D < v[1] ) 		OK = 0;
  }
  else if ( v[0] == 4 ) 
	if ( do_ABCD < 1 ) 			OK = 0;
  
  int check = v[0]*v[1] + 2;
  for ( int c = 2;  c < check;  c++ )
  {
	if ( v[0] == 3 ) { if ( ( c == 2 || c == 5 ) && ( v[c] >= cov->N_part || v[c] < 0 ) )  OK = 0; }
	else if ( v[c] >= dimension || v[c] < 0 )		OK = 0;
  }
  if ( !OK ) Error( 109, "likelihood:GetOK", (void*) &v[0] );
  return OK;
}


//#define alloc_test(name)  if (!name) Error( 101, func_name )

template <typename realnum>
void lapack_workspace<realnum>::DeAlloc()	
{
  delete [] iwork;
  if ( inout1_c )
  {
	delete [] rwork;
	delete [] inout1_c;
	delete [] inout2_c;
	delete [] work_c;
  }
  if ( inout1_r )
  {
	delete [] inout1_r;
	delete [] inout2_r;
	delete [] work_r;
  }
}


#define alloc_test_init(var, type, dim, func) \
{ var = new type[dim]; alloc_test(var,func); memset(var,'\0',sizeof(type)*dim); }

template <typename realnum>
void lapack_workspace<realnum>::Alloc( int n1, int n2, int complexxx )
{ 
  char* func_name = (char*) "lapack_workspace::Alloc";  
  num1 = n1; num2 = n2; 
  cmplx = complexxx;
  dim_iw = 3 + 5*num1; 
  alloc_test_init( iwork, int, dim_iw, func_name )
//  iwork = new int[dim_iw];								alloc_test(iwork,func_name);
//	memset(iwork,'\0',sizeof(int)*dim_iw);
  inout1_r = 0;
  inout1_c = 0;
  if ( cmplx )
  {
	dim_w = 2*num1 + num1*num1; 
	dim_rw = 1 + 5*num1 + 2*num1*num1; 
	alloc_test_init( rwork, realnum, dim_rw, func_name )
//	rwork = new realnum[dim_rw];						alloc_test(rwork,func_name);
//	  memset(rwork,'\0',sizeof(realnum)*dim_rw);
	alloc_test_init( inout1_c, complex<realnum>, num1*num1, func_name )
//	inout1_c = new complex<realnum>[num1*num1];			alloc_test(inout1_c,func_name);
//	  memset(inout1_c,'\0',sizeof(realnum)*num1*num1);
	alloc_test_init( inout2_c, complex<realnum>, num2*num2, func_name )
//	inout2_c = new complex<realnum>[num2*num2];			alloc_test(inout2_c,func_name);
//	  memset(inout2_c,'\0',sizeof(realnum)*num2*num2);
	alloc_test_init( work_c, complex<realnum>, dim_w, func_name )
//	work_c = new complex<realnum>[dim_w];				alloc_test(work_c,func_name);
//	  memset(work_c,'\0',sizeof(realnum)*dim_w);
  }
  else
  {
	dim_w = 1 + 6*num1 + 2*num1*num1; 
	alloc_test_init( inout1_r, realnum, num1*num1, func_name )
//	inout1_r = new realnum[num1*num1];					alloc_test(inout1_r,func_name);
//	  memset(inout1_r,'\0',sizeof(realnum)*num1*num1);
	alloc_test_init( inout2_r, realnum, num2*num2, func_name )
//	inout2_r = new realnum[num2*num2];					alloc_test(inout2_r,func_name);
//	  memset(inout2_r,'\0',sizeof(realnum)*num2*num2);
	alloc_test_init( work_r, realnum, dim_w, func_name )
//	work_r = new realnum[dim_w];						alloc_test(work_r,func_name);
//	  memset(work_r,'\0',sizeof(realnum)*dim_w);
  }
}


template<typename realnum>
void likelihood<realnum>::DeAlloc()
{
  char* func_name = (char*) "likelihood::DeAlloc";
  int N_mod = Num_p - N_meas_p;
		    // precomputing of weights and sin(node), cos(node) of Gauss_number of nodes for integration
  if ( N_meas_p >= 2 )
  {
	for ( int i=0; i<num_Gauss_struct; i++ )
	{
  	  delete [] Gauss_struct[i].sin;
  	  delete [] Gauss_struct[i].cos;
	  if ( N_meas_p == 3 ) {
  		delete [] Gauss_struct[i].sin2;
  		delete [] Gauss_struct[i].cos2;
	  }
	  delete [] Gauss_struct[i].weight;
	}
	delete [] Gauss_struct;
	if ( N_meas_p == 3 )
	{
  	  for ( int k=0;  k<Gauss_number_2;  k++ )	delete [] cos_last12_del[k];
  	  delete [] cos_last12_del;
  	  for ( int k=0;  k<work.dim6*Gauss_number_2;  k++ )	delete [] work.space6[k];
  	  delete [] work.space6;
  	  for ( int k=0;  k<work.dim10*Gauss_number_2;  k++ )
  	  {
		delete [] work.space10[k];
	  }
	  delete [] work.space10;
	}
	delete [] work.space3;
	if ( N_meas_p == 3 )
	  delete [] work.space9;
  	for ( int k=0;  k<work.dim7*Gauss_number;  k++ )
  	{
	  delete [] work.space7[k];
	}
	delete [] work.space7;
	if ( do_der_D == 2 )
	{
	  for ( int p = 0;  p < work.dim5*cov_ppi.N_part;  p++ )
	  {
  		for ( int k=0;  k<Gauss_number;  k++ )
  		{
		  for ( int i=0;  i<Num_p;  i++ )
		  {  
			delete [] work.space5[p][k][i];
		  }
		  delete [] work.space5[p][k];	
		}
		delete [] work.space5[p];
	  }
	  delete [] work.space5;
	}
  }

  if ( do_der_Fph )
  {
    delete [] der_F;
    delete [] der_ph;
    delete [] auxder_F;
    delete [] auxder_ph;
  }
  if ( do_der_Fph == 2 )
  {
    for ( int i = 0;  i < Num_p;  i++ ) 
    {
      delete [] der2_F2[i] ;
      delete [] der2_ph2[i];
      delete [] der2_phF[i];
      delete [] auxder2_F2[i] ;
      delete [] auxder2_ph2[i];
      delete [] auxder2_phF[i];
    }
    delete [] auxder2_F2;
    delete [] auxder2_ph2;
    delete [] auxder2_phF;
    delete [] der2_F2 ;
    delete [] der2_ph2;
    delete [] der2_phF;
  }

  if ( do_der_D )
  {
	for ( int p = 0;  p < work.dim8*cov_ppi.N_part;  p++ )
	{
	  for ( int i = 0;  i < Num_p;  i++ ) 
	  {
    	if ( p<cov_ppi.N_part ) delete [] der_D[p][i];
    	if ( p<cov_ppi.N_part ) delete [] auxder_D[p][i];
    	delete [] work.space8[p][i];
	  }
  	  if ( p<cov_ppi.N_part ) delete [] der_D[p];
  	  if ( p<cov_ppi.N_part ) delete [] auxder_D[p];
	  delete [] work.space8[p]; 
	}
    delete [] der_D;
    delete [] auxder_D;
	delete [] work.space8;
  }

  if ( do_der_D == 2 )
  {
	for ( int p = 0;  p < cov_ppi.N_part;  p++ )
	{
	  for ( int i = 0;  i < Num_p;  i++ )
	  {
		for ( int j = 0;  j < Num_p;  j++ )
		{
		  for ( int p2 = 0;  p2 < cov_ppi.N_part;  p2++ ) 
		  {
			for ( int k = 0;  k < Num_p;  k++ ) 
			{
			  delete [] der2_D[p][i][j][p2][k];
			  delete [] auxder2_D[p][i][j][p2][k];
			}
			delete [] der2_D[p][i][j][p2];
			delete [] auxder2_D[p][i][j][p2];
		  }
		  delete [] der2_D[p][i][j];
		  delete [] auxder2_D[p][i][j];
		}
    	delete [] der2_D[p][i];
    	delete [] auxder2_D[p][i];
	  }
  	  delete [] der2_D[p];
  	  delete [] auxder2_D[p];
	}
    delete [] der2_D;
    delete [] auxder2_D;

    for ( int p = 0;  p < work.dim4*cov_ppi.N_part;  p++ ) 
    {
	  for ( int i = 0;  i < Num_p;  i++ ) 
	  {
		for ( int j = 0;  j < Num_p;  j++ ) 
		{
		  delete [] work.space4[p][i][j];
		}
		delete [] work.space4[p][i];
	  }
	  delete [] work.space4[p];
    }
	delete [] work.space4;
  }

  delete [] work.space1;
  for ( int i = 0;  i < work.dim2*Num_p;  i++ ) 
  {
    delete [] work.space2[i];
  }    
  delete [] work.space2;
  for ( int i = 0;  i < Num_p;  i++ )   
  {
    delete [] a[i];    delete [] b[i];
    delete [] c[i];    delete [] d[i];
  }
  delete [] a;    delete [] b;
  delete [] c;    delete [] d;

  delete [] eigenvalues1[0] ;
  delete [] eigenvalues2[0] ;
  delete [] eigenvalues1[1] ;
  delete [] eigenvalues2[1] ;
  delete [] inv_eigenval ;

  if ( Tab_func )
	delete tab;

}



template<typename realnum>
void likelihood<realnum>::Alloc( Gauss_structure<realnum>* G_s_copy )
{
  char* func_name = (char*) "likelihood::Alloc";
  int N_mod = Num_p - N_meas_p;
		    // precomputing of weights and sin(node), cos(node) of Gauss_number of nodes for integration
  if ( N_meas_p >= 2 )
  {
// alocating all the Gauss_struct
	int Gauss_num_points_record[18] = { Gauss_number, Gauss_number_2, 8, 12, 15, 20, 25, 35, 50, 70, 100, 135, 180, 250, 350, 500, 700, 1000 };
	num_Gauss_struct = 1;   
	if ( N_meas_p == 3 ) num_Gauss_struct++;
	int Gauss_number_max = Gauss_number;
	if ( N_meas_p == 3 && Gauss_number_2 > Gauss_number )   Gauss_number_max = Gauss_number_2;
	if ( improve_integ >= 1 ) 
	  for (int i=num_Gauss_struct; Gauss_num_points_record[i]<Gauss_number_max && i<18; i++)  num_Gauss_struct++;
	Gauss_struct = new Gauss_structure<realnum>[num_Gauss_struct];
	for ( int i=0; i<num_Gauss_struct; i++ )
	{
	  Gauss_struct[i].num_points = Gauss_num_points_record[i];
	  alloc_test_init( Gauss_struct[i].sin, realnum, Gauss_num_points_record[i], func_name )
//  	  Gauss_struct[i].sin = new realnum[Gauss_num_points_record[i]]; 	alloc_test(Gauss_struct[i].sin,func_name);
//		memset(Gauss_struct[i].sin,'\0',sizeof(realnum)*Gauss_num_points_record[i]);
	  alloc_test_init( Gauss_struct[i].cos, realnum, Gauss_num_points_record[i], func_name )
//  	  Gauss_struct[i].cos = new realnum[Gauss_num_points_record[i]];   	alloc_test(Gauss_struct[i].cos,func_name);
//		memset(Gauss_struct[i].cos,'\0',sizeof(realnum)*Gauss_num_points_record[i]);
	  if ( N_meas_p == 3 ) {
	  alloc_test_init( Gauss_struct[i].sin2, realnum, Gauss_num_points_record[i], func_name )
//  		Gauss_struct[i].sin2 = new realnum[Gauss_num_points_record[i]]; 	alloc_test(Gauss_struct[i].sin2,func_name);
//		  memset(Gauss_struct[i].sin2,'\0',sizeof(realnum)*Gauss_num_points_record[i]);
	  alloc_test_init( Gauss_struct[i].cos2, realnum, Gauss_num_points_record[i], func_name )
//  		Gauss_struct[i].cos2 = new realnum[Gauss_num_points_record[i]];   	alloc_test(Gauss_struct[i].cos2,func_name);
//		  memset(Gauss_struct[i].cos2,'\0',sizeof(realnum)*Gauss_num_points_record[i]);
	  }
	  alloc_test_init( Gauss_struct[i].weight, realnum, Gauss_num_points_record[i], func_name )
//  	  Gauss_struct[i].weight = new realnum[Gauss_num_points_record[i]]; alloc_test(Gauss_struct[i].weight,func_name);
	  if (!G_s_copy) {
// calculate all the nodes and weights using Gauss class
		Gauss<realnum> gausslegendre("Legendre", Gauss_num_points_record[i]);
  		Gauss_struct[i].node  = gausslegendre.Getnode();    
		vector<realnum> weight = gausslegendre.Getweight();
  		for ( int k=0;  k<Gauss_num_points_record[i];  k++ ) {
    	  realnum m = Gauss_struct[i].node[k]*PI + PI ;	// < 0, 2*PI >
    	  Gauss_struct[i].sin[k] = sin(m); Gauss_struct[i].cos[k] = cos(m);
		  if ( N_meas_p == 3 ) {
    		m = Gauss_struct[i].node[k]*HALF_PI + PI ;  // < 1/2*PI, 3/2*PI >
    		Gauss_struct[i].sin2[k] = sin(m); Gauss_struct[i].cos2[k] = cos(m);
		  }
		  Gauss_struct[i].weight[k] = weight[k];
		}	
  	  } 
  	  else {
// copy all the nodes and weights from passed G_s_copy
  		for ( int k=0;  k<Gauss_num_points_record[i];  k++ ) {
    	  Gauss_struct[i].sin[k] = G_s_copy[i].sin[k]; Gauss_struct[i].cos[k] = G_s_copy[i].cos[k];
		  if ( N_meas_p == 3 ) {
    		Gauss_struct[i].sin2[k] = G_s_copy[i].sin2[k]; Gauss_struct[i].cos2[k] = G_s_copy[i].cos2[k];
		  }
		  Gauss_struct[i].weight[k] = G_s_copy[i].weight[k];
		}
	  }
	}
	work.dim3 = 4;		//if ( (do_der_D && N_meas_p==2) )  work.dim3 = 4;   
	alloc_test_init( work.space3, realnum, work.dim3*Gauss_number, func_name )
//	work.space3 = new realnum[work.dim3*Gauss_number];   			alloc_test(work.space3,func_name);
	if ( N_meas_p == 3 )
	{
	  work.dim9 = 3;
	  alloc_test_init( work.space9, realnum, work.dim9*Gauss_number_2, func_name )
//	  work.space9 = new realnum[work.dim9*Gauss_number_2];   		alloc_test(work.space9,func_name);
	}
	work.dim7 = 2;
	work.space7 = new realnum*[work.dim7*Gauss_number];				alloc_test(work.space7,func_name);
  	for ( int k=0;  k<work.dim7*Gauss_number;  k++ )
  	{
	  alloc_test_init( work.space7[k], realnum, Num_p, func_name )
//	  work.space7[k] = new realnum[Num_p];							alloc_test(work.space7[k],func_name);
	}
  }

  if ( N_meas_p == 3 )
  {
	work.dim10 = 2;
	work.space10 = new realnum*[work.dim7*Gauss_number_2];			alloc_test(work.space10,func_name);
  	for ( int k=0;  k<work.dim10*Gauss_number_2;  k++ )
  	{
	  alloc_test_init( work.space10[k], realnum, Num_p, func_name )
//	  work.space10[k] = new realnum[Num_p];							alloc_test(work.space10[k],func_name);
	}
    cos_last12_del = new realnum*[Gauss_number_2]; 					alloc_test(cos_last12_del,func_name);
  	for ( int k=0;  k<Gauss_number_2;  k++ )
	{
	  alloc_test_init( cos_last12_del[k], realnum, Gauss_number, func_name )
//	  cos_last12_del[k] = new realnum[Gauss_number];				alloc_test(cos_last12_del[k],func_name);
	}
	work.dim6 = 2;			if ( do_der_D )   work.dim6 = 4;
	work.space6 = new realnum*[work.dim6*Gauss_number_2];   		alloc_test(work.space6,func_name);
  	for ( int k=0;  k<work.dim6*Gauss_number_2;  k++ )
  	{
	  alloc_test_init( work.space6[k], realnum, Gauss_number, func_name )
//	  work.space6[k] = new realnum[Gauss_number];					alloc_test(work.space6[k],func_name);
	}
  }

  if ( N_meas_p >= 2 )
  {
	if ( do_der_D == 2 )
	{
	  work.dim5 = 1;
	  work.space5 = new realnum***[work.dim5*cov_ppi.N_part];		alloc_test(work.space5,func_name);
	  for ( int p = 0;  p < work.dim5*cov_ppi.N_part;  p++ )
	  {
		work.space5[p] = new realnum**[Gauss_number];				alloc_test(work.space5[p],func_name);
  		for ( int k=0;  k<Gauss_number;  k++ )
  		{
		  work.space5[p][k] = new realnum*[Num_p];					alloc_test(work.space5[p][k],func_name);
		  for ( int i=0;  i<Num_p;  i++ )
		  {  
			alloc_test_init( work.space5[p][k][i], realnum, Num_p, func_name);
//			work.space5[p][k][i] = new realnum[Num_p];				alloc_test(work.space5[p][k][i],func_name);
		  }
		}
	  }
	}
  }

  if ( do_der_Fph )
  {
	alloc_test_init( der_F, realnum, Num_p, func_name )
//    der_F = new realnum[Num_p]; 									alloc_test(der_F,func_name);
	alloc_test_init( der_ph, realnum, Num_p, func_name )
//    der_ph = new realnum[Num_p];  									alloc_test(der_ph,func_name);
//	if (no_integ_last)
	{
	  alloc_test_init( auxder_F, realnum, Num_p, func_name )
//  	  auxder_F = new realnum[Num_p]; 								alloc_test(auxder_F,func_name);
	  alloc_test_init( auxder_ph, realnum, Num_p, func_name )
//  	  auxder_ph = new realnum[Num_p];  								alloc_test(auxder_ph,func_name);
	}
  }
  if ( do_der_Fph == 2 )
  {
    der2_F2 = new realnum*[Num_p]; 									alloc_test(der2_F2,func_name);
    der2_ph2 = new realnum*[Num_p];  								alloc_test(der2_ph2,func_name);
    der2_phF = new realnum*[Num_p]; 								alloc_test(der2_phF,func_name);
//	if (no_integ_last)
	{
  	  auxder2_F2 = new realnum*[Num_p]; 							alloc_test(auxder2_F2,func_name);
  	  auxder2_ph2 = new realnum*[Num_p];  							alloc_test(auxder2_ph2,func_name);
  	  auxder2_phF = new realnum*[Num_p]; 							alloc_test(auxder2_phF,func_name);
	}
    for ( int i = 0;  i < Num_p;  i++ ) 
    {
	  alloc_test_init( der2_F2[i], realnum, Num_p, func_name )
//  	  der2_F2[i] = new realnum[Num_p]; 								alloc_test(der2_F2[i],func_name);
	  alloc_test_init( der2_ph2[i], realnum, Num_p, func_name )
//  	  der2_ph2[i] = new realnum[Num_p];  							alloc_test(der2_ph2[i],func_name);
	  alloc_test_init( der2_phF[i], realnum, Num_p, func_name )
//  	  der2_phF[i] = new realnum[Num_p]; 							alloc_test(der2_phF[i],func_name); 
//	  if (no_integ_last)
	  {
		alloc_test_init( auxder2_F2[i], realnum, Num_p, func_name )
//    	auxder2_F2[i] = new realnum[Num_p]; 						alloc_test(auxder2_F2[i],func_name);
		alloc_test_init( auxder2_ph2[i], realnum, Num_p, func_name )
//    	auxder2_ph2[i] = new realnum[Num_p];  						alloc_test(auxder2_ph2[i],func_name);
		alloc_test_init( auxder2_phF[i], realnum, Num_p, func_name )
//    	auxder2_phF[i] = new realnum[Num_p]; 						alloc_test(auxder2_phF[i],func_name); 
	  }
    }
  }// maybe separate struct for der variables?

  if ( do_der_D )
  {
    der_D = new realnum**[cov_ppi.N_part];							alloc_test(der_D,func_name);
    auxder_D = new realnum**[cov_ppi.N_part];						alloc_test(auxder_D,func_name);
	work.dim8 = 2;
	work.space8 = new realnum**[work.dim8*cov_ppi.N_part];			alloc_test(work.space8,func_name);
	for ( int p = 0;  p < work.dim8*cov_ppi.N_part;  p++ )
	{
  	  if ( p<cov_ppi.N_part ) { der_D[p] = new realnum*[Num_p]; 	alloc_test(der_D[p],func_name);	}
  	  if ( p<cov_ppi.N_part ) { auxder_D[p] = new realnum*[Num_p]; 	alloc_test(auxder_D[p],func_name);	}
	  work.space8[p] = new realnum*[Num_p]; 						alloc_test(work.space8[p],func_name);
	  for ( int i = 0;  i < Num_p;  i++ ) 
	  {
    	if ( p<cov_ppi.N_part ) alloc_test_init( der_D[p][i], realnum, Num_p, func_name )
//    	{ der_D[p][i] = new realnum[Num_p]; alloc_test(der_D[p][i],func_name);	}
    	if ( p<cov_ppi.N_part ) alloc_test_init( auxder_D[p][i], realnum, Num_p, func_name )
//    	{ auxder_D[p][i] = new realnum[Num_p]; alloc_test(auxder_D[p][i],func_name);	}
		alloc_test_init( work.space8[p][i], realnum, Num_p, func_name )
//    	work.space8[p][i] = new realnum[Num_p]; 					alloc_test(work.space8[p][i],func_name);
	  }
	}
  }

  if ( do_der_D == 2 )
  {
    der2_D = new realnum*****[cov_ppi.N_part]; 						alloc_test(der2_D,func_name);
    auxder2_D = new realnum*****[cov_ppi.N_part]; 					alloc_test(auxder2_D,func_name);
	for ( int p = 0;  p < cov_ppi.N_part;  p++ )
	{
  	  der2_D[p] = new realnum****[Num_p]; 							alloc_test(der2_D[p],func_name);
  	  auxder2_D[p] = new realnum****[Num_p]; 						alloc_test(auxder2_D[p],func_name);
	  for ( int i = 0;  i < Num_p;  i++ )
	  {
    	der2_D[p][i] = new realnum***[Num_p]; 						alloc_test(der2_D[p][i],func_name);
    	auxder2_D[p][i] = new realnum***[Num_p]; 					alloc_test(auxder2_D[p][i],func_name);
		for ( int j = 0;  j < Num_p;  j++ )
		{
		  der2_D[p][i][j] = new realnum**[cov_ppi.N_part]; 			alloc_test(der2_D[p][i][j],func_name);
		  auxder2_D[p][i][j] = new realnum**[cov_ppi.N_part]; 		alloc_test(auxder2_D[p][i][j],func_name);
		  for ( int p2 = 0;  p2 < cov_ppi.N_part;  p2++ ) 
		  {
			der2_D[p][i][j][p2] = new realnum*[Num_p];				alloc_test(der2_D[p][i][j][p2],func_name);
			auxder2_D[p][i][j][p2] = new realnum*[Num_p];			alloc_test(auxder2_D[p][i][j][p2],func_name);
			for ( int k = 0;  k < Num_p;  k++ ) 
			{
			  alloc_test_init( der2_D[p][i][j][p2][k], realnum, Num_p, func_name )
//			  der2_D[p][i][j][p2][k] = new realnum[Num_p];			alloc_test(der2_D[p][i][j][p2][k],func_name);
			  alloc_test_init( auxder2_D[p][i][j][p2][k], realnum, Num_p, func_name )
//			  auxder2_D[p][i][j][p2][k] = new realnum[Num_p];		alloc_test(auxder2_D[p][i][j][p2][k],func_name);
			}
		  }
		}
	  }
	}

	work.dim4 = 4; 	if ( N_meas_p >= 2 ) work.dim4 = 6;
	work.space4 = new realnum***[work.dim4*cov_ppi.N_part];			alloc_test(work.space4,func_name);
    for ( int p = 0;  p < work.dim4*cov_ppi.N_part;  p++ ) 
    {
	  work.space4[p] = new realnum**[Num_p]; 						alloc_test(work.space4[p],func_name);
	  for ( int i = 0;  i < Num_p;  i++ ) 
	  {
		work.space4[p][i] = new realnum*[Num_p]; 					alloc_test(work.space4[p][i],func_name);
		for ( int j = 0;  j < Num_p;  j++ ) 
		{
		  alloc_test_init( work.space4[p][i][j], realnum, Num_p, func_name )
//		  work.space4[p][i][j] = new realnum[Num_p]; 				alloc_test(work.space4[p][i][j],func_name);
		}
	  }
    }
  }

  work.dim1 = 12; if ( N_meas_p>=2 ) work.dim1 = 24; else if ( do_der_D ) work.dim1 = 16;
  work.dim2 = 7; //if ( do_der_D ) work.dim2 = 7 ;
  alloc_test_init( work.space1, realnum, work.dim1*Num_p, func_name )
//  work.space1 = new realnum[work.dim1*Num_p];						alloc_test(work.space1,func_name);
  work.space2 = new realnum*[work.dim2*Num_p];						alloc_test(work.space2,func_name);
  for ( int i = 0;  i < work.dim2*Num_p;  i++ ) 
  {
	alloc_test_init( work.space2[i], realnum, Num_p, func_name )
//    work.space2[i] = new realnum[Num_p];							alloc_test(work.space2[i],func_name);
  }    
  a = new realnum*[Num_p];    	b = new realnum*[Num_p];			alloc_test(a,func_name);	alloc_test(b,func_name);
  c = new realnum*[Num_p]; 	 	d = new realnum*[Num_p];			alloc_test(c,func_name);	alloc_test(d,func_name);
  for ( int i = 0;  i < Num_p;  i++ )   
  {
	alloc_test_init( a[i], realnum, Num_p, func_name )
	alloc_test_init( b[i], realnum, Num_p, func_name )
	alloc_test_init( c[i], realnum, Num_p, func_name )
	alloc_test_init( d[i], realnum, Num_p, func_name )
//    a[i] = new realnum[Num_p];    b[i] = new realnum[Num_p];		alloc_test(a[i],func_name); alloc_test(b[i],func_name);
//    c[i] = new realnum[Num_p];    d[i] = new realnum[Num_p];		alloc_test(c[i],func_name); alloc_test(d[i],func_name);
  }

//  for ( int i = 0;  i < Num_p;  i++ )
//  {
	alloc_test_init( eigenvalues1[0], realnum, Num_p, func_name )
//	eigenvalues1[0] = new realnum[Num_p];
	alloc_test_init( eigenvalues2[0], realnum, N_mod, func_name )
//	eigenvalues2[0] = new realnum[N_mod];
	alloc_test_init( eigenvalues1[1], realnum, 2, func_name )
//	eigenvalues1[1] = new realnum[2];
	alloc_test_init( eigenvalues2[1], realnum, 1, func_name )
//	eigenvalues2[1] = new realnum[1];
	alloc_test_init( inv_eigenval, realnum, Num_p, func_name )
//	inv_eigenval = new realnum[Num_p];
//  }

  if ( Tab_func )
	tab = new TabFunc<realnum>;

}

#undef alloc_test_init


template<typename realnum>
void likelihood<realnum>::SetNoIntegLast ( int nointlast )
{
  char* func_name = (char*) "likelihood::SetNoIntegLast";
  no_integ_last = nointlast;
}

template<typename realnum>
void likelihood<realnum>::SetCentRice( int centricity, int rice )
{
  cent = centricity;
  Rice = rice;
  if ( cov_ppi.Dtype == "sad" || cov_ppi.Dtype == "sadh"  
//	|| cov_ppi.Dtype == "sras" || cov_ppi.Dtype == "sras3d" || cov_ppi.Dtype == "sras4d" ) 
  )
	  if ( cent ) Rice = 1;
  if ( cov_ppi.Dtype == "no" ) Rice = 1;
  if ( Rice == -1 ) Rice = 0;

  if ( Rice )
  {
	Num = 2;
	N_meas = 1;
	cov = &cov_rice;
  }
  else
  {
	Num = Num_p;
	N_meas = N_meas_p;
	cov = &cov_ppi;
  }
}


// Basic constructor

template<typename realnum>
likelihood<realnum>::likelihood( int N_measur, int Number, int dFph, int dD, const string ppiDtype, 
  int N_part, int Gauss_Number_Of_Nodes,  int Gauss_Number_Of_Nodes_2,
  int no_imaginary, const int Tabul_functions, Gauss_structure<realnum>* G_s_copy )
  : Num_p(Number), N_meas_p(N_measur), cov_ppi(this, N_part, Number, N_measur, ppiDtype, dD, no_imaginary), 
	cov_rice(this, N_part, 2, 1, "single", dD, no_imaginary)
{ 
  char* func_name = (char*) "likelihood::likelihood";
  SetCentRice( 0, -1 );
  Tab_func = Tabul_functions;
  if ( ppiDtype == "sras" || ppiDtype == "sras3d" || ppiDtype == "mad" ) 
	no_imag = 1; // imaginary terms not supported for SIRAS or MAD
  else no_imag = no_imaginary;
  no_integ_last = 0;
  if ( ( ppiDtype == "sras" ) && N_part < 2 ) {
	N_part = 2;
	likelihood<realnum>::Error( 202, func_name, (void*) &N_part );
  }
  if ( ( ppiDtype == "sras3d" ) && N_part < 3 ) {
	N_part = 3;
	likelihood<realnum>::Error( 202, func_name, (void*) &N_part );
  }
  checkOK = 0;
  max_num_points_allowed = 1e100;
  safer_num_points = 0;
  do_der_Fph = dFph;
  do_der_D = dD;
  do_ABCD = 0;
  if ( do_der_Fph > 2 || do_der_D > 2 || do_der_Fph < 0 || do_der_D < 0 )  Error( 107, func_name ); 
  if (  ( do_der_D == 2 ) && 
		( ppiDtype == "sras" || ppiDtype == "sras3d" || ppiDtype == "mad" ) )
		  Error( 113, func_name, (void*) &ppiDtype ); 
  min_acc_eigen = 0.001;
  min_nonzero = 1e-14;
  max_hyparg = 50;
  Bess_prec1 = 13;
  Bess_prec2 = 19;
  if ( no_integ_last ) nil_m = -1;
  int N_mod = Num_p - N_meas_p;
  if ( N_mod<1 ) Error( 102, func_name );
  if ( N_meas_p >= 2 )	Gauss_number = Gauss_Number_Of_Nodes; 
  if ( N_meas_p >= 3 )	Gauss_number_2 = Gauss_Number_Of_Nodes_2; 
  if ( N_meas_p > 2 ) improve_integ = 1;   else improve_integ = 0;
  num4phib = -1;
  lpck.Alloc(Num_p,N_mod,abs(no_imag-1));
  Alloc( G_s_copy );
}


// Copy constructor

template<typename realnum>
likelihood<realnum>::likelihood( likelihood<realnum>& l_in ) 
  : Num_p(l_in.Num), N_meas_p(l_in.N_meas),
	cov_ppi(this, l_in.cov_ppi.N_part, l_in.cov_ppi.Num, l_in.cov_ppi.N_meas, l_in.cov_ppi.Dtype, l_in.do_der_D, l_in.no_imag), 
	cov_rice(this, l_in.cov_rice.N_part, l_in.cov_rice.Num, l_in.cov_rice.N_meas, l_in.cov_rice.Dtype, l_in.do_der_D, l_in.no_imag)
{ 
  char* func_name = (char*) "likelihood::likelihood(likelihood)";
  SetCentRice( l_in.cent, l_in.Rice );
  Tab_func = l_in.Tab_func;
  no_imag = l_in.no_imag;
  no_integ_last =  l_in.no_integ_last;
  checkOK = l_in.checkOK;
  max_num_points_allowed = l_in.max_num_points_allowed;
  safer_num_points = l_in.safer_num_points;
  do_der_Fph = l_in.do_der_Fph;
  do_der_D = l_in.do_der_D;
  do_ABCD = l_in.do_ABCD;
  if ( do_der_Fph > 2 || do_der_D > 2 || do_der_Fph < 0 || do_der_D < 0 )  Error( 107, func_name ); 
  min_acc_eigen = l_in.min_acc_eigen;
  min_nonzero = l_in.min_nonzero;
  max_hyparg = l_in.max_hyparg;
  Bess_prec1 = l_in.Bess_prec1;
  Bess_prec2 = l_in.Bess_prec2;
  if ( no_integ_last ) nil_m = l_in.nil_m;
  int N_mod = Num_p - N_meas_p;
  if ( N_mod<1 ) Error( 102, func_name );
  if ( N_meas_p >= 2 )	Gauss_number = l_in.Gauss_number;
  if ( N_meas_p >= 3 )	Gauss_number_2 = l_in.Gauss_number_2;
  improve_integ = l_in.improve_integ;
  num4phib = l_in.num4phib;
  lpck = l_in.lpck;
  Alloc( l_in.Gauss_struct );
}


// Assignment operator

template<typename realnum>
likelihood<realnum>& likelihood<realnum>::operator=( likelihood<realnum>& l_in ) 
{ 
  char* func_name = (char*) "likelihood::operator=";
  if ( a != 0 )
	DeAlloc();
  cov_ppi = l_in.cov_ppi;
  cov_ppi.llh = this;
  cov_rice = l_in.cov_rice;
  cov_rice.llh = this;
  Num = l_in.Num;
  Num_p = l_in.Num_p;
  N_meas_p = l_in.N_meas_p;
  SetCentRice( l_in.cent, l_in.Rice );
  Tab_func = l_in.Tab_func;
  no_imag = l_in.no_imag;
  no_integ_last =  l_in.no_integ_last;
  checkOK = l_in.checkOK;
  max_num_points_allowed = l_in.max_num_points_allowed;
  safer_num_points = l_in.safer_num_points;
  do_der_Fph = l_in.do_der_Fph;
  do_der_D = l_in.do_der_D;
  do_ABCD = l_in.do_ABCD;
  if ( do_der_Fph > 2 || do_der_D > 2 || do_der_Fph < 0 || do_der_D < 0 )  Error( 107, func_name ); 
  min_acc_eigen = l_in.min_acc_eigen;
  min_nonzero = l_in.min_nonzero;
  max_hyparg = l_in.max_hyparg;
  Bess_prec1 = l_in.Bess_prec1;
  Bess_prec2 = l_in.Bess_prec2;
  if ( no_integ_last ) nil_m = l_in.nil_m;
  int N_mod = Num_p - N_meas_p;
  if ( N_mod<1 ) Error( 102, func_name );
  if ( N_meas_p >= 2 )	Gauss_number = l_in.Gauss_number;
  if ( N_meas_p >= 3 )	Gauss_number_2 = l_in.Gauss_number_2;
  improve_integ = l_in.improve_integ;
  lpck = l_in.lpck;
  num4phib = l_in.num4phib;
  Alloc( l_in.Gauss_struct );
  return *this;
}


//#undef alloc_test





// Computes inverse matrix to Hermitian complex or real symmetrix matrix Inp using precomputed eigenvalues
// result is stored in matrices 'a' (real part) and 'b' (imag. part) 
// only upper triangle of Inp is required as input

template<typename realnum>
inline void likelihood<realnum>::InverseMatrix( unsigned N, realnum** Inp_Re, realnum** Inp_Im, 
	realnum **Inver_Re, realnum **Inver_Im, realnum* Eigenvalues, realnum* Inv_Eigenval )
{

  if (cov->no_imag)
  {
	for (int  i = 0;  i < N;  i++) {
	  for (int  j = i;  j < N;  j++) {
    	realnum Aux_Re = 0;
    	for (int  k = 0;  k < N;  k++) 	if ( Eigenvalues[k] > min_acc_eigen )
      	  Aux_Re += Inp_Re[k][i]*Inp_Re[k][j] * Inv_Eigenval[k];
    	Inver_Re[i][j] = Inver_Re[j][i] = Aux_Re;
    	Inver_Im[i][j] = 0;
  	  }
	}
  }
  else
  {
	for (int  i = 0;  i < N;  i++) {
	  for (int  j = i;  j < N;  j++) {
    	realnum Aux_Re = 0, Aux_Im = 0;
    	for (int  k = 0;  k < N;  k++) 	if ( Eigenvalues[k] > min_acc_eigen )
    	{
      	  Aux_Re += ( Inp_Re[k][i]*Inp_Re[k][j] + Inp_Im[k][i]*Inp_Im[k][j] ) * Inv_Eigenval[k];
      	  Aux_Im += ( Inp_Re[k][i]*Inp_Im[k][j] - Inp_Im[k][i]*Inp_Re[k][j] ) * Inv_Eigenval[k];
    	}
    	Inver_Re[i][j] = Inver_Re[j][i] = Aux_Re;
    	Inver_Im[i][j] = Aux_Im;   
		Inver_Im[j][i] = -Aux_Im;
  	  }
	}
  }

}


#define exchange(x,y,aux) { aux=x; x=y; y=aux; }
#define make_w_pos(i,j,m,n,c) { w[i][j]*=-1; w[m][n]*=-1; if (i==m) c=i; else if (j==n) c=j; else c=1; phw[c]+=PI; }


// calculate the w coefficients (if improv_integ is set) and use them to initialize the values for integration improvement

template<typename realnum>
void likelihood<realnum>::InitIntegImpr( realnum *F, realnum *phase )
{
  sum_w = w_crit = w_crit_2 = num_changing = 0;
  realnum *Wc = work.space1 + 4*Num, *Ws = work.space1 + 5*Num; // temporary use of allocated space
  for (int j=0; j<N_meas; j++)
  {
	Wc[j] = Ws[j] = 0; 
	for (int i=j+1; i<Num; i++)
	{
	  w[j][i] = -2*F[i]*F[j]*a[j][i]; 
	  if (i>=N_meas) {
		if (Tab_func) Wc[j] += w[j][i]*tab->Cos(phase[i]);
		else			Wc[j] += w[j][i]*cos(phase[i]);
		if (Tab_func) Ws[j] += w[j][i]*tab->Sin_charged(phase[i]);
		else			Ws[j] += w[j][i]*sin(phase[i]);
	  }
	}
	w[j][j] = sqrt(Wc[j]*Wc[j]+Ws[j]*Ws[j]);
	if ( w[j][j] != 0. )  phw[j] = acos(Wc[j]/w[j][j]);
	if ( Ws[j] < 0. )  phw[j] *= -1;
  }
// make w's positive - all up to the smallest w[i][j],i!=j (all w[i][i]'s already are positive!)
// and find the max. w and calculate sum_w
  realnum maxw = 0, minw = 1e20, midw;
  realnum W_max=0;
  int W_max_ind=0;
  int minw_ind[2];
  int num_neg=0;
  int neg[3][2];
  for (int j=0; j<N_meas; j++)  
	for (int i=j; i<N_meas; i++)
	{ 
	  realnum w_abs = w[j][i] ;   
	  if ( w_abs<0. ) { w_abs *= -1;  neg[num_neg][0]=j; neg[num_neg][1]=i; num_neg++; } 
	  sum_w += w_abs;
	  if ( i!=j ) {
		if ( w_abs > maxw ) { maxw_ind[0] = j; maxw_ind[1] = i; maxw = w_abs; }
		if ( w_abs < minw ) { minw_ind[0] = j; minw_ind[1] = i; minw = w_abs; }
	  }
	  else if ( w_abs > W_max ) { W_max_ind = j; W_max = w_abs; }
	}
// remember if wmax was negative  
  int neg_wmax = 0;
  if ( w[maxw_ind[0]][maxw_ind[1]]<0 ) neg_wmax = 1;
  if ( debug_info == 1 ) {
	cout << "W's: " << w[0][0] << " " << w[1][1] << " "  << w[2][2] << endl;
	cout << "phi's: " << phw[0] << " " << phw[1] << " "  << phw[2] << endl;
	cout << "w's: " << w[0][1] << " " << w[0][2] << " "  << w[1][2] << endl;
  }	
// here making the w's positive (up to the smallest w[i][j])
  int common_i=0;
  if ( num_neg==1 && (neg[0][0]!=minw_ind[0]||neg[0][1]!=minw_ind[1]) )  {
    make_w_pos(neg[0][0],neg[0][1],minw_ind[0],minw_ind[1],common_i)
  } else if ( num_neg==2 )  {
	make_w_pos(neg[0][0],neg[0][1],neg[1][0],neg[1][1],common_i)
  } else if ( num_neg==3 )  {
	for (int i=0; i<2; i++)   
	  if (neg[i][0]==minw_ind[0]&&neg[i][1]==minw_ind[1]) 
		{ exchange(neg[i][0],neg[2][0],common_i)  exchange(neg[i][1],neg[2][1],common_i) }
	make_w_pos(neg[0][0],neg[0][1],neg[1][0],neg[1][1],common_i)
  }
  if ( w[minw_ind[0]][minw_ind[1]] < 0 ) sum_w -= 2*minw;
  realnum Wmi0 = w[maxw_ind[0]][maxw_ind[0]], Wmi1 = w[maxw_ind[1]][maxw_ind[1]];
  realnum vert_term = 0;
// find out whether we can localize maximum, if so then do that (in case that improve_integ>=2)
// and estimate the height of the peak and do the automatic resampling for n (the last variable)
  if ( improve_integ_now && improve_integ_now < 100 )
  {
	for (int j=0; j<N_meas; j++)  
	  for (int i=j; i<N_meas; i++)
	  { 
		if ( i==j && j!=W_max_ind ) 
		  if (Tab_func) 	w_crit += w[j][i]*( 1-tab->Cos(phw[W_max_ind]-phw[j]) ); // don't need fabs since w[i][i] >= 0 
		  else 				w_crit += w[j][i]*( 1-cos(phw[W_max_ind]-phw[j]) );
		else if ( w[j][i]!=maxw && w[j][i]!=minw )   midw=w[j][i];
	  }
	w_crit_2 = midw*minw + maxw*(minw+midw);
	w_crit_2 = w_crit_2*w_crit_2 - 2*maxw*midw*minw*w_crit ;
    if ( w_crit_2 > 0 )  w_crit_2=sqrt(w_crit_2);  else  w_crit_2=1e100;
    if ( maxw>0 ) w_crit_2 = ((maxw+midw)*(maxw+minw)-w_crit*maxw-w_crit_2)/(maxw*maxw);
	if ( debug_info == 1 ) cout << "w_crit_2,acos(wcrit_2): " << w_crit_2 << " " << acos(w_crit_2) << endl;
	if ( w[minw_ind[0]][minw_ind[1]] < 0 ) w_crit += 2*minw;
	if ( maxw>0 ) w_crit = 1-w_crit/maxw;
	if ( w_crit_2<-0.5 || w[maxw_ind[0]][maxw_ind[1]]<4 ) improve_integ_now = 1;
	if ( improve_integ_now >= 3 )
	  if (Tab_func) 	vert_term = sqrt(Wmi0*Wmi0+Wmi1*Wmi1+2*Wmi0*Wmi1*tab->Cos(phw[maxw_ind[0]]-phw[maxw_ind[1]])) ;
	  else				vert_term = sqrt(Wmi0*Wmi0+Wmi1*Wmi1+2*Wmi0*Wmi1*cos(phw[maxw_ind[0]]-phw[maxw_ind[1]])) ;
	if ( w_crit_2 > 0.75 && improve_integ_now >= 3 )   sum_w += vert_term -Wmi0 -Wmi1;
	if ( debug_info == 1 ) cout << "w_crit,acos(wcrit),sum_w: " << w_crit << " " << acos(w_crit) << " " << sum_w << endl;
  }
// now when we know improve_integ_now, reinitialize Gauss2 cos and sin if we it is needed
  if ( improve_integ_now > 1 ) {
	cos_last2 = Gauss_struct[1].cos2;
	sin_last2 = Gauss_struct[1].sin2;	  
  }
// automatic resampling for n (the last variable)
  int new_num_points = 0;
  if ( improve_integ_now >= 3 || improve_integ_now == 1 ) 
  {	
	if ( sum_w < 4 ) new_num_points = 8;
	else if ( sum_w < 6 ) new_num_points = 12;
	else if ( sum_w < 7.5 ) new_num_points = 15;
	else if ( sum_w < 12 ) new_num_points = 20;
	else if ( sum_w < 20 ) new_num_points = 25;
	else if ( sum_w < 35 ) new_num_points = 35;
	else if ( sum_w < 70 ) new_num_points = 50;
	else if ( sum_w < 125 ) new_num_points = 70;
	else if ( sum_w < 235 ) new_num_points = 100;
	else if ( sum_w < 410 ) new_num_points = 135;
	else if ( sum_w < 800 ) new_num_points = 180;
	else if ( sum_w < 1700 ) new_num_points = 250;
	else if ( sum_w < 4000 ) new_num_points = 350;
	else if ( sum_w < 10000 ) new_num_points = 500;
// increasing the number of points as given by safer_num_points
	if ( safer_num_points > 0 )
	  for (int i=num_Gauss_struct-1-safer_num_points; i>=0; i--)
		if ( new_num_points == Gauss_struct[i].num_points ) 
		  new_num_points = Gauss_struct[i+safer_num_points].num_points;
	for (int i=2; i<num_Gauss_struct; i++)  
	  if (new_num_points==Gauss_struct[i].num_points && new_num_points<Gauss_number_2) {
		Gauss_number_now_2 = Gauss_struct[i].num_points;
		if ( improve_integ_now == 1 ) {
		  cos_last2 = Gauss_struct[i].cos;
		  sin_last2 = Gauss_struct[i].sin;
		} else {
		  cos_last2 = Gauss_struct[i].cos2;
		  sin_last2 = Gauss_struct[i].sin2;
		}
		Gauss_weight_2 = Gauss_struct[i].weight;
	  }
// if negative wmax then shift the n values so that the maximum is in the middle of n-interval
// (this is not ideal and actully very rarely used and it could be omitted in principle)

	if (neg_wmax && improve_integ_now >= 3) {
	  realnum *aux_c = work.space9, *aux_s = work.space9 + Gauss_number_now_2, *aux_G = work.space9 + Gauss_number_now_2*2;
	  int G2_half=Gauss_number_now_2/2, G2_mod=Gauss_number_now_2%2;
	  for (int n=0; n<G2_half; n++) { 
		aux_c[n+G2_half+G2_mod]=cos_last2[n];  aux_s[n+G2_half+G2_mod]=sin_last2[n];  aux_G[n+G2_half+G2_mod]=Gauss_weight_2[n];
	  }
	  for (int n=G2_half; n<Gauss_number_now_2; n++) { 
		aux_c[n-G2_half]=cos_last2[n];  aux_s[n-G2_half]=sin_last2[n];  aux_G[n-G2_half]=Gauss_weight_2[n];  
	  }
	  cos_last2=aux_c; sin_last2=aux_s; Gauss_weight_2=aux_G;
	}
// estimate the height of the peak and do the automatic resampling for m (the second variable) 
	new_num_points = 0;
	sum_wm = 1e30;
	if ( improve_integ_now >= 3 && w_crit>0 /*&& w[maxw_ind[0]][maxw_ind[1]]>5.*/ ) 
	{  
	  realnum diag_term=0, horiz_term=0, vert_ph=0, horiz_ph;
	  if ( vert_term != 0. )  {
		if (Tab_func) 	vert_ph = (Wmi0*tab->Cos(phw[maxw_ind[0]])+Wmi1*tab->Cos(phw[maxw_ind[1]]))/vert_term;
		else 			vert_ph = (Wmi0*cos(phw[maxw_ind[0]])+Wmi1*cos(phw[maxw_ind[1]]))/vert_term;
		if ( vert_ph < -1+1e-4 ) vert_ph=PI;
		else if ( vert_ph > 1-1e-4 ) vert_ph=0;
		else vert_ph = acos(vert_ph);
	  }
	  if ( Tab_func && Wmi0*tab->Sin(phw[maxw_ind[0]])+Wmi1*tab->Sin(phw[maxw_ind[1]]) < 0. )  vert_ph *= -1;
	  else if ( !Tab_func && Wmi0*sin(phw[maxw_ind[0]])+Wmi1*sin(phw[maxw_ind[1]]) < 0. )  vert_ph *= -1;
		  
	  for (int i=0; i<N_meas; i++) 
		for (int j=i+1; j<N_meas; j++) 
		  if ( i!=maxw_ind[0] || j!=maxw_ind[1] ) 
			diag_term += w[i][j];
	  if (diag_term<0) diag_term*=-1;
	  for (int i=0; i<N_meas; i++) 
		if ( i!=maxw_ind[0] && i!=maxw_ind[1] ) {
		  horiz_term = w[i][i];  horiz_ph=phw[i]; 
		}
	  if ( diag_term>vert_term ) { //sum_wm = vert_term + horiz_term;
		sum_wm = sqrt(vert_term*vert_term+horiz_term*horiz_term+2*vert_term*horiz_term*cos(vert_ph-horiz_ph)) ;
		realnum diff_here = fabs(vert_term-horiz_term);
		if (diff_here/horiz_term<0.04) sum_wm*=1.5;
		else if (diff_here/horiz_term<0.06) sum_wm*=1.3;
		else if (diff_here/horiz_term<0.10) sum_wm*=1.1;
	  } 
	  else {
		if ( fabs(diag_term-horiz_term)/horiz_term < 0.2 )  
		  sum_wm = diag_term + horiz_term;
		else
		  sum_wm = sqrt(diag_term*diag_term+horiz_term*horiz_term+2*diag_term*horiz_term*cos(vert_ph-horiz_ph)) ;
	  }
	  if ( (w_crit<.5&&sum_wm*10<Wmi0+Wmi1) || (w_crit<.9&&sum_wm*30<Wmi0+Wmi1) ) sum_wm *= 2;
	  if ( debug_info == 1 ) {
		cout << "terms(diag,vert,horiz)" << diag_term << " " << vert_term << " " << horiz_term << endl;
		cout << "ph(vert, horiz): " << vert_ph << " " << horiz_ph << endl;
	  }
	  if ( sum_wm < 0.75 ) new_num_points = 12;
	  else if ( sum_wm < 2.2 ) new_num_points = 15;
	  else if ( sum_wm < 5 ) new_num_points = 20;
	  else if ( sum_wm < 9.5 ) new_num_points = 25;
	  else if ( sum_wm < 20 ) new_num_points = 35;
	  else if ( sum_wm < 50 ) new_num_points = 50;
	  else if ( sum_wm < 120 ) new_num_points = 70;
	}
	else { // if improve_integ_now == 1 || w_crit<0
	  sum_wm = sum_w - w[maxw_ind[0]][maxw_ind[1]];
	  if ( sum_wm < 2 ) new_num_points = 8;
	  else if ( sum_wm < 3.5 ) new_num_points = 12;
	  else if ( sum_wm < 5 ) new_num_points = 15;
	  else if ( sum_wm < 8 ) new_num_points = 20;
	  else if ( sum_wm < 13 ) new_num_points = 25;
	  else if ( sum_wm < 28 ) new_num_points = 35;
	  else if ( sum_wm < 60 ) new_num_points = 50;
	  else if ( sum_wm < 150 ) new_num_points = 70;  
	}
// increasing the number of points as given by safer_num_points
	if ( safer_num_points > 0 )
	  for (int i=num_Gauss_struct-1-safer_num_points; i>=0; i--)
		if ( new_num_points == Gauss_struct[i].num_points ) 
		  new_num_points = Gauss_struct[i+safer_num_points].num_points;
	if ( debug_info == 1 ) cout << "sum_wm: " << sum_wm << endl;
	for (int i=2; i<num_Gauss_struct; i++)  
	  if (new_num_points==Gauss_struct[i].num_points && new_num_points<Gauss_number) {
		Gauss_number_now = Gauss_struct[i].num_points;
		cos_last = Gauss_struct[i].cos;
		sin_last = Gauss_struct[i].sin;
		Gauss_weight = Gauss_struct[i].weight;
	  }
  }
}

#undef make_w_pos


// changes the order of the variables in order to improve the integration

template<typename realnum>
void likelihood<realnum>::ChangeIntegOrder( realnum *F, realnum *phase, int back )
{
  int k_start=0, k_stop=num_changing;
  if (back) { k_start=num_changing-1; k_stop=0; }
  for (int k=k_start; back ? k>=k_stop : k<k_stop ; back ? k-- : k++)
	if ( maxw_ind[k] > k )
	{
	  realnum aux;
	  exchange(F[k],F[maxw_ind[k]],aux)
	  exchange(phase[k],phase[maxw_ind[k]],aux)
//	  exchange(eigenvalues1[0][k],eigenvalues1[0][maxw_ind[k]],aux)
//	  exchange(eigenvalues2[0][k],eigenvalues2[0][maxw_ind[k]],aux)
//	  exchange(inv_eigenval[k],inv_eigenval[maxw_ind[k]],aux)
	  exchange(a[k][k],a[maxw_ind[k]][maxw_ind[k]],aux)
	  for (int i=0; i<Num; i++)  
		if (i!=k && i!=maxw_ind[k]) { exchange(a[k][i],a[maxw_ind[k]][i],aux)   exchange(a[i][k],a[i][maxw_ind[k]],aux) }
/*	  exchange(cov->re[k][k],cov->re[maxw_ind[k]][maxw_ind[k]],aux)
	  for (int i=0; i<Num; i++)  
		if (i!=k && i!=maxw_ind[k]) { exchange(cov->re[k][i],cov->re[maxw_ind[k]][i],aux) 
									  exchange(cov->re[i][k],cov->re[i][maxw_ind[k]],aux) 
							}
	  exchange(cov->im[k][k],cov->im[maxw_ind[k]][maxw_ind[k]],aux)
*/	  
	  if ( do_der_D)
	  for (int p=0; p<cov->N_part; p++)
	  {
		// checking whether the bottom triangle deriv. terms were assigned (once is enough so only done for k==0)
		if ( k==0 && !back )
		for (int i=0; i<Num; i++)  
		  for (int j=i+1; j<Num; j++)  
			if ( cov->part[p].deriv_D_r[i][j] != cov->part[p].deriv_D_r[j][i] )
			  cov->part[p].deriv_D_r[j][i] = cov->part[p].deriv_D_r[i][j];
		exchange(cov->part[p].deriv_D_r[k][k],cov->part[p].deriv_D_r[maxw_ind[k]][maxw_ind[k]],aux)
		for (int i=0; i<Num; i++)  
		{
		  if (i!=k && i!=maxw_ind[k]) { exchange(cov->part[p].deriv_D_r[k][i],cov->part[p].deriv_D_r[maxw_ind[k]][i],aux)
										exchange(cov->part[p].deriv_D_r[i][k],cov->part[p].deriv_D_r[i][maxw_ind[k]],aux)
									  }
		}
#define to(i,j) to_array[2*i+j]
		int aux_int;
		for (int i=0; i<cov->Dass[p][0][1].toNum; i++)
		{
		  for (int j=0; j<2; j++)
		  {
			if (cov->Dass[p][0][1].to(i,j)==k)  cov->Dass[p][0][1].to(i,j)=maxw_ind[k];
			else if (cov->Dass[p][0][1].to(i,j)==maxw_ind[k])  cov->Dass[p][0][1].to(i,j)=k;
		  }
		  if ( cov->Dass[p][0][1].to(i,0) > cov->Dass[p][0][1].to(i,1) ) 
			exchange(cov->Dass[p][0][1].to(i,0),cov->Dass[p][0][1].to(i,1),aux_int);
		}
#undef to_array
	  }
	}
}


// Calculates -ln of conditional probability for a single reflection

// Actually for both N_meas=1,2 calculates the -ln llhood target, it's 1. and 2. derivatives wrt F,phase
// and 1. and 2. derivatives wrt D, only first derivatives for N_meas>=3


template<typename realnum>
realnum likelihood<realnum>::EvaluateSR( realnum *F, realnum *phase )
{

  realnum *cos_ph = work.space1,              	*sin_ph = work.space1 + Num;  
  realnum *temp3 = work.space1 +2*Num,          *temp4 = work.space1 + 3*Num; // only used if N_meas == 3 && improve_integ>=2
  realnum *dRdF = work.space1 + 4*Num,   		*dEdF = work.space1 + 5*Num;
  realnum *dRdph = work.space1 + 6*Num,  		*dEdph = work.space1 + 7*Num;
  realnum *F_inv = work.space1 + 8*Num;         phw = work.space1 + 9*Num;  // only used if improve_integ
//  realnum *temp7 = work.space1 + 10*Num,		*temp8 = work.space1 + 11*Num;  // only used if N_meas == 3 && improve_integ>=2
  realnum *dRdF_cos2 = work.space1 + 12*Num,    *dRdph_cos2 = work.space1 + 13*Num;	//only used if improve_integ > 1
  realnum *dRdF_sin2 = work.space1 + 14*Num,    *dRdph_sin2 = work.space1 + 15*Num; //only used if improve_integ > 1
  realnum **storage1 = work.space2,				**storage2 = work.space2 + Num;
  realnum **storage3 = work.space2 + 2*Num,     **storage4 = work.space2 + 3*Num;
  realnum **storage5 = work.space2 + 4*Num,     **storage6 = work.space2 + 5*Num;
  w = work.space2 + 6*Num;
  realnum **sinlast_del = work.space7,			**coslast_del = work.space7 + Gauss_number;
  realnum **sinlast2_del = work.space10,		**coslast2_del = work.space10 + Gauss_number_2;	// only used if N_meas>2

// this is only about storage
#define PCpMS(i,j) storage1[i][j]
#define PSmMC(i,j) storage1[j][i]
#define ACmBS(i,j) storage2[i][j]
#define ASpBC(i,j) storage2[j][i]
#define sin_del(i,j) storage3[i][j]
#define cos_del(i,j) storage3[j][i]
#define secder_cos1(i,j) storage4[i][j]
#define secder_sin1(i,j) storage4[j][i]
#define secder_cos2(i,j) storage5[i][j]
#define secder_sin2(i,j) storage5[j][i]
#define secdermix_cos2(i,j) storage6[i][j]
#define secdermix_sin2(i,j) storage6[j][i]

  
// some precalculations and initializations
  debug_info=0;
  if ( N_meas >= 2 && !no_integ_last ) 
  {	
	Gauss_number_now=Gauss_number; Gauss_number_now_2=Gauss_number_2;
	cos_last = Gauss_struct[0].cos;
	sin_last = Gauss_struct[0].sin;
	Gauss_weight = Gauss_struct[0].weight;
	if ( N_meas > 2 ) {
	  cos_last2 = Gauss_struct[1].cos;
	  sin_last2 = Gauss_struct[1].sin;
	  Gauss_weight_2 = Gauss_struct[1].weight;
	}
  }

// calculate the w coefficients in integrand used for optimization of the integration (if improv_integ is set)
  improve_integ_now = improve_integ;
  if ( N_meas < 3 )  improve_integ_now = 0;	
  if ( improve_integ_now )
  {
	InitIntegImpr(F,phase);
// changing the analytically integrated variable (reverse change is applied at the end of the procedure!)
	num_changing = 1;	// to integrate one important phase out analytically (changing indices maxw_ind[0] and 0)
	if ( improve_integ_now >= 2 )  num_changing ++;	 // also to put second important one to second position (changing maxw_ind[1] and 1)
// forcing first variable to be used for phib calculation
// leads to less efficient calculation so should only be used when necessary
	int save_maxw_ind[2] = {maxw_ind[0],maxw_ind[1]};
	if ( num4phib == 0 )  {
	  maxw_ind[0]=N_meas-1; // aka last
	  if ( maxw_ind[1] == N_meas-1 ) maxw_ind[1] = save_maxw_ind[0];
	}
// forcing second variable to be used for phib calculation
	if ( num4phib == 1 )  {
	  if ( num_changing==1 )  num_changing ++;
	  maxw_ind[1]=N_meas-1; // aka last
	  if ( maxw_ind[0] == 1 ) maxw_ind[0] = save_maxw_ind[1];
	}
	ChangeIntegOrder( F, phase, 0 );
  }
  else if ( num4phib >= 0 && N_meas > 1 )  {
	num_changing = 1;
	maxw_ind[0]=N_meas-1; // aka last
	ChangeIntegOrder( F, phase, 0 );
  }

  int last = N_meas-1;
  int last2 = N_meas-2;		// only used for N_meas > 2
// cos and sin of last phase can be precomputed (acc. to Gauss. integ.) in case of no_integ_last
  if ( no_integ_last )
  { 
	if ( nil_m < 0 || nil_m >= Gauss_number )
	{
	  if (Tab_func) cos_ph[last] = tab->Cos(phase[last]);
	  else			cos_ph[last] = cos(phase[last]);
	  if (Tab_func) sin_ph[last] = tab->Sin_charged(phase[last]);
	  else			sin_ph[last] = sin(phase[last]);
	}
	else
	{
	  cos_ph[last] = cos_last[nil_m];
	  sin_ph[last] = sin_last[nil_m];
	}
  }
  if ( do_ABCD )  HLA = HLB = HLC = HLD = 0.;
  for (int i=N_meas; i<Num; i++) 	
  {
    if (Tab_func) 	{ cos_ph[i] = tab->Cos(phase[i]);   sin_ph[i] = tab->Sin_charged(phase[i]); }
	else			{ cos_ph[i] = cos(phase[i]);   sin_ph[i] = sin(phase[i]); }
    if ( do_der_Fph ) 
    {
      der_F[i] = 0;				der_ph[i] = 0;
      dRdF[i] = 0;           	dEdF[i] = 0;     
      dRdph[i] = 0;          	dEdph[i] = 0;
    }  
    if ( do_der_Fph == 2 ) 
    {
      for (int j=N_meas; j<Num; j++) 
      {
        der2_F2[i][j] = 0;
        der2_phF[i][j] = 0;
		der2_ph2[i][j] = 0;
      }
    }
  }
  for (int i=0; i<Num; i++) 
  {
	if ( F[i] > min_nonzero || F[i] < -min_nonzero ) 	F_inv[i] = 1./F[i]; 	
	  else	F_inv[i] = 1.;
    for (int j=0; j<Num; j++) 
    {
      storage1[i][j] = 0;
      storage2[i][j] = 0;
    }
  }    


// the calculation of llhood value (and derivatives) begins here
  realnum Aux_exp = 0, Aux_Bessel = 0;
  realnum temp1, temp2, temp5, temp6, temp7;
  for (int i=0; i<N_meas; i++)
    Aux_exp += a[i][i]*F[i]*F[i];
  if ( no_integ_last && N_meas == 2 )	Aux_Bessel += F[last]*F[last]*( a[0][last]*a[0][last] + b[0][last]*b[0][last] ) ;
  realnum F2sq_div_F1sq, F2_div_F1, Aux_Bessel_cos2, Aux_Bessel_sin2;	
  if ( improve_integ_now >= 2 ) {
	F2sq_div_F1sq = F[1]*F[1]*F_inv[0]*F_inv[0];
	F2_div_F1 = F[1]*F_inv[0];
	Aux_Bessel_cos2 =  Aux_Bessel_sin2 = 0;
	for (int  i = N_meas;  i < Num;  i++) 
	  dRdF_cos2[i] = dRdF_sin2[i] = dRdph_sin2[i] = dRdph_cos2[i] = 0.;
  }

  for (int  i = N_meas;  i < Num;  i++)
  {
    realnum Fi = F[i] ;
    realnum two_Fi = 2*Fi ;
    int i_mod = i - N_meas ; 
    ACmBS(i,i) = a[i][i] - c[i_mod][i_mod];
  	temp2 = ACmBS(i,i) * Fi;
    Aux_exp += temp2 * Fi;     
	if ( !no_integ_last || N_meas >= 2 )
	{
  	  PCpMS(i,i) = a[0][i]*a[0][i] + b[0][i]*b[0][i];
  	  temp1 = PCpMS(i,i) * Fi;
  	  Aux_Bessel += temp1 * Fi;
	  if ( improve_integ_now >= 2 )
	  {
		temp5 = F2sq_div_F1sq*a[1][i]*a[1][i] *Fi;
		temp6 = F2_div_F1*a[0][i]*a[1][i] *two_Fi;
		Aux_Bessel +=  temp5 * Fi;
		Aux_Bessel_cos2 +=  temp6 *Fi;
	  }
	}
    if ( do_der_Fph ) 
    {    
  	  dEdF[i] += temp2 ; 
  	  if ( !no_integ_last || N_meas >= 2 ) dRdF[i] += temp1 ;
	  if ( improve_integ_now >= 2 )
	  {
		if ( do_der_Fph == 2 ) {
		  secder_cos1(i,i) = 2*temp5*F_inv[i];
		  secder_cos2(i,i) = 2*temp6*F_inv[i];
        }
		dRdF[i] += temp5;
		dRdF_cos2[i] += temp6;
	  }
    }
	if ( no_integ_last )	// no_integ_last part begins
	{
      realnum cos_delta = cos_ph[i]*cos_ph[last] + sin_ph[i]*sin_ph[last] ;
	  realnum sin_delta = sin_ph[i]*cos_ph[last] - cos_ph[i]*sin_ph[last] ;
	  if (!cov->no_imag)
	  {	temp1 = F[last]*(a[last][i]*sin_delta+b[last][i]*cos_delta);
		temp2 = F[last]*(a[last][i]*cos_delta-b[last][i]*sin_delta);
	  }
	  else
	  {	temp1 = F[last]*a[last][i]*sin_delta;
		temp2 = F[last]*a[last][i]*cos_delta;
	  }
  	  Aux_exp += temp2 * two_Fi;
  	  if ( do_der_Fph ) 
  	  {    
    	dEdF[i] += temp2 ; 
        dEdph[i] -= temp1 * Fi  ;
	  }
	  if ( N_meas == 2 )
	  {
		if (!cov->no_imag)
		{
    	  temp1 = F[last] * ( ( a[0][last]*a[0][i] + b[0][last]*b[0][i] ) *cos_delta +
							  ( b[0][last]*a[0][i] - a[0][last]*b[0][i] ) *sin_delta );
		  temp2 = F[last] * ( ( a[0][last]*a[0][i] + b[0][last]*b[0][i] ) *sin_delta -
							  ( b[0][last]*a[0][i] - a[0][last]*b[0][i] ) *cos_delta );
		}
		else
		{
    	  temp1 = F[last] * a[0][last]*a[0][i] * cos_delta ;
    	  temp2 = F[last] * a[0][last]*a[0][i] * sin_delta ;
		}
  		Aux_Bessel += two_Fi * temp1 ;	
  		if ( do_der_Fph ) 
  		{    
    	  dRdF[i] += temp1 ;
    	  dRdph[i] -= Fi * temp2 ;
		}
	  }
	}		// end of no_integ_last part
    for (int  j = i+1;  j < Num;  j++)
    {
//      realnum sin_delta = sin( phase[j] - phase[i] ) ;
//      realnum cos_delta = cos( phase[j] - phase[i] ) ;
	  realnum Fi_Fj = Fi*F[j];
	  realnum two_Fi_Fj = two_Fi*F[j];
      realnum sin_delta = sin_ph[j]*cos_ph[i] - cos_ph[j]*sin_ph[i] ;
      realnum cos_delta = cos_ph[j]*cos_ph[i] + sin_ph[j]*sin_ph[i] ;
	  if ( do_der_D )
	  {
		sin_del(i,j) = sin_delta;
		cos_del(i,j) = cos_delta;
	  }
      int j_mod = j - N_meas ;
      temp1 = ( a[i][j] - c[i_mod][j_mod] ) ;
      if (!cov->no_imag) 
	  {
	    temp2 = ( b[i][j] - d[i_mod][j_mod] ) ;
        // only the lower triangle of ACmBS, ASpBC, PCpMS, PSmMC is stored
        ACmBS(i,j) = temp1*cos_delta - temp2*sin_delta;
        ASpBC(i,j) = temp1*sin_delta + temp2*cos_delta;
		if ( !no_integ_last || N_meas >= 2 )
		{
      	  temp1 = ( a[0][i]*a[0][j] + b[0][i]*b[0][j] ) ;
      	  temp2 = ( b[0][i]*a[0][j] - a[0][i]*b[0][j] ) ;
      	  PCpMS(i,j) = temp1*cos_delta + temp2*sin_delta ;
      	  PSmMC(i,j) = temp1*sin_delta - temp2*cos_delta ; 
		}
      }
	  else // if no_imag
	  {	
	    // only the lower triangle of ACmBS, ASpBC, PCpMS, PSmMC is stored
        ACmBS(i,j) = temp1*cos_delta;
        ASpBC(i,j) = temp1*sin_delta;
		if ( !no_integ_last || N_meas >= 2 )
		{
      	  temp1 = ( a[0][i]*a[0][j] ) ;
      	  PCpMS(i,j) = temp1*cos_delta ;
      	  PSmMC(i,j) = temp1*sin_delta ;      
		}
	  }
      Aux_exp += two_Fi_Fj * ACmBS(i,j) ;
	  if ( !no_integ_last || N_meas >= 2 ) 
		Aux_Bessel += two_Fi_Fj * PCpMS(i,j) ;
	  if ( improve_integ_now >= 2 )  
	  {
	    temp5 = F2sq_div_F1sq * a[1][i]*a[1][j] ;
	    temp6 = F2_div_F1 * (a[0][i]*a[1][j]+a[0][j]*a[1][i]) ;
	    temp7 = F2_div_F1 * (a[0][i]*a[1][j]-a[0][j]*a[1][i]) ;
		Aux_Bessel += two_Fi_Fj * temp5 * cos_delta;
		Aux_Bessel_cos2 += two_Fi_Fj * temp6 * cos_delta ;
		Aux_Bessel_sin2 -= two_Fi_Fj * temp7 * sin_delta ;
		if ( do_der_Fph ) {
		  dRdF[i] += F[j]*temp5*cos_delta;
		  dRdF[j] += Fi  *temp5*cos_delta;
		  dRdph[i] += Fi_Fj*temp5*sin_delta;
		  dRdph[j] -= Fi_Fj*temp5*sin_delta;
		  dRdF_cos2[i] += F[j]*temp6*cos_delta;
		  dRdF_cos2[j] += Fi  *temp6*cos_delta;
		  dRdF_sin2[i] -= F[j]*temp7*sin_delta;
		  dRdF_sin2[j] -= Fi  *temp7*sin_delta;
		  dRdph_cos2[i] -= Fi_Fj*temp6*sin_delta;
		  dRdph_cos2[j] += Fi_Fj*temp6*sin_delta;
		  dRdph_sin2[i] -= Fi_Fj*temp7*cos_delta;
		  dRdph_sin2[j] += Fi_Fj*temp7*cos_delta;		
		  if ( do_der_Fph == 2 ) {
		    secder_cos1(i,j) = 2*temp5*cos_delta;
		    secder_sin1(i,j) = 2*temp5*sin_delta;
		    secder_cos2(i,j) = 2*temp6*cos_delta;
		    secder_sin2(i,j) = 2*temp7*sin_delta;
		    secdermix_cos2(i,j) = 2*temp6*sin_delta;
		    secdermix_sin2(i,j) = -2*temp7*cos_delta;
		  }
		}
	  }
    }  
    if ( do_der_Fph ) 
	{
      for (int  j = N_meas;  j < Num;  j++)   
	  if ( j != i )
  	  {
    	    // here the upper triangle terms of ACmBS, ASpBC, PCpMS, PSmMC are used
      	if (i<j) 	{ dEdF[i] += F[j] * ACmBS(i,j) ;		dEdph[i] += Fi * F[j] * ASpBC(i,j) ;  }
		else		{ dEdF[i] += F[j] * ACmBS(j,i) ;		dEdph[i] -= Fi * F[j] * ASpBC(j,i) ;  }
		if ( !no_integ_last || N_meas >= 2 ) 
		{
      	  if (i<j) 	{ dRdF[i] += F[j] * PCpMS(i,j) ;	dRdph[i] += Fi * F[j] * PSmMC(i,j) ;  }
		  else 		{ dRdF[i] += F[j] * PCpMS(j,i) ;	dRdph[i] -= Fi * F[j] * PSmMC(j,i) ;  }
		}
// the following commented out code is already implemented above in the j=i+1..Num loop
/*		if ( improve_integ_now >= 2 ) 
		{
		  realnum sin_delta = sin_ph[j]*cos_ph[i] - cos_ph[j]*sin_ph[i] ;
		  realnum cos_delta = cos_ph[j]*cos_ph[i] + sin_ph[j]*sin_ph[i] ;
      	  dRdF[i] += F[j]*F2sq_div_F1sq*a[1][i]*a[1][j]*cos_delta ;
		  dRdF_cos2[i] += F[j]*F2_div_F1*(a[0][i]*a[1][j]+a[0][j]*a[1][i])*cos_delta;
		  dRdF_sin2[i] -= F[j]*F2_div_F1*(a[0][i]*a[1][j]-a[0][j]*a[1][i])*sin_delta;
		  dRdph_cos2[i] -= Fi*F[j]*F2_div_F1*(a[0][i]*a[1][j]+a[0][j]*a[1][i])*sin_delta;
		  dRdph_sin2[i] -= Fi*F[j]*F2_div_F1*(a[0][i]*a[1][j]-a[0][j]*a[1][i])*cos_delta;
      	  dRdph[i] +=  Fi*F[j]*F2sq_div_F1sq*a[1][i]*a[1][j]*sin_delta ;
		}
*/
	  }
	}
  }

  realnum Bef_term = 2;   
  for (int  i = 0;  i < Num-N_meas;  i++) 
    if ( eigenvalues1[Rice][i] > min_acc_eigen && eigenvalues2[Rice][i] > min_acc_eigen ) 
      Bef_term *= eigenvalues2[Rice][i]*inv_eigenval[i]; 
//	else {
//	  if ( eigenvalues1[Rice][i] > min_acc_eigen )   Bef_term *= inv_eigenval[i];
//	  if ( eigenvalues2[Rice][i] > min_acc_eigen )   Bef_term *= eigenvalues2[Rice][i];
//	}

  for (int i=0; i<N_meas; i++) 
  { 
    if ( F[i] > 0 ) 
  	  Bef_term *= F[i];  
	if ( eigenvalues1[Rice][Num-i-1] > min_acc_eigen ) 
	  Bef_term *= inv_eigenval[Num-i-1]; 
    if (i>0) Bef_term *= PI_INV; 
  }
  if ( cent && Rice && !no_integ_last ) 
  {
	Bef_term *= PI_INV*F_inv[0] ;  
	Bef_term = sqrt(Bef_term);  
  }
  if ( no_integ_last && N_meas == 1 ) Bef_term *= 0.5; 
  
// when N_meas == 1 (or not integrating out last phase) *****
  realnum sim, F0_isqrtR_sim, F0F0_iR_simder, twoR_inv;
  if ( N_meas == 1 || no_integ_last ) 
  {
	realnum BessI0_term , BessI1_term , Exp_term , cosh_term = 1;
    realnum Bess_arg = 0;
	if ( !no_integ_last || N_meas >= 2 )
	{
	  if ( Aux_Bessel < min_nonzero ) Aux_Bessel = min_nonzero ;
  	  Aux_Bessel = sqrt(Aux_Bessel);
  	  Bess_arg = F[0] * Aux_Bessel ;
	  if ( !cent || no_integ_last )  Bess_arg *= 2;
	  if ( Bess_arg < min_nonzero ) Bess_arg = min_nonzero ;
  	  if ( !cent || no_integ_last ) {
  		if (Tab_func) BessI0_term = tab->I0e( Bess_arg ) ;
		else BessI0_term = i0e( Bess_arg, Bess_prec1, Bess_prec2 );
	  }
	  else	if ( Bess_arg < max_hyparg )	cosh_term = cosh( Bess_arg ); 
	}
    Exp_term = - Aux_exp ;
	if ( cent && !no_integ_last ) Exp_term *= .5;	
    if ( ( do_der_Fph || do_der_D ) && ( !no_integ_last || N_meas >= 2 ) ) 
    {
	  sim = 1;
	  if ( !cent || no_integ_last )
	  {
    	if (Tab_func) 	BessI1_term = tab->I1e( Bess_arg ) ;
		else			BessI1_term = i1e( Bess_arg, Bess_prec1, Bess_prec2 ) ;
    	if (BessI0_term > min_nonzero) sim = BessI1_term/BessI0_term;
	  }	
	  else
	  {
		if ( Bess_arg < max_hyparg )  sim = tanh(Bess_arg);
	  }
      F0_isqrtR_sim = F[0]/Aux_Bessel*sim;
// three equivalent ways to compute FOM are here (the most simple is uncommented)	 
//      realnum integ = 0;
//      realnum FOM2 = 0, FOMsin = 0, FOMcos = 0;
//      realnum PHIB1=0,PHIB2=0;
//      for ( int m = 0;  m < Gauss_number;  m++ ) 
//      { 
//        realnum sin_delta = sin_ph[1]*cos_last[m] - cos_ph[1]*sin_last[m] ;
//        realnum cos_delta = cos_ph[1]*cos_last[m] + sin_ph[1]*sin_last[m] ;
//        realnum temp = exp( - 2*F[0]*F[1]*(a[0][1]*cos_delta-b[0][1]*sin_delta) ) * Gauss_struct[0].weight[m]; 
//        integ += temp;
//        FOMsin += temp * sin_last[m];
//        FOMcos += temp * cos_last[m];
//        FOM2 += temp * cos_delta;
//      }
//      FOM = sqrt( FOMsin*FOMsin + FOMcos*FOMcos ) /integ ;  
//      FOM2 = FOM2 / integ;  
//      cout << "FOM1 = " << FOM << " " << endl;  
//      cout << "FOM2 = " << FOM2 << endl;  
//      cout << "FOM from sim = " << BessI1_term/BessI0_term << endl;
//      PHIB1 = atan2(sin_ph[1],cos_ph[1]);
//      cout << "PHIB1 = " << PHIB1 << " " << endl;  
//	  PHIB2 = atan2(FOMsin, FOMcos);
//      cout << "PHIB2 = " << PHIB2 << endl;  
      FOM = sim;
      PHIB = phase[1];
      if ( do_ABCD ) {
    	realnum atanhfom = atanh(FOM);
    	if (Tab_func) {	
    	  HLA = atanhfom * tab->Cos(PHIB) ;
    	  HLB = atanhfom * tab->Sin_charged(PHIB) ;
    	}
    	else {			
    	  HLA = atanhfom * cos(PHIB) ;
    	  HLB = atanhfom * sin(PHIB) ;
    	}
      }
	}
	if ( do_der_Fph )
	{
      for (int i = N_meas;  i < Num;  i++) 
      {
	    dEdF[i] *= -2;		dEdph[i] *= -2;
        if ( no_integ_last && N_meas == 1 )
		{
		  der_F[i] =  - dEdF[i] ;
      	  der_ph[i] = - dEdph[i] ;
		}
		else
		{
		  dRdF[i] *= 2;		dRdph[i] *= 2;      
		  der_F[i] =  - F0_isqrtR_sim * dRdF[i] - dEdF[i] ;
      	  der_ph[i] = - F0_isqrtR_sim * dRdph[i] - dEdph[i] ;
		}
		if ( cent && !no_integ_last ) 
		{
		  der_F[i] *= .5;
		  der_ph[i] *= .5;
		}
      }
	}
    if ( ( do_der_Fph == 2 || do_der_D == 2 ) && ( !no_integ_last || N_meas >= 2 ) )
    {
	  F0F0_iR_simder = F[0]*F[0]/Aux_Bessel/Aux_Bessel;
	  if ( !cent || no_integ_last )
		F0F0_iR_simder *= ( 1 - sim*sim - sim/Bess_arg );
	  else 
	  {
		F0F0_iR_simder *= .5;
		if ( Bess_arg < max_hyparg ) 	F0F0_iR_simder /= ( cosh_term*cosh_term );
		else F0F0_iR_simder = 0;
	  }  
	  twoR_inv = 1./(2*Aux_Bessel*Aux_Bessel);
	}
    if ( do_der_Fph == 2 )
	{
	  realnum d2E, d2R;
      for (int i = N_meas;  i < Num;  i++) 
      {
        for (int j = N_meas;  j < Num;  j++) 
	  	{		// derivatives wrt ph,F
	  	  if ( i == j ) 
	  	  {  
	    	if ( F[i] < min_nonzero  &&  F[i] > -min_nonzero ) 
		  	  d2E = d2R = 0;
			else 
			{
			  d2E = dEdph[i]*F_inv[i];
			  if ( !no_integ_last || N_meas >= 2 )
				d2R = dRdph[i]*F_inv[i];
			}	
	  	  }
	  	  else 
	  	  {	    	// upper or lower triangle terms 
	    	if (i>j) d2E = 2*F[i]*ASpBC(j,i); else d2E = -F[i]*ASpBC(i,j);
			if ( !no_integ_last || N_meas >= 2 )	
	    	  if (i>j) d2R = -2*F[i]*PSmMC(j,i); else d2R = 2*F[i]*PSmMC(i,j);
	  	  }  
		  if ( no_integ_last && N_meas == 1 )
	  		der2_phF[i][j] = - d2E;
		  else
	  		der2_phF[i][j] = 
	    	  - F0_isqrtR_sim * (-dRdF[j]*dRdph[i]*twoR_inv+d2R) -
	    	  d2E - F0F0_iR_simder * dRdph[i]*dRdF[j];		     

	  	  if ( i <= j ) 
	  	  {		// derivatives wrt ph,ph
	    	if ( i == j )
	    	{
	      	  d2E = - F[i]*( dEdF[i] + 2*F[i]*ACmBS(i,i) );
			  if ( !no_integ_last || N_meas >= 2 )	
	      		d2R = - F[i]*( dRdF[i] - 2*F[i]*PCpMS(i,i) );
	    	}
	    	else
	    	{
	      	  d2E = -2*F[i]*F[j]*ACmBS(i,j);
			  if ( !no_integ_last || N_meas >= 2 )	
	      		d2R = 2*F[i]*F[j]*PCpMS(i,j);
	    	}
			if ( no_integ_last && N_meas == 1 )
	    	  der2_ph2[i][j] = - d2E;
			else
	    	  der2_ph2[i][j] = 
	      		- F0_isqrtR_sim * (-dRdph[i]*dRdph[j]*twoR_inv+d2R) -
				d2E - F0F0_iR_simder * dRdph[i]*dRdph[j] ;	
				// derivatives wrt F,F  (no need to distinguish diagonal terms here)
	    	d2E = -2*ACmBS(i,j);
			if ( no_integ_last && N_meas == 1 )
	    	  der2_F2[i][j] = - d2E;
			else
			{
	    	  d2R = 2*PCpMS(i,j);	
	  		  der2_F2[i][j] = 
				- F0_isqrtR_sim * (-dRdF[i]*dRdF[j]*twoR_inv+d2R) -
				d2E - F0F0_iR_simder * dRdF[i]*dRdF[j] ;	
			}
          	if ( cent && !no_integ_last )   
			{
			  der2_F2[i][j] *= .5;
              der2_phF[i][j] *= .5;
              der2_ph2[i][j] *= .5;
			}
	  		der2_ph2[j][i] = der2_ph2[i][j];
	    	der2_F2[j][i] = der2_F2[i][j];
		  }    
		}      
	  }
    }

	if (cent)	// this is only for 1. function !
	{
  	  fval = - Exp_term - log( cosh_term * Bef_term ) ;
	  if ( Bess_arg >= max_hyparg ) {  fval += - Bess_arg + log(2.) ;  }
	}
	else
	{
  	  if ( !no_integ_last || N_meas >= 2 ) 
		fval = - Exp_term - log( BessI0_term * Bef_term ) - Bess_arg;	 
	  else
		fval = - Exp_term - log( Bef_term );
	}
  }


// when N_meas == 2 or 3 *****
  realnum* sqrtR = work.space3 + Gauss_number;		// not used if N_meas == 3 - sqrtRtot used instead
  realnum* R3 = sqrtR;								// only used if N_meas == 3
  realnum *for_derD___F0_isqrtR_weight_E2_BessI1;	// not used if N_meas == 3 - *_siras used instead
  realnum *for_derD___weight_E2_BessI0;				// not used if N_meas == 3 - *_siras used instead
  realnum** for_derD___F0_isqrtR_weight_E2_BessI1_siras;		// used only when N_meas == 3
  realnum** for_derD___weight_E2_BessI0_siras;					// used only when N_meas == 3
  realnum *Aux_Bessel_cos3, *Aux_Bessel_sin3;					// only used if N_meas == 3
  int n_top=1, n_bottom=-1; 	// the boundaries in integration over n
  int m_top=1, m_bottom=0; 		// the boundaries in integration over m
  int m_rej_num=0, m_rej_bottom=Gauss_number; 		// defines the "hole" in integration range over m
  realnum integ = 0, integ_inv = 0;
  if ( N_meas >= 2 && !no_integ_last ) 
  {	
	m_top = Gauss_number_now;
	if (N_meas>2)
	{
	  for ( int k=0;  k<Gauss_number_now_2;  k++ )
  		for ( int l=0;  l<Gauss_number_now;  l++ )
		  cos_last12_del[k][l] = cos_last2[k]*cos_last[l] + sin_last2[k]*sin_last[l] ;
	}
	if ( do_der_D )
	{
	  if ( N_meas == 2 )
	  {
		for_derD___F0_isqrtR_weight_E2_BessI1 = work.space3 + 2*Gauss_number;
		for_derD___weight_E2_BessI0 = work.space3 + 3*Gauss_number ; 
	  }
	  else // N_meas == 3
	  {
		for_derD___F0_isqrtR_weight_E2_BessI1_siras = work.space6 + 2*Gauss_number_2;
		for_derD___weight_E2_BessI0_siras = work.space6 + 3*Gauss_number_2 ; 
	  }
	}
	if ( improve_integ_now >= 2 && N_meas>2 ) 
	{
	  Aux_Bessel_cos3 = work.space3 + 2*Gauss_number;		
	  Aux_Bessel_sin3 = work.space3 + 3*Gauss_number;		
	}
    
    realnum Bess_arg = 0;  
	realnum FOMsin = 0, FOMcos = 0;
//	realnum forFOMmax = -100000.;	FOM = 0;
	if ( improve_integ_now < 2 )
  	  Aux_Bessel += F[last]*F[last]*
                  ( a[0][last]*a[0][last] + b[0][last]*b[0][last] ) ;
	if ( N_meas == 3 ) 	
	  if ( improve_integ_now >= 2 )   
	  {
		Aux_Bessel += F[2]*F[2]*( a[0][2]*a[0][2] + F2sq_div_F1sq*a[1][2]*a[1][2] ) ;
		Aux_Bessel_cos2 += 2*F2_div_F1*F[2]*F[2]*a[0][2]*a[1][2] ;
	  }
	  else   Aux_Bessel += F[last2]*F[last2]*a[0][last2]*a[0][last2] ;
    realnum Aux_Bess_int = Aux_Bessel;
    realnum *dRdF_con = work.space1 + 16*Num,       *dE2dF = work.space1 + 17*Num;
    realnum *dRdph_con = work.space1 + 18*Num,      *dE2dph = work.space1 + 19*Num;
    realnum *temp1 = work.space1 + 20*Num,          *temp2 = work.space1 + 21*Num;
    realnum *alastj_Flast = work.space1 + 22*Num,   *blastj_Flast = work.space1 + 23*Num;
	realnum *alast2j_Flast2 = blastj_Flast;	// only used if N_meas == 3
    for ( int j = N_meas;  j < Num;  j++ )
    {
      if ( do_der_Fph )
      { 
        dEdF[j] *= -2;		dEdph[j] *= -2;
		dRdF[j] *= 2;		dRdph[j] *= 2;      	  
		if ( improve_integ_now >= 2 )  {
		  dRdF_cos2[j] *= 2;		dRdph_cos2[j] *= 2; 
		  dRdF_sin2[j] *= 2;		dRdph_sin2[j] *= 2; 
		}
        dRdF_con[j] = dRdF[j] ;
        dE2dF[j] = 0 ;
        dRdph_con[j] = dRdph[j] ;
        dE2dph[j] = 0 ;
      }
      alastj_Flast[j] =  2 * a[last][j] * F[last] ;
      if ( N_meas == 3 )  alast2j_Flast2[j] =  2 * a[last2][j] * F[last2] ;
      if (!cov->no_imag)
      {
        blastj_Flast[j] =  2 * b[last][j] * F[last] ;
        temp1[j] = 2 * F[last] * ( a[0][last]*a[0][j] + b[0][last]*b[0][j] ) ;
        temp2[j] = 2 * F[last] * ( b[0][last]*a[0][j] - a[0][last]*b[0][j] ) ;       
      }
      else
	  {
        temp1[j] = 2 * F[last] * a[0][last]*a[0][j] ;
		if ( N_meas == 3 )  
		  if ( improve_integ_now >= 2 ) temp2[j] = 2 * F2sq_div_F1sq*F[2] * a[1][j]*a[1][2] ;
		  else							temp2[j] = 2 * F[last2] * a[0][last2]*a[0][j] ;							
        else temp2[j] = 0 ;       
		if ( N_meas == 3 && improve_integ_now >= 2 )
		  temp3[j] = 2*F2_div_F1*F[2]*(a[0][2]*a[1][j]+a[0][j]*a[1][2]);
		  temp4[j] = 2*F2_div_F1*F[2]*(a[0][2]*a[1][j]-a[0][j]*a[1][2]);
	  }
    }

	realnum* E2arg = work.space3;			// not used if N_meas == 3 - E23argtot used instead
	realnum* E3arg = E2arg;					// only used if N_meas == 3
	realnum twoF1 = 2*F[0];
realnum integ_m[1000];
	realnum E2argmax = -1e100;      // in order to avoid abnormally great numbers in integration
// index m is used for the integration over the last phase, n over the last-1 phase (if applicable)
    for ( int m = m_bottom;  m < m_top;  m++ ) 
    {
      realnum Aux_int = 0;
      Aux_Bess_int = Aux_Bessel;
  	  if ( do_ABCD && N_meas == 3 )  integ_m[m] = 0.;
	  if ( improve_integ_now >= 2 && N_meas == 3 )   Aux_Bessel_cos3[m] = Aux_Bessel_sin3[m] = 0;
	  if ( improve_integ_now >= 2 && N_meas == 2 )
	  {
		realnum cos_2last=cos_last[m]*cos_last[m]-sin_last[m]*sin_last[m];
		realnum sin_2last=2*sin_last[m]*cos_last[m];
		Aux_Bess_int += cos_2last*Aux_Bessel_cos2 + sin_2last*Aux_Bessel_sin2;
		Aux_int += 2*F[0]*F[1]*a[0][1]*cos_2last;
	  }
	  else
	  {
    	for ( int j = N_meas;  j < Num;  j++ )
  		{
      	  sinlast_del[m][j] = sin_ph[j]*cos_last[m] - cos_ph[j]*sin_last[m] ;
      	  coslast_del[m][j] = cos_ph[j]*cos_last[m] + sin_ph[j]*sin_last[m] ;
      	  if (!cov->no_imag)
      	  {
      		Aux_Bess_int += F[j] * ( temp1[j]*coslast_del[m][j] + temp2[j]*sinlast_del[m][j] ) ;
	  		Aux_int += F[j] * ( alastj_Flast[j]*coslast_del[m][j] - blastj_Flast[j]*sinlast_del[m][j] ) ;
      	  }
		  else
		  {
      		Aux_Bess_int += F[j] * ( temp1[j]*coslast_del[m][j] ) ;
	  		Aux_int += F[j] * ( alastj_Flast[j]*coslast_del[m][j] ) ;	
		  }
		  if ( improve_integ_now >= 2 )
		  {
			Aux_Bess_int += 2 * F2sq_div_F1sq * F[last]*F[j] * a[1][last]*a[1][j] * coslast_del[m][j] ;
			Aux_Bessel_cos3[m] += 2 * F2_div_F1 * F[last]*F[j] * (a[0][last]*a[1][j]+a[1][last]*a[0][j]) * coslast_del[m][j] ;
			Aux_Bessel_sin3[m] -= 2 * F2_div_F1 * F[last]*F[j] * (a[0][last]*a[1][j]-a[1][last]*a[0][j]) * sinlast_del[m][j] ;
		  }
  		} 
	  }
	  if ( N_meas == 2 ) 
	  {
  		if ( Aux_Bess_int < min_nonzero ) { Aux_Bess_int = min_nonzero ;  }
		else Aux_Bess_int = sqrt(Aux_Bess_int);
  		Bess_arg = twoF1 * Aux_Bess_int ;
		if ( Bess_arg < min_nonzero ) { Bess_arg = min_nonzero ;  }
		sqrtR[m] = Aux_Bess_int ;
		E2arg[m] = fabs(Bess_arg) - Aux_int ;
		if ( E2argmax < E2arg[m] ) E2argmax = E2arg[m];
	  }
	  else    // if N_meas == 3
	  {
		E3arg[m] = - Aux_int ;
		R3[m] = Aux_Bess_int ;
	  }
	}

	realnum** sqrtRtot = work.space6 ;					// used only when N_meas == 3
	realnum** E23argtot = work.space6 + Gauss_number_2 ;	// used only when N_meas == 3
    if ( N_meas == 3 )				// only when N_meas==3 
	{
	  realnum twoF2F3a23 = 2*F[1]*F[2]*a[1][2];
	  realnum twoF2F3a12a13 = 2*F[1]*F[2]*a[0][1]*a[0][2];
	  realnum cos_2last2, sin_2last2;	// only used if  improve_integ >= 2
	  realnum Aux_int;
// loops over n - if improve_integ>=2 then only while the contributions to integral are significant
	  bool end_n_loop_top=0, end_n_loop_bottom=0;
	  int Gauss_num_2_half=Gauss_number_now_2/2, n;
	  if ( improve_integ_now >= 2 )	{ n_top = n_bottom = Gauss_num_2_half; }
	  else   						{ n_top = Gauss_number_now_2; n_bottom = -1; n = 0; }
	  realnum contrib_top=0, contrib_top_prev=1, contrib_bottom=0, contrib_bottom_prev=1, contrib_max=0, contrib_n=0; // only used if improve_integ >= 2
	  while ( !end_n_loop_top || !end_n_loop_bottom ) 		// looping over n - getting the terms dependent on second phase
  	  {
    	realnum Aux_int_E2 = 0;
    	realnum Aux_Bess_int_R2 = 0;
		if ( improve_integ_now >= 2 ) 
		{
		  if ( contrib_bottom>=contrib_top ) n=n_bottom; else n=n_top;
		  contrib_n = 0;
		  cos_2last2=cos_last2[n]*cos_last2[n]-sin_last2[n]*sin_last2[n];
		  sin_2last2=2*sin_last2[n]*cos_last2[n];
		  Aux_Bess_int_R2 = cos_2last2*Aux_Bessel_cos2 + sin_2last2*Aux_Bessel_sin2;
		  Aux_int_E2 = 2*F[0]*F[1]*a[0][1]*cos_2last2;
		}
		else
    	  for ( int j = N_meas;  j < Num;  j++ )
  		  {
      		sinlast2_del[n][j] = sin_ph[j]*cos_last2[n] - cos_ph[j]*sin_last2[n] ;
      		coslast2_del[n][j] = cos_ph[j]*cos_last2[n] + sin_ph[j]*sin_last2[n] ;
      		Aux_Bess_int_R2 += F[j] * ( temp2[j]*coslast2_del[n][j] ) ;
	  		Aux_int_E2 += F[j] * ( alast2j_Flast2[j]*coslast2_del[n][j] ) ;	
  		  } 
  		for ( int m = m_bottom;  m < m_top;  (m==m_rej_bottom-1) ? m+=m_rej_num+1 : m++ )  // getting the terms dependent on both second and third phase
  		{ 
		  if ( improve_integ_now >= 2 )
		  {
			Aux_Bess_int = Aux_Bess_int_R2 + R3[m] + cos_2last2*Aux_Bessel_cos3[m] + sin_2last2*Aux_Bessel_sin3[m];
			Aux_int = - Aux_int_E2 + E3arg[m] ;
		  }
		  else
		  {
			Aux_Bess_int = Aux_Bess_int_R2 + R3[m] + twoF2F3a12a13*cos_last12_del[n][m] ;
			Aux_int = - Aux_int_E2 + E3arg[m] - twoF2F3a23*cos_last12_del[n][m] ;
		  }
  		  if ( Aux_Bess_int < min_nonzero ) { Aux_Bess_int = min_nonzero ;  }
		  else Aux_Bess_int = sqrt(Aux_Bess_int);
  		  Bess_arg = twoF1 * Aux_Bess_int ;
		  if ( Bess_arg < min_nonzero ) { Bess_arg = min_nonzero ;  }
		  sqrtRtot[n][m] = Aux_Bess_int ;
		  E23argtot[n][m] = Bess_arg + Aux_int ;
		  if ( E2argmax < E23argtot[n][m] ) E2argmax = E23argtot[n][m];
		  if ( improve_integ_now >= 2 ) 
			if (Tab_func) contrib_n += tab->ExpM(-E23argtot[n][m]+E23argtot[Gauss_num_2_half][0]);
			else			contrib_n += exp(E23argtot[n][m]-E23argtot[Gauss_num_2_half][0]);
//if (improve_integ >= 2 ) cout << n << " "<< m << " : " << E23argtot[n][m] << endl;
		}
		if ( improve_integ_now >= 2 ) 
		{
		  if ( contrib_n > contrib_max )   contrib_max=contrib_n;
		  if ( n == n_bottom )   { 
			contrib_bottom_prev=contrib_bottom; contrib_bottom=contrib_n; n_bottom--; 
			if ( contrib_bottom*contrib_bottom/contrib_bottom_prev<0.0006*contrib_max && contrib_bottom<contrib_bottom_prev ) end_n_loop_bottom=1;
			if ( n_bottom < 0 )  { end_n_loop_bottom = 1; contrib_bottom = 0; }
		  }
		  if ( n == n_top )   { 
			contrib_top_prev=contrib_top; contrib_top=contrib_n; n_top++; 
			if ( contrib_top*contrib_top/contrib_top_prev<0.0006*contrib_max && contrib_top<contrib_top_prev ) end_n_loop_top=1;
			if ( n_top >= Gauss_number_now_2 )  { end_n_loop_top = 1; contrib_top = 0; }
		  }
  
		  if (n==Gauss_num_2_half&&improve_integ_now>3) 
		  { 
			int new_m_bottom=Gauss_number_now, new_m_top=0, m_rej_top=Gauss_number_now;
			for ( int m = m_bottom;  m < m_top;  m++ ) 
			  if (E2argmax-E23argtot[n][m]<12.5) {
				if (m<new_m_bottom) new_m_bottom=m;
				new_m_top=m; 
			  }
			  else {
				if (new_m_bottom<m && m<m_rej_bottom) { m_rej_bottom=m_rej_top=m; }
				if (m==m_rej_top+1) m_rej_top++;
			  } 
			m_bottom=new_m_bottom; m_top=new_m_top+1;
			if (m_rej_top<m_top) m_rej_num=m_rej_top-m_rej_bottom+1;
			if (m_rej_num<0)  m_rej_num=0; 
		  }
		}
		else 
		  if (++n==n_top)   end_n_loop_bottom = end_n_loop_top = 1;
		if ( debug_info == 1 ) {
		  realnum testmin=1e30;
		  if (n==Gauss_num_2_half) for ( int m = 0;  m < Gauss_number_now;  m++ )  
			if (testmin>E23argtot[n][m]) testmin=E23argtot[n][m];
		  if (n==Gauss_num_2_half) cout << "diff: " << E2argmax - testmin << endl;
		}
//	if (n==Gauss_num_2_half)  for ( int m = m_bottom;  m < m_top;  (m==m_rej_bottom-1) ? m+=m_rej_num+1 : m++ )
//	  cout << m << " : " << E23argtot[n][m] << endl;
	  }
	  if ( debug_info == 1 ) cout << "Number of points used for integration: " << n_top-n_bottom -1 << " (from " << n_bottom+1 <<  " to " << n_top-1 << " ; sampling "<< Gauss_number_now_2<<")" <<endl;
	  if ( debug_info == 1 ) cout << "num_m_points:" << m_top-m_bottom-m_rej_num << " ; sampling "<< Gauss_number_now << endl;
	  int act_num_points = (n_top-n_bottom-1)*(m_top-m_bottom-m_rej_num);
	  if ( max_num_points_allowed < act_num_points && improve_integ >= 1 ) {  
		fval = -act_num_points; return fval;
	  }
	}	// end of N_meas == 3 if


	realnum  E2, weight_E2_BessI0=0, weight_E2_BessI1=0, F0_isqrtR_weight_E2_BessI1=0;
	realnum weight_log_E2_BessI0=0;
	realnum exp_arg;
	for ( int n = n_bottom+1;  n < n_top;  n++ )
	{
	realnum contrib_n = 0;
	realnum cos_2last2, sin_2last2;
	if ( N_meas > 2 && improve_integ_now >= 2 ) {
	  cos_2last2=cos_last2[n]*cos_last2[n]-sin_last2[n]*sin_last2[n];
	  sin_2last2=2*sin_last2[n]*cos_last2[n];
    }
    for ( int m = m_bottom;  m < m_top;  (m==m_rej_bottom-1) ? m+=m_rej_num+1 : m++ ) 
    {
	  if ( N_meas == 2 )
	  {
		exp_arg = E2argmax - E2arg[m];
		Bess_arg = twoF1 * sqrtR[m] ;
	  }
	  else 		// if N_meas == 3
	  {
		exp_arg = E2argmax - E23argtot[n][m] ;
		Bess_arg = twoF1 * sqrtRtot[n][m] ;
	  }
	  if (exp_arg <= 50. || do_ABCD) 
	  {
	    if ( Tab_func && exp_arg >= -50. && exp_arg<= 50. ) 
	  	  E2 = tab->ExpM_nocheck( exp_arg );
	    else
	  	  E2 = exp( -exp_arg ); 
	    if ( Tab_func && Bess_arg <= 200000. )
      	  weight_E2_BessI0 = Gauss_weight[m] * E2 * tab->I0e_nocheck( Bess_arg ) ;
	    else
  	  	  weight_E2_BessI0 = Gauss_weight[m] * E2 * i0e( Bess_arg, Bess_prec1, Bess_prec2 ) ;
  	    if ( do_ABCD ) {
  	  	  if ( N_meas == 2 ) {
  	  	    if ( Tab_func && Bess_arg <= 200000. )
  	  		  weight_log_E2_BessI0 = Gauss_weight[m] * log( E2 * tab->I0e_nocheck( Bess_arg ) ) ;
  	  	    else
  	  		  weight_log_E2_BessI0 = Gauss_weight[m] * log( E2 * i0e( Bess_arg, Bess_prec1, Bess_prec2 ) ) ;
	      } 
	      else { // N_meas == 3
  	  	    if ( Tab_func && Bess_arg <= 200000. )
  	  		  integ_m[m] += Gauss_weight_2[n] * E2 * tab->I0e_nocheck( Bess_arg ) ;
  	  	    else
  	  		  integ_m[m] += Gauss_weight_2[n] * E2 * i0e( Bess_arg, Bess_prec1, Bess_prec2 ) ;
		  }
  	    }
	    contrib_n += E2;
	    if ( N_meas == 3 )
	  	  weight_E2_BessI0 *= Gauss_weight_2[n] ;
  	    if ( do_der_Fph || do_der_D )
	    {
	  	  if ( Tab_func && Bess_arg <= 200000. )
  	  		weight_E2_BessI1 = Gauss_weight[m] * E2 * tab->I1e_nocheck_charged( Bess_arg ) ;
	  	  else
  	  		weight_E2_BessI1 = Gauss_weight[m] * E2 * i1e( Bess_arg, Bess_prec1, Bess_prec2 ) ;
	  	  if ( N_meas == 2 )
      		F0_isqrtR_weight_E2_BessI1 = F[0] / sqrtR[m] * weight_E2_BessI1 ;
	  	  else		// if N_meas == 3
	  		F0_isqrtR_weight_E2_BessI1 = F[0] / sqrtRtot[n][m] * Gauss_weight_2[n] * weight_E2_BessI1 ;
	    }
  	    integ += weight_E2_BessI0;
	    if ( do_der_Fph || do_der_D ) {
  	  	  FOMsin += weight_E2_BessI0 * sin_last[m];
  	  	  FOMcos += weight_E2_BessI0 * cos_last[m];
	    }
	    if ( do_ABCD && N_meas == 2 ) {
  	  	  HLA += weight_log_E2_BessI0 * cos_last[m] ;
  	  	  HLB += weight_log_E2_BessI0 * sin_last[m] ;
  	  	  HLC += weight_log_E2_BessI0 * (cos_last[m]*cos_last[m]-sin_last[m]*sin_last[m]) ;
  	  	  HLD += 2*weight_log_E2_BessI0 * sin_last[m] * cos_last[m] ;
  	    }
//	    FOM += weight_E2_BessI0;
//	    if ( weight_E2_BessI0 > forFOMmax ) forFOMmax = weight_E2_BessI0;
      
  	    if ( do_der_D ) 
  	    {
	  	  if ( N_meas == 2 )
	  	  {
        	for_derD___F0_isqrtR_weight_E2_BessI1[m] = F0_isqrtR_weight_E2_BessI1 ;
	    	for_derD___weight_E2_BessI0[m] = weight_E2_BessI0 ;
	  	  }
	  	  else	// N_meas == 3
	  	  {
        	for_derD___F0_isqrtR_weight_E2_BessI1_siras[n][m] = F0_isqrtR_weight_E2_BessI1 ;
	    	for_derD___weight_E2_BessI0_siras[n][m] = weight_E2_BessI0 ;
	  	  }
  	    }
      
  	    if ( do_der_Fph )
	    {
	  	  for (int  i = N_meas;  i < Num;  i++) 
      	  {
        	realnum sin_delta = sinlast_del[m][i];
        	realnum cos_delta = coslast_del[m][i];
        	if (!cov->no_imag)
        	{
        	  dRdF[i] = dRdF_con[i] + temp1[i]*cos_delta + temp2[i]*sin_delta ;
	    	  dRdph[i] = dRdph_con[i] - F[i]*( temp1[i]*sin_delta - temp2[i]*cos_delta ) ;
        	  dE2dF[i] = - ( alastj_Flast[i]*cos_delta - blastj_Flast[i]*sin_delta ) ;	
        	  dE2dph[i] = F[i]*( alastj_Flast[i]*sin_delta + blastj_Flast[i]*cos_delta ) ;
        	}
	  		else 
        	{
	  		  if ( improve_integ_now >= 2 && N_meas == 2 ) 
	  		  {
	  			realnum cos_2last=cos_last[m]*cos_last[m]-sin_last[m]*sin_last[m];
	  			realnum sin_2last=2*sin_last[m]*cos_last[m];
	  			dRdF[i] = dRdF_con[i] + cos_2last* dRdF_cos2[i] + sin_2last* dRdF_sin2[i] ;
	  			dE2dF[i] = dE2dph[i] = 0;
	  			dRdph[i] = dRdph_con[i] - cos_2last* dRdph_cos2[i] - sin_2last* dRdph_sin2[i] ;
	  		  }
	  		  else
	  		  {
        		dRdF[i] = dRdF_con[i] + temp1[i]*cos_delta ;
        		dE2dF[i] = - alastj_Flast[i]*cos_delta ;
	    		dRdph[i] = dRdph_con[i] - F[i]*temp1[i]*sin_delta ;
        		dE2dph[i] = F[i]*alastj_Flast[i]*sin_delta ;
	  		  }
        	}
	  		if ( N_meas == 3 ) 
	  		{
	  		  if ( improve_integ_now >= 2 ) {
	  			dRdF[i] += temp2[i]*cos_delta + cos_2last2*( dRdF_cos2[i] + temp3[i]*cos_delta )
	  										  + sin_2last2*( dRdF_sin2[i] - temp4[i]*sin_delta ) ;
	  			dRdph[i] -= temp2[i]*F[i]*sin_delta + cos_2last2*( dRdph_cos2[i] + temp3[i]*F[i]*sin_delta )
	  												+ sin_2last2*( dRdph_sin2[i] + temp4[i]*F[i]*cos_delta ) ;
	  		  }
	  		  else {
        		dRdF[i] += temp2[i]*coslast2_del[n][i] ;
	    		dRdph[i] += - F[i]*temp2[i]*sinlast2_del[n][i] ;
        		dE2dF[i] += - alast2j_Flast2[i]*coslast2_del[n][i] ;
        		dE2dph[i] += F[i]*alast2j_Flast2[i]*sinlast2_del[n][i] ;
	  		  }
	  		}
        	der_F[i] +=  F0_isqrtR_weight_E2_BessI1*dRdF[i] +
  	    	                      weight_E2_BessI0*dE2dF[i] ;
        	der_ph[i] += F0_isqrtR_weight_E2_BessI1*dRdph[i] +
  	  		                      weight_E2_BessI0*dE2dph[i] ;
	  	  }
        }
	    if ( do_der_Fph == 2 )
	    {
	  	  realnum d2R, d2E2;
	  	  realnum R_inv; 
	  	  if (N_meas==2)  R_inv = 1./(sqrtR[m]*sqrtR[m]) ;
	  	  else            R_inv = 1./(sqrtRtot[n][m]*sqrtRtot[n][m]) ;
	      for (int  i = N_meas;  i < Num;  i++) 
      	  {
	    	for (int  j = N_meas;  j < Num;  j++) 
	  		{			// for derivatives wrt ph,F
	    	  if ( i == j ) 
	    	  {  
	      		if ( F[i] < min_nonzero  &&  F[i] > -min_nonzero ) 
	  		  	  d2R = d2E2 = 0.;
	  			else 
	  			{
	  			  d2E2 = dE2dph[i]*F_inv[i];
	      		  d2R = dRdph[i]*F_inv[i];
	  			}	
	    	  }
	    	  else 
	    	  {	    	// upper or lower triangle terms 
	  			d2E2 = 0.;
	      		if (i>j) d2R = -2*F[i]*PSmMC(j,i); else d2R = 2*F[i]*PSmMC(i,j);
	    	  }  
	    	  if ( N_meas>2 && improve_integ_now>=2 ) {
	    		if (i<j) d2R += F[i]*( secder_sin1(i,j) + secdermix_cos2(i,j)*cos_2last2 - secdermix_sin2(i,j)*sin_2last2 );
	    		if (i>j) d2R -= F[i]*( secder_sin1(j,i) + secdermix_cos2(j,i)*cos_2last2 - secdermix_sin2(j,i)*sin_2last2 );
	    	  }
	    	  der2_phF[i][j] += 
	  			F0_isqrtR_weight_E2_BessI1*( dE2dph[i]*dRdF[j] + dE2dF[j]*dRdph[i] - dRdph[i]*dRdF[j]*R_inv + d2R ) +
	  			      	  weight_E2_BessI0*( dE2dph[i]*dE2dF[j] + d2E2 + dRdph[i]*dRdF[j]*F[0]*F[0]*R_inv ) ;
	  		  
	    	  if ( i <= j )
	    	  {		// for derivatives wrt ph,ph
	      		if ( i == j )
	      		{
	        	  d2R = - F[i]*( dRdF[i] - 2*F[i]*PCpMS(i,i) );
	  			  d2E2 = - F[i]*dE2dF[i];
	      		}
	      		else
	      		{
	        	  d2R = 2*F[i]*F[j]*PCpMS(i,j);
	  			  d2E2 = 0.;
	      		}
	        	if ( N_meas>2 && improve_integ_now>=2 ) {
                  if (i==j) d2R += F[i]*F[i]*( secder_cos1(i,j) + secder_cos2(i,j)*cos_2last2 );
                  else      d2R += F[i]*F[j]*( secder_cos1(i,j) + secder_cos2(i,j)*cos_2last2 - secder_sin2(i,j)*sin_2last2 );
                }
	      		der2_ph2[i][j] +=
	  			  F0_isqrtR_weight_E2_BessI1*( dE2dph[i]*dRdph[j] + dE2dph[j]*dRdph[i] - dRdph[i]*dRdph[j]*R_inv + d2R ) +
	  			      		weight_E2_BessI0*( dE2dph[i]*dE2dph[j] + d2E2 + dRdph[i]*dRdph[j]*F[0]*F[0]*R_inv ) ;
      
	      		// for derivatives wrt F,F  (no need to distinguish diagonal terms here)		  
	        	d2R = 2*PCpMS(i,j);
	        	if ( N_meas>2 && improve_integ_now>=2 ) {
                  if (i==j) d2R += secder_cos1(i,j) + secder_cos2(i,j)*cos_2last2;
                  else      d2R += secder_cos1(i,j) + secder_cos2(i,j)*cos_2last2 - secder_sin2(i,j)*sin_2last2;
//	        	  if (i==j) d2R += 2*F2sq_div_F1sq*a[1][i]*a[1][i] + 4*F2_div_F1*cos_2last2*a[0][i]*a[1][i];
//                  else      d2R += 2*F2sq_div_F1sq*a[1][i]*a[1][j]*cos(phase[j]-phase[i]) + 2*F2_div_F1*(a[0][j]*a[1][i]*cos(phase[j]-phase[i])+a[0][i]*a[1][j]*cos(phase[i]-phase[j]))*cos_2last2 - 2*F2_div_F1*(a[0][j]*a[1][i]*sin(phase[j]-phase[i])+a[0][i]*a[1][j]*sin(phase[i]-phase[j]))*sin_2last2 ;
                }
	  			d2E2 = 0;
	  			der2_F2[i][j] += 
	  			  F0_isqrtR_weight_E2_BessI1*( dE2dF[i]*dRdF[j] + dE2dF[j]*dRdF[i] - dRdF[i]*dRdF[j]*R_inv + d2R ) +
	  			          	weight_E2_BessI0*( dE2dF[i]*dE2dF[j] + d2E2 + dRdF[i]*dRdF[j]*F[0]*F[0]*R_inv ) ;
	  		  }
	  		}
	  	  }	
	    }
	  }
	  else  // if E2 is effectively 0.
	  {
  	    if ( do_der_D ) 
	  	  if ( N_meas == 2 )
        	for_derD___F0_isqrtR_weight_E2_BessI1[m] = for_derD___weight_E2_BessI0[m] = 0. ;
	  	  else	// N_meas == 3
        	for_derD___F0_isqrtR_weight_E2_BessI1_siras[n][m] = for_derD___weight_E2_BessI0_siras[n][m] = 0. ;
	  }
	}		// end of m loop
/////if (improve_integ >= 2 ) cout << n << "  real: " << contrib_n << endl;
    }  		// end of n loop

	integ_inv = 1./integ;
	if ( N_meas == 2 )
  	  fval = - log( Bef_term*integ*PI ) + Aux_exp - E2argmax;
	else {    // N_meas == 3
  	  fval = - log( Bef_term*integ*PI*PI ) + Aux_exp - E2argmax;
  	  if (do_ABCD) {
  		for ( int m = m_bottom;  m < m_top;  (m==m_rej_bottom-1) ? m+=m_rej_num+1 : m++ ) 
  		{
  		  realnum log_int_m_gaus = log(integ_m[m]) * Gauss_weight[m];
  	  	  HLA += log_int_m_gaus * cos_last[m] ;
  	  	  HLB += log_int_m_gaus * sin_last[m] ;
  	  	  HLC += log_int_m_gaus * (cos_last[m]*cos_last[m]-sin_last[m]*sin_last[m]) ;
  	  	  HLD += 2*log_int_m_gaus * sin_last[m] * cos_last[m] ;  		  
  	    }
  	  }
  	}
//cout << Bef_term << " " << integ*PI << " " <<  Aux_exp << " " << E2argmax << endl;
	FOM = sqrt( FOMsin*FOMsin + FOMcos*FOMcos ) *integ_inv;
//	FOM = 1 - FOM/(forFOMmax*Gauss_number_now);//*FOM/(forFOMmax*Gauss_number_now);
	PHIB = atan2(FOMsin, FOMcos);

    if ( do_der_Fph )
	{
      if ( do_der_Fph == 2 ) 
	  {		// 1. derivatives of I are now in der_*;  2. deriv. of I in der2_* variables
        realnum d2E;
        for (int i = N_meas;  i < Num;  i++) 
        {
          for (int j = N_meas;  j < Num;  j++) 
	      {		// derivatives wrt ph,F
	  		if ( i == j ) 
	  		{  
	    	  if ( F[i] < min_nonzero  &&  F[i] > -min_nonzero ) 
		  		d2E = 0;
			  else 
				d2E = dEdph[i]*F_inv[i];
	  		}
	  		else 
	  		{	    	// upper or lower triangle terms 
	    	  if (i>j) d2E = 2*F[i]*ASpBC(j,i); else d2E = -2*F[i]*ASpBC(i,j);
	  		}  
	  		der2_phF[i][j] = - d2E - ( der2_phF[i][j]*integ - der_ph[i]*der_F[j] )*integ_inv*integ_inv ;
		  
	  		if ( i <= j ) 
	  		{		// derivatives wrt ph,ph
	    	  if ( i == j )
	      		d2E = - F[i]*( dEdF[i] + 2*F[i]*ACmBS(i,i) );
	    	  else
	      		d2E = -2*F[i]*F[j]*ACmBS(i,j);
	  		  der2_ph2[i][j] = - d2E - ( der2_ph2[i][j]*integ - der_ph[i]*der_ph[j] )*integ_inv*integ_inv ;
	    	  der2_ph2[j][i] = der2_ph2[i][j];
		  
	    		// derivatives wrt F,F  (no need to distinguish diagonal terms here)
	  		  d2E = -2*ACmBS(i,j);
	  		  der2_F2[i][j] = - d2E - ( der2_F2[i][j]*integ - der_F[i]*der_F[j] )*integ_inv*integ_inv ;
	    	  der2_F2[j][i] = der2_F2[i][j];
			}
		  }
		}
	  }
      for (int i=N_meas; i<Num; i++) 
      {
        der_F[i] = - dEdF[i] - der_F[i]*integ_inv ;
        der_ph[i] = - dEdph[i] - der_ph[i]*integ_inv ;
      }
	}	  
  
  }	    // end of N_meas >=2 if


// this is only about storage
#undef PCpMS
#undef PSmMC
#undef ACmBS
#undef ASpBC
  



// computation of derivatives wrt D[k][l]

  if ( do_der_D ) 
  {

	realnum aux1, aux2, aux3, aux4;
// compute derivatives of [i][j] element of inverse to complex matrix mr+i*mi wrt D[k,l] (resp. D[pm,k,l])
#define comp_dinvD_r( mr, mi, i, j, pm, k, l )           						 	 - (    \
    cov->part[pm].deriv_D_r[k][l] * ( mr[i][k]*mr[l][j] - mi[i][k]*mi[l][j] + mr[i][l]*mr[k][j] - mi[i][l]*mi[k][j] ) +        \
    cov->part[pm].deriv_D_i[k][l] * (-mi[i][k]*mr[l][j] - mr[i][k]*mi[l][j] + mi[i][l]*mr[k][j] + mr[i][l]*mi[k][j] )     ) 
#define comp_dinvD_i( mr, mi, i, j, pm, k, l )           						 	 - (    \
    cov->part[pm].deriv_D_r[k][l] * ( mi[i][k]*mr[l][j] + mr[i][k]*mi[l][j] + mi[i][l]*mr[k][j] + mr[i][l]*mi[k][j] ) +        \
    cov->part[pm].deriv_D_i[k][l] * ( mr[i][k]*mr[l][j] - mi[i][k]*mi[l][j] - mr[i][l]*mr[k][j] + mi[i][l]*mi[k][j] )     ) 
// only to save some time for diagonal terms    
#define comp_dinvdiagD_r( mr, mi, i, pm, k, l )                     -2 * (   \
    cov->part[pm].deriv_D_r[k][l] * ( mr[i][k]*mr[l][i] - mi[i][k]*mi[l][i] ) +         \
    cov->part[pm].deriv_D_i[k][l] * (-mi[i][k]*mr[l][i] - mr[i][k]*mi[l][i] )       ) 
#define comp_dinvdiagD_i( mr, mi, i, pm, k, l )                	0 

// and now for case of real terms only
#define comp_dinvD_r_no_imag( mr, mi, i, j, pm, k, l )           					 - (    \
    cov->part[pm].deriv_D_r[k][l] * ( mr[i][k]*mr[l][j] + mr[i][l]*mr[k][j] )          )
// only to save some time for diagonal terms    
#define comp_dinvdiagD_r_no_imag( mr, mi, i, pm, k, l )                     -2 * (   \
    cov->part[pm].deriv_D_r[k][l] * mr[i][k]*mr[l][i]                            )

// if k = l (i.e. der. wrt D[k,k] ) - this MUST be used if k=l
#define comp_dinvD_r_kk( mr, mi, i, j, pm, k )           						 	    \
	cov->no_imag ? 																		\
	- cov->part[pm].deriv_D_r[k][k] *   mr[i][k]*mr[k][j] :								\
    - cov->part[pm].deriv_D_r[k][k] * ( mr[i][k]*mr[k][j] - mi[i][k]*mi[k][j] ) 
#define comp_dinvD_i_kk( mr, mi, i, j, pm, k )           						 	 0  \


#define comp_aux4inv( mr, mi, dmr, dmi, i, j, p, q, k, l )			\
	if ( p == q )  																					{	\
	  if (cov->no_imag) 	{																			\
		aux1 = dmr[i][p]*mr[q][j] + mr[i][p]*dmr[q][j];													\
		aux3 = 0;			}																			\
	  else 					{																			\
		aux1 = dmr[i][p]*mr[q][j] - dmi[i][p]*mi[q][j] + mr[i][p]*dmr[q][j] - mi[i][p]*dmi[q][j] ;		\
		aux3 = dmi[i][p]*mr[q][j] + dmr[i][p]*mi[q][j] + mi[i][p]*dmr[q][j] + mr[i][p]*dmi[q][j] ;		\
							}																		}	\
	else if ( i == p || j == q )															{	\
	  aux1 = dmr[i][p]*mr[q][j]	+ mr[i][p]*dmr[q][j] ;	\
	  if ( i == p && j == q ) 									{	\
		aux3 = 0;	\
		aux2 = 2*( dmr[i][q]*mr[p][j] - dmi[i][q]*mi[p][j]	) ; 	\
		aux4 = 2*( dmi[i][q]*mr[p][j] + dmr[i][q]*mi[p][j]	) ;	}	\
	  else																	{	\
		aux2 = dmr[i][q]*mr[p][j] - dmi[i][q]*mi[p][j] + mr[i][q]*dmr[p][j] - mi[i][q]*dmi[p][j] ;	\
		aux4 = dmi[i][q]*mr[p][j] + dmr[i][q]*mi[p][j] + mi[i][q]*dmr[p][j] + mr[i][q]*dmi[p][j] ; 	\
		if ( i == p )	aux3 = dmr[i][p]*mi[q][j] + mr[i][p]*dmi[q][j] ;		\
		else			aux3 = dmi[i][p]*mr[q][j] + mi[i][p]*dmr[q][j] ; 	}				}	\
	else																									{	\
	  aux1 = dmr[i][p]*mr[q][j] - dmi[i][p]*mi[q][j] + mr[i][p]*dmr[q][j] - mi[i][p]*dmi[q][j] ;	\
	  aux3 = dmi[i][p]*mr[q][j] + dmr[i][p]*mi[q][j] + mi[i][p]*dmr[q][j] + mr[i][p]*dmi[q][j] ; 	\
	  if ( i == q )											{	\
		aux2 = dmr[i][q]*mr[p][j] + mr[i][q]*dmr[p][j] ;		\
		aux4 = dmr[i][q]*mi[p][j] + mr[i][q]*dmi[p][j] ; 	}	\
	  else if ( p == j ) 									{	\
		aux2 = dmr[i][q]*mr[p][j] + mr[i][q]*dmr[p][j] ;		\
		aux4 = dmi[i][q]*mr[p][j] + mi[i][q]*dmr[p][j] ; 	}	\
	  else																							{	\
		aux2 = dmr[i][q]*mr[p][j] - dmi[i][q]*mi[p][j] + mr[i][q]*dmr[p][j] - mi[i][q]*dmi[p][j] ;	\
		aux4 = dmi[i][q]*mr[p][j] + dmr[i][q]*mi[p][j] + mi[i][q]*dmr[p][j] + mr[i][q]*dmi[p][j] ; 	}		}

// and now for case of real terms only
#define comp_aux4inv_no_imag( mr, mi, dmr, dmi, i, j, p, q, k, l )			\
	aux1 = dmr[i][p]*mr[q][j]	+ mr[i][p]*dmr[q][j] ;						\
	if (p!=q) aux2 = dmr[i][q]*mr[p][j] + mr[i][q]*dmr[p][j] ;		


#define add_d2invD_r( res, mr, mi, i, j, pm2, p, q, pm, k, l )  	  	res -= 	 	(	\
    cov->partpart[pm][pm2].deriv_D_r[p][q] * ( mr[i][p]*mr[q][j] - mi[i][p]*mi[q][j] + mr[i][q]*mr[p][j] - mi[i][q]*mi[p][j] ) +        	\
    cov->partpart[pm][pm2].deriv_D_i[p][q] * (-mi[i][p]*mr[q][j] - mr[i][p]*mi[q][j] + mi[i][q]*mr[p][j] + mr[i][q]*mi[p][j] ) +    	 	\
    cov->part[pm2].deriv_D_r[p][q] * ( aux1 + aux2 )  +   \
    cov->part[pm2].deriv_D_i[p][q] * ( aux4 - aux3 )		);


#define add_d2invD_i( res, mr, mi, i, j, pm2, p, q, pm, k, l )          res	-=	 	(  	\
    cov->partpart[pm][pm2].deriv_D_r[p][q] * ( mi[i][p]*mr[q][j] + mr[i][p]*mi[q][j] + mi[i][q]*mr[p][j] + mr[i][q]*mi[p][j] ) +        	\
    cov->partpart[pm][pm2].deriv_D_i[p][q] * ( mr[i][p]*mr[q][j] - mi[i][p]*mi[q][j] - mr[i][q]*mr[p][j] + mi[i][q]*mi[p][j] ) +    		\
    cov->part[pm2].deriv_D_r[p][q] * ( aux3 + aux4 ) +   	\
    cov->part[pm2].deriv_D_i[p][q] * ( aux1 - aux2 )    	);


// only to save some time for diagonal terms ( has real impact only for small matrix dimensions )	; comp_aux not required for this macro
#define add_d2invdiagD_r( res, mr, mi, dmr, dmi, i, pm2, p, q, pm, k, l )           	\
	if ( i == p )		\
	  res	-= 	2*( cov->partpart[pm][pm2].deriv_D_r[p][q] * mr[i][p]*mr[q][i] -        	\
  					cov->partpart[pm][pm2].deriv_D_i[p][q] * mr[i][p]*mi[q][i] +   	 	\
  					cov->part[pm2].deriv_D_r[p][q] * ( dmr[i][p]*mr[q][i] + mr[i][p]*dmr[q][i] ) +        \
  					cov->part[pm2].deriv_D_i[p][q] * (-dmr[i][p]*mi[q][i] - mr[i][p]*dmi[q][i] )		);  \
	else if ( i == q )	\
	  res	-= 	2*( cov->partpart[pm][pm2].deriv_D_r[p][q] * mr[i][p]*mr[q][i] -        	\
  					cov->partpart[pm][pm2].deriv_D_i[p][q] * mi[i][p]*mr[q][i] +   	 	\
  					cov->part[pm2].deriv_D_r[p][q] * ( dmr[i][p]*mr[q][i] + mr[i][p]*dmr[q][i] ) +        \
  					cov->part[pm2].deriv_D_i[p][q] * (-dmi[i][p]*mr[q][i] - mi[i][p]*dmr[q][i] )		) ; \
	else	\
	  res	-= 	2*( cov->partpart[pm][pm2].deriv_D_r[p][q] * ( mr[i][p]*mr[q][i] - mi[i][p]*mi[q][i] ) +        	\
  					cov->partpart[pm][pm2].deriv_D_i[p][q] * (-mi[i][p]*mr[q][i] - mr[i][p]*mi[q][i] ) +   	 	\
  					cov->part[pm2].deriv_D_r[p][q] * ( dmr[i][p]*mr[q][i] - dmi[i][p]*mi[q][i] +        	\
					    				 mr[i][p]*dmr[q][i] - mi[i][p]*dmi[q][i] ) +        \
  					cov->part[pm2].deriv_D_i[p][q] * (-dmi[i][p]*mr[q][i] - dmr[i][p]*mi[q][i] +     		\
										-mi[i][p]*dmr[q][i] - mr[i][p]*dmi[q][i] )		) ;
#define add_d2invdiagD_i( res, mr, mi, i, p, q, k, l )           	;

// and now for case of real terms only
#define add_d2invD_r_no_imag( res, mr, mi, i, j, pm2, p, q, pm, k, l )  	  	res -= 	 	(	\
    cov->partpart[pm][pm2].deriv_D_r[p][q] * ( mr[i][p]*mr[q][j] + mr[i][q]*mr[p][j] ) +        \
    cov->part[pm2].deriv_D_r[p][q] * ( aux1 + aux2 ) 		);
// only to save some time for diagonal terms ( has real impact only for small matrix dimensions ) ; comp_aux not required for this macro
#define add_d2invdiagD_r_no_imag( res, mr, dmr, i, pm2, p, q, pm, k, l )           	\
	res	-= 	2*( cov->partpart[pm][pm2].deriv_D_r[p][q] * mr[i][p]*mr[q][i] +        	\
  				cov->part[pm2].deriv_D_r[p][q] * ( dmr[i][p]*mr[q][i] + mr[i][p]*dmr[q][i] ) );  \


// if p = q (i.e. der. wrt D[k,l],D[p,p] ) - this MUST be used if p=q !
#define add_d2invD_r_pp( res, mr, mi, i, j, pm2, p, pm, k, l )  	  {	res -= 	 	(	\
    cov->partpart[pm][pm2].deriv_D_r[p][p] * mr[i][p]*mr[p][j]  + 						\
    cov->part[pm2].deriv_D_r[p][p] * aux1   );											\
  	if (!cov->no_imag)													res +=		(	\
    cov->partpart[pm][pm2].deriv_D_r[p][p] *  mi[i][p]*mi[p][j]   );  }
#define add_d2invD_i_pp( res, mr, mi, i, j, pm2, p, pm, k, l )          res	-=	 	(  	\
    cov->partpart[pm][pm2].deriv_D_r[p][p] * ( mi[i][p]*mr[p][j] + mr[i][p]*mi[p][j] ) +      \
    cov->part[pm2].deriv_D_r[p][p] * aux3   );   	


// store 1. raw of matrix of derivatives of inverse cov matrix terms wrt D[k][l] in derD_raw1	  
  realnum *derD_raw1_a = work.space1 + 10*Num, 	*derD_raw1_b = work.space1 + 11*Num; 
  realnum *tmp4dElast = work.space1 + 12*Num,	*tmp4dElast2 = work.space1 + 13*Num;
  realnum *tmp4dRlast = work.space1 + 14*Num,	*tmp4dRlast2 = work.space1 + 15*Num;
  realnum *derD_raw2_a = derD_raw1_b ;		// only used if improve_integ >= 2
  realnum ****derDpq_raw1_a, ****derDpq_raw1_b;
  realnum ****der2D_raw1_a, ****der2D_raw1_b;
  realnum ****tmp4d2E, ****tmp4d2R;
  if ( do_der_D == 2 )	
  { 
	derDpq_raw1_a = work.space4;				derDpq_raw1_b = work.space4 + cov->N_part;
	der2D_raw1_a = work.space4 + 2*cov->N_part;	der2D_raw1_b = work.space4 + 3*cov->N_part;
	tmp4d2E = work.space4 + 4*cov->N_part;		tmp4d2R = work.space4 + 5*cov->N_part;
  }
  realnum dEdD, dE2dD, dRdD_con;
  realnum dRdD_cos2, dRdD_sin2;  	// only used if improve_integ >= 2
  realnum **derD_a = storage1, 					**derD_b = storage2;
  realnum **derD_c = work.space2 + 4*Num,		**derD_d = work.space2 + 5*Num;
  int *kD = cov->workspace,						*lD = cov->workspace + cov->maxtoNum;
  realnum ***d2 = work.space8,					***d_dD = work.space8 + cov->N_part;

#define dIdD(pm,i,j) d_dD[pm][i][j]
#define dRdD(pm,i,j) d_dD[pm][j][i]
#define d2E(pm,i,j) d2[pm][i][j]
#define d2R(pm,i,j) d2[pm][j][i]
#define int_dE2dD(pm,m,i,j) work.space5[pm][m][i][j]
#define int_dRdD(pm,m,i,j) work.space5[pm][m][j][i]
#define tmp4d2E1(pm,i,p,q) tmp4d2E[pm][i][p][q]
#define tmp4d2E2(pm,i,p,q) tmp4d2E[pm][i][q][p]
#define tmp4d2R1(pm,i,p,q) tmp4d2R[pm][i][p][q]
#define tmp4d2R2(pm,i,p,q) tmp4d2R[pm][i][q][p]

#define to(i,j) to_array[2*i+j]

  for ( int i = Num-1;  i >= N_meas;  i-- )		// moving inverse cov2 matrix terms for convenience 
  {
	for ( int j = Num-1;  j >= N_meas;  j-- )
	{
	  c[i][j] = c[i-N_meas][j-N_meas]; 
	  d[i][j] = d[i-N_meas][j-N_meas]; 
	}
  }


  realnum *R = sqrtR;
  if ( do_der_D == 2 && N_meas == 2 && !no_integ_last )
  {
	for ( int m = m_bottom;  m < m_top;  (m==m_rej_bottom-1) ? m+=m_rej_num+1 : m++ )
	{
	  R[m] *= R[m];
	}
  }


  for (int  k = 0;  k < cov->DNum;  k++)
  {
    for (int  l = k+1;  l < cov->DNum;  l++)
    {  
	  for ( int pm = 0;  pm < cov->N_part;  pm++ )			// now looping over k,l,pm to get derivatives wrt D[pm,k,l]
	  {

	  int NDact = cov->Dass[pm][k][l].toNum;
	  for ( int nD = 0;  nD < NDact;  nD++ )
	  {
		kD[nD] = cov->Dass[pm][k][l].to(nD,0);
		lD[nD] = cov->Dass[pm][k][l].to(nD,1);
	  }
	  
      der_D[pm][k][l] = 0;
      // supposing that 2. derivative of cov[i][j] wrt D[i][j]'s is 0 !!!
//	  realnum deriv2_cov_r_klpq = 0, deriv2_cov_i_klpq = 0;
  
// precompute all derivatives of cov. matrix terms wrt D params
	  for (int  i = 0;  i < Num;  i++)
	  {
		derD_a[i][i] = derD_b[i][i] = 0;
		for ( int nD = 0;  nD < NDact;  nD++ ) 
		{
		  if ( kD[nD] == lD[nD] ) 
			derD_a[i][i] += comp_dinvD_r_kk( a, b, i, i, pm, kD[nD] );
		  else {
			if (!cov->no_imag) {
			  derD_a[i][i] += comp_dinvdiagD_r( a, b, i, pm, kD[nD], lD[nD] );
			  derD_b[i][i] += comp_dinvdiagD_i( a, b, i, pm, kD[nD], lD[nD] );
			}
			else   
			  derD_a[i][i] += comp_dinvdiagD_r_no_imag( a, b, i, pm, kD[nD], lD[nD] );
		  }
		}
		if ( i >= N_meas ) 
		{
		  derD_c[i][i] = derD_d[i][i] = 0;
		  for ( int nD = 0;  nD < NDact;  nD++ ) 
		  {
			if ( kD[nD] >= N_meas ) 
			{ 
			  if ( kD[nD] == lD[nD] ) 
				derD_c[i][i] += comp_dinvD_r_kk( c, d, i, i, pm, kD[nD] );
			  else {
				if (!cov->no_imag) {
				  derD_c[i][i] += comp_dinvdiagD_r( c, d, i, pm, kD[nD], lD[nD] );
				  derD_d[i][i] += comp_dinvdiagD_i( c, d, i, pm, kD[nD], lD[nD] );
				}
				else
				  derD_c[i][i] += comp_dinvdiagD_r_no_imag( c, d, i, pm, kD[nD], lD[nD] );
			  }
			}
		  }
		}
		if ( i == 0 )
		{	
		  derD_raw1_a[i] = derD_a[i][i];
		  derD_raw1_b[i] = derD_b[i][i];
		  if ( do_der_D == 2 )  {
			derDpq_raw1_a[pm][k][l][i] = derD_raw1_a[i] ;
			derDpq_raw1_b[pm][k][l][i] = derD_raw1_b[i] ;
		  }
		}
		if ( improve_integ_now >= 2 && i == 1 )  derD_raw2_a[i] = derD_a[i][i];
		for (int  j = i+1;  j < Num;  j++)
		{
		  derD_a[i][j] = derD_b[i][j] = 0;
		  for ( int nD = 0;  nD < NDact;  nD++ )
		  {
			if ( kD[nD] == lD[nD] ) 
			  derD_a[i][j] += comp_dinvD_r_kk( a, b, i, j, pm, kD[nD] );
			else {
			  if (!cov->no_imag) {
				derD_a[i][j] += comp_dinvD_r( a, b, i, j, pm, kD[nD], lD[nD] );
				derD_b[i][j] += comp_dinvD_i( a, b, i, j, pm, kD[nD], lD[nD] );
			  }
			  else
				derD_a[i][j] += comp_dinvD_r_no_imag( a, b, i, j, pm, kD[nD], lD[nD] );
			}
		  }
		  derD_a[j][i] = derD_a[i][j];
		  derD_b[j][i] = - derD_b[i][j];
		  if ( i >= N_meas )
		  {
			derD_c[i][j] = derD_d[i][j] = 0;
			for ( int nD = 0;  nD < NDact;  nD++ )
			{
			  if ( kD[nD] >= N_meas )
			  {
				if ( kD[nD] == lD[nD] ) 
				  derD_c[i][j] += comp_dinvD_r_kk( c, d, i, j, pm, kD[nD] );
				else {
				  if (!cov->no_imag) {
					derD_c[i][j] += comp_dinvD_r( c, d, i, j, pm, kD[nD], lD[nD] );
					derD_d[i][j] += comp_dinvD_i( c, d, i, j, pm, kD[nD], lD[nD] );
				  }
				  else
					derD_c[i][j] += comp_dinvD_r_no_imag( c, d, i, j, pm, kD[nD], lD[nD] );
				}
			  }
			}
			derD_c[j][i] = derD_c[i][j];
			derD_d[j][i] = - derD_d[i][j];
		  }
		  if ( i == 0 ) 
		  {
			derD_raw1_a[j] = derD_a[i][j]; 	
			derD_raw1_b[j] = derD_b[i][j];
			if ( do_der_D == 2 )  {
			  derDpq_raw1_a[pm][k][l][j] = derD_raw1_a[j] ;
			  derDpq_raw1_b[pm][k][l][j] = derD_raw1_b[j] ;
			}
		  }
		  if ( improve_integ_now >= 2 && i == 1 )  derD_raw2_a[j] = derD_a[i][j]; 	
		}
	  }
// end of precalculation of derivatives of cov. matrix terms wrt D params

      dRdD(pm,k,l) = 0;
      dEdD = - F[0]*F[0]*derD_raw1_a[0];
	  if ( no_integ_last && N_meas >= 2 )
	  {
		dEdD -= F[last]*F[last]*derD_a[last][last];
      	dRdD(pm,k,l) += 2*F[last]*F[last]*( a[0][last]*derD_raw1_a[last] + b[0][last]*derD_raw1_b[last] ) ;
	  }		// end of no_integ_last part
	  if ( improve_integ_now >= 2 )  dRdD_cos2 = dRdD_sin2 = 0;
      for (int  i = N_meas;  i < Num;  i++)
      {
        realnum temp1 = ( a[0][i]*derD_raw1_a[i] + b[0][i]*derD_raw1_b[i] ) *2*F[i]*F[i] ;
        realnum temp2 = ( derD_a[i][i] - derD_c[i][i] ) *F[i]*F[i] ;
        dEdD -= temp2 ;
		if ( !no_integ_last || N_meas >= 2 )
		  dRdD(pm,k,l) += temp1 ;
		if ( improve_integ_now >= 2 )
		{
		  dRdD(pm,k,l) +=  F2sq_div_F1sq*a[1][i]*derD_raw2_a[i] *2*F[i]*F[i] ;
		  dRdD_cos2 += F2_div_F1*(a[0][i]*derD_raw2_a[i]+a[1][i]*derD_raw1_a[i]) *2*F[i]*F[i] ;
		}
		if ( no_integ_last )
		{
		  realnum derD_a_last_i = 0, derD_b_last_i = 0;
    	  realnum cos_delta = cos_ph[i]*cos_ph[last] + sin_ph[i]*sin_ph[last] ;
		  realnum sin_delta = sin_ph[i]*cos_ph[last] - cos_ph[i]*sin_ph[last] ;
		  if ( N_meas >= 2 )
		  {
			if (!cov->no_imag)
        	  dRdD(pm,k,l) += 2*F[last]*F[i] * (
				( derD_raw1_a[i]*a[0][1] + a[0][i]*derD_raw1_a[1] + derD_raw1_b[i]*b[0][1] + b[0][i]*derD_raw1_b[1] )*cos_delta +
        		( derD_raw1_a[i]*b[0][1] + a[0][i]*derD_raw1_b[1] - derD_raw1_b[i]*a[0][1] - b[0][i]*derD_raw1_a[1] )*sin_delta );
			else
        	  dRdD(pm,k,l) += 2*F[last]*F[i] * 
				( derD_raw1_a[i]*a[0][1] + a[0][i]*derD_raw1_a[1] )*cos_delta ;
			if ( do_der_D == 2 )	
			{
			  derD_a_last_i = derD_a[last][i];
			  derD_b_last_i = derD_b[last][i];
			}
		  }	
		  else	// if N_meas == 1
		  {	
			derD_a_last_i = derD_raw1_a[i];
		  	derD_b_last_i = derD_raw1_b[i];
		  }
		  if (!cov->no_imag)
			dEdD -= 2*F[i]*F[last]*(derD_a_last_i*cos_delta-derD_b_last_i*sin_delta);
		  else
			dEdD -= 2*F[i]*F[last]*derD_a_last_i*cos_delta;
		}	// end of no_integ_last part
        for (int  j = i+1;  j < Num;  j++)
        {
		  temp1 = ( derD_a[i][j] - derD_c[i][j] ) ;
    	  if (!cov->no_imag)
    	  {
        	temp2 = ( derD_b[i][j] - derD_d[i][j] ) ;
        	dEdD -= 2*F[i]*F[j]* ( temp1*cos_del(i,j) - temp2*sin_del(i,j) ) ;
			if ( !no_integ_last || N_meas >= 2 )	
			{
        	  temp1 = ( derD_raw1_a[i]*a[0][j] + a[0][i]*derD_raw1_a[j] + derD_raw1_b[i]*b[0][j] + b[0][i]*derD_raw1_b[j] ) ;
        	  temp2 = ( derD_raw1_b[i]*a[0][j] + b[0][i]*derD_raw1_a[j] - derD_raw1_a[i]*b[0][j] - a[0][i]*derD_raw1_b[j] ) ;
        	  dRdD(pm,k,l) += 2*F[i]*F[j]* ( temp1*cos_del(i,j) + temp2*sin_del(i,j) ) ;
			}
		  }
		  else
    	  {
        	dEdD -= 2*F[i]*F[j]*temp1*cos_del(i,j) ;
			if ( !no_integ_last || N_meas >= 2 )	
			{
        	  temp1 = derD_raw1_a[i]*a[0][j] + a[0][i]*derD_raw1_a[j] ;
        	  dRdD(pm,k,l) += 2*F[i]*F[j]*temp1*cos_del(i,j) ;
			}
			if ( improve_integ_now >= 2 )  
			{
			  dRdD(pm,k,l) += 2*F[i]*F[j]*F2sq_div_F1sq*(derD_raw2_a[i]*a[1][j] + a[1][i]*derD_raw2_a[j])*cos_del(i,j) ;
			  dRdD_cos2 += 2*F[i]*F[j]*F2_div_F1*(  derD_raw1_a[i]*a[1][j] + a[0][i]*derD_raw2_a[j] +
													derD_raw2_a[i]*a[0][j] + a[1][i]*derD_raw1_a[j] ) * cos_del(i,j) ;
			  dRdD_sin2 -= 2*F[i]*F[j]*F2_div_F1*(  derD_raw1_a[i]*a[1][j] + a[0][i]*derD_raw2_a[j] -
													derD_raw2_a[i]*a[0][j] - a[1][i]*derD_raw1_a[j] ) * sin_del(i,j) ;
			}
		  }
        }
      }		// end of i = N_meas .. Num loop

	  if ( do_der_D == 2 )
	  {
		for (int  p = 0;  p <= k;  p++) 
		{
		  for (int  q = p+1;  (q < cov->DNum && p < k) || ( q <= l && p == k );  q++)
		  {
			for (int pm2 = 0;  ( pm2 <= pm && p == k && q == l ) || ( pm2 < cov->N_part && (p < k || q < l) );  pm2++) 
			{
			  for (int  i = 0;  i < Num;  i++)
			  {
				der2D_raw1_a[pm2][p][q][i] = der2D_raw1_b[pm2][p][q][i] = 0;
			  }
			  for ( int nD = 0;  nD < cov->Dass[pm2][p][q].toNum;  nD++ )
			  {	
  				int pD = cov->Dass[pm2][p][q].to(nD,0);
  				int qD = cov->Dass[pm2][p][q].to(nD,1);
				for (int  i = 0;  i < Num;  i++)
				{
    			  if ( pD == qD ) {
					comp_aux4inv( a, b, derD_a, derD_b, 0, i, pD, qD, k, l );
					add_d2invD_r_pp( der2D_raw1_a[pm2][p][q][i], a, b, 0, i, pm2, pD, pm, k, l );
					if (!cov->no_imag) add_d2invD_i_pp( der2D_raw1_b[pm2][p][q][i], a, b, 0, i, pm2, pD, pm, k, l );
				  }
				  else 
    				if (!cov->no_imag) {
					  comp_aux4inv( a, b, derD_a, derD_b, 0, i, pD, qD, k, l );
					  add_d2invD_r( der2D_raw1_a[pm2][p][q][i], a, b, 0, i, pm2, pD, qD, pm, k, l );
					  add_d2invD_i( der2D_raw1_b[pm2][p][q][i], a, b, 0, i, pm2, pD, qD, pm, k, l );
					}
					else 
					{
					  comp_aux4inv_no_imag( a, b, derD_a, derD_b, 0, i, pD, qD, k, l );
					  add_d2invD_r_no_imag( der2D_raw1_a[pm2][p][q][i], a, b, 0, i, pm2, pD, qD, pm, k, l ); 
					}
				}
			  }
			  d2E(pm2,p,q) = - F[0]*F[0]*der2D_raw1_a[pm2][p][q][0] ;
			  d2R(pm2,p,q) = 0;
			  if ( no_integ_last && N_meas >= 2 )
			  {
				realnum der2D_a_last_last = 0.;
				for ( int nD = 0;  nD < cov->Dass[pm2][p][q].toNum;  nD++ )
				{	
  				  int pD = cov->Dass[pm2][p][q].to(nD,0);
  				  int qD = cov->Dass[pm2][p][q].to(nD,1);
    			  if ( pD == qD ) {
					comp_aux4inv( a, b, derD_a, derD_b, last, last, pD, qD, k, l );
					add_d2invD_r_pp( der2D_a_last_last, a, b, last, last, pm2, pD, pm, k, l );
				  }
				  else {
    				if (!cov->no_imag) 
					  { add_d2invdiagD_r( der2D_a_last_last, a, b, derD_a, derD_b, last, pm2, pD, qD, pm, k, l ) }
    				else 
					  { add_d2invdiagD_r_no_imag( der2D_a_last_last, a, derD_a, last, pm2, pD, qD, pm, k, l ) }
				  }
				}
				d2R(pm2,p,q) += ( derDpq_raw1_a[pm2][p][q][last]*derD_raw1_a[last] + a[0][last]*der2D_raw1_a[pm2][p][q][last] +
						  		  derDpq_raw1_b[pm2][p][q][last]*derD_raw1_b[last] + b[0][last]*der2D_raw1_b[pm2][p][q][last] ) *2*F[last]*F[last] ;
				d2E(pm2,p,q) -= F[last]*F[last]*der2D_a_last_last ;
			  }
    		  for (int  i = N_meas;  i < Num;  i++)
    		  {
				realnum der2D_aii = 0., der2D_cii = 0.;
				for ( int nD = 0;  nD < cov->Dass[pm2][p][q].toNum;  nD++ )
				{	
  				  int pD = cov->Dass[pm2][p][q].to(nD,0);
  				  int qD = cov->Dass[pm2][p][q].to(nD,1);
    			  if ( pD == qD ) {
					comp_aux4inv( a, b, derD_a, derD_b, i, i, pD, qD, k, l );
					add_d2invD_r_pp( der2D_aii, a, b, i, i, pm2, pD, pm, k, l );
					if (pD>=N_meas)  {
					  comp_aux4inv( c, d, derD_c, derD_d, i, i, pD, qD, k, l );
					  add_d2invD_r_pp( der2D_cii, c, d, i, i, pm2, pD, pm, k, l );
					}
				  }
				  else 
    				if (!cov->no_imag) {
					  add_d2invdiagD_r( der2D_aii, a, b, derD_a, derD_b, i, pm2, pD, qD, pm, k, l );
					  if (pD>=N_meas)  add_d2invdiagD_r( der2D_cii, c, d, derD_c, derD_d, i, pm2, pD, qD, pm, k, l );
					}
    				else
    				{
					  add_d2invdiagD_r_no_imag( der2D_aii, a, derD_a, i, pm2, pD, qD, pm, k, l );
					  if (pD>=N_meas)  add_d2invdiagD_r_no_imag( der2D_cii, c, derD_c, i, pm2, pD, qD, pm, k, l );
					}
				}
				if ( !no_integ_last || N_meas >= 2 )	
				  d2R(pm2,p,q) += ( derDpq_raw1_a[pm2][p][q][i]*derD_raw1_a[i] + a[0][i]*der2D_raw1_a[pm2][p][q][i] +
						  			derDpq_raw1_b[pm2][p][q][i]*derD_raw1_b[i] + b[0][i]*der2D_raw1_b[pm2][p][q][i] ) *2*F[i]*F[i] ;
				d2E(pm2,p,q) -= ( der2D_aii - der2D_cii ) *F[i]*F[i] ;
				if ( no_integ_last )
				{
				  realnum der2D_a_last_i = 0, der2D_b_last_i = 0;
    			  realnum cos_delta = cos_ph[i]*cos_ph[last] + sin_ph[i]*sin_ph[last] ;
				  realnum sin_delta = sin_ph[i]*cos_ph[last] - cos_ph[i]*sin_ph[last] ;
				  if ( N_meas >= 2 )
				  {
					if (!cov->no_imag)
					{
          			  d2R(pm2,p,q) += 2*F[last]*F[i] * ( 
						( derD_raw1_a[i]*derDpq_raw1_a[pm2][p][q][1] + derDpq_raw1_a[pm2][p][q][i]*derD_raw1_a[1] + der2D_raw1_a[pm2][p][q][i]*a[0][1] + a[0][i]*der2D_raw1_a[pm2][p][q][1] + 
						  derD_raw1_b[i]*derDpq_raw1_b[pm2][p][q][1] + derDpq_raw1_b[pm2][p][q][i]*derD_raw1_b[1] + der2D_raw1_b[pm2][p][q][i]*b[0][1] + b[0][i]*der2D_raw1_b[pm2][p][q][1] )*cos_delta +
						( derD_raw1_a[i]*derDpq_raw1_b[pm2][p][q][1] + derDpq_raw1_a[pm2][p][q][i]*derD_raw1_b[1] + der2D_raw1_a[pm2][p][q][i]*b[0][1] + a[0][i]*der2D_raw1_b[pm2][p][q][1] 
						- derD_raw1_b[i]*derDpq_raw1_a[pm2][p][q][1] - derDpq_raw1_b[pm2][p][q][i]*derD_raw1_a[1] - der2D_raw1_b[pm2][p][q][i]*a[0][1] - b[0][i]*der2D_raw1_a[pm2][p][q][1] )*sin_delta ) ;
					}
    				else
          			  d2R(pm2,p,q) += 2*F[last]*F[i] * (
						derD_raw1_a[i]*derDpq_raw1_a[pm2][p][q][1] + derDpq_raw1_a[pm2][p][q][i]*derD_raw1_a[1] + der2D_raw1_a[pm2][p][q][i]*a[0][1] + a[0][i]*der2D_raw1_a[pm2][p][q][1] )*cos_delta ; 
					for ( int nD = 0;  nD < cov->Dass[pm2][p][q].toNum;  nD++ )
					{					  
  					  int pD = cov->Dass[pm2][p][q].to(nD,0);
  					  int qD = cov->Dass[pm2][p][q].to(nD,1);
					  if ( pD == qD )  {
						comp_aux4inv( a, b, derD_a, derD_b, last, i, pD, qD, k, l );
						add_d2invD_r_pp( der2D_a_last_i, a, b, last, i, pm2, pD, pm, k, l );
    				  }
    				  else 
						if (!cov->no_imag) {
						  comp_aux4inv( a, b, derD_a, derD_b, last, i, pD, qD, k, l );
						  add_d2invD_r( der2D_a_last_i, a, b, last, i, pm2, pD, qD, pm, k, l );
						  add_d2invD_i( der2D_b_last_i, a, b, last, i, pm2, pD, qD, pm, k, l );
						}
						else
						{
						  comp_aux4inv_no_imag( a, b, derD_a, derD_b, last, i, pD, qD, k, l );
						  add_d2invD_r_no_imag( der2D_a_last_i, a, b, last, i, pm2, pD, qD, pm, k, l );
						}
					}
				  }	
				  else	// if N_meas == 1
				  {	
					der2D_a_last_i = der2D_raw1_a[pm2][p][q][i];
		  			der2D_b_last_i = der2D_raw1_b[pm2][p][q][i];
				  }
				  if (!cov->no_imag)
					d2E(pm2,p,q) -= 2*F[i]*F[last]*(der2D_a_last_i*cos_delta-der2D_b_last_i*sin_delta);
				  else
					d2E(pm2,p,q) -= 2*F[i]*F[last]*der2D_a_last_i*cos_delta;
				}	// end of no_integ_last part
				
				for (int  j = i+1;  j < Num;  j++)
      			{
				  realnum der2D_aij = 0.,	der2D_bij = 0.;
				  realnum der2D_cij = 0.,	der2D_dij = 0.;
				  for ( int nD = 0;  nD < cov->Dass[pm2][p][q].toNum;  nD++ )
				  {	
  					int pD = cov->Dass[pm2][p][q].to(nD,0);
  					int qD = cov->Dass[pm2][p][q].to(nD,1);
    				if ( pD == qD ) {
					  comp_aux4inv( a, b, derD_a, derD_b, i, j, pD, qD, k, l );
					  add_d2invD_r_pp( der2D_aij, a, b, i, j, pm2, pD, pm, k, l );
					  if (!cov->no_imag) add_d2invD_i_pp( der2D_bij, a, b, i, j, pm2, pD, pm, k, l );
					}
					else 
    				  if (!cov->no_imag) {
						comp_aux4inv( a, b, derD_a, derD_b, i, j, pD, qD, k, l );
						add_d2invD_r( der2D_aij, a, b, i, j, pm2, pD, qD, pm, k, l );
						add_d2invD_i( der2D_bij, a, b, i, j, pm2, pD, qD, pm, k, l );
    				  }
					  else
					  {
						comp_aux4inv_no_imag( a, b, derD_a, derD_b, i, j, pD, qD, k, l );
						add_d2invD_r_no_imag( der2D_aij, a, b, i, j, pm2, pD, qD, pm, k, l );
					  }
					if (pD>=N_meas)
					{
    				  if ( pD == qD )
    				  {
						comp_aux4inv( c, d, derD_c, derD_d, i, j, pD, qD, k, l );
				  		add_d2invD_r_pp( der2D_cij, c, d, i, j, pm2, pD, pm, k, l );
				  		if (!cov->no_imag) add_d2invD_i_pp( der2D_dij, c, d, i, j, pm2, pD, pm, k, l );
					  }
					  else 
    					if (!cov->no_imag) {
						  comp_aux4inv( c, d, derD_c, derD_d, i, j, pD, qD, k, l );
				  		  add_d2invD_r( der2D_cij, c, d, i, j, pm2, pD, qD, pm, k, l );
				  		  add_d2invD_i( der2D_dij, c, d, i, j, pm2, pD, qD, pm, k, l );
						}
    					else
						{
						  comp_aux4inv_no_imag( c, d, derD_c, derD_d, i, j, pD, qD, k, l );
				  		  add_d2invD_r_no_imag( der2D_cij, c, d, i, j, pm2, pD, qD, pm, k, l );
						}
					}
				  }
				  realnum temp1 = ( der2D_aij - der2D_cij ) ;
    			  if (!cov->no_imag)
    			  {
        			realnum temp2 = ( der2D_bij - der2D_dij ) ;
        			d2E(pm2,p,q) -= 2*F[i]*F[j]* ( temp1*cos_del(i,j) - temp2*sin_del(i,j) ) ;
					if ( !no_integ_last || N_meas >= 2 )	
					{
        			  temp1 = derD_raw1_a[i]*derDpq_raw1_a[pm2][p][q][j] + derDpq_raw1_a[pm2][p][q][i]*derD_raw1_a[j] + 
							  derD_raw1_b[i]*derDpq_raw1_b[pm2][p][q][j] + derDpq_raw1_b[pm2][p][q][i]*derD_raw1_b[j] + 
							  der2D_raw1_a[pm2][p][q][i]*a[0][j] + a[0][i]*der2D_raw1_a[pm2][p][q][j] +
							  der2D_raw1_b[pm2][p][q][i]*b[0][j] + b[0][i]*der2D_raw1_b[pm2][p][q][j] ;
        			  temp2 = derD_raw1_b[i]*derDpq_raw1_a[pm2][p][q][j] + derDpq_raw1_b[pm2][p][q][i]*derD_raw1_a[j] + 
							- derD_raw1_a[i]*derDpq_raw1_b[pm2][p][q][j] - derDpq_raw1_a[pm2][p][q][i]*derD_raw1_b[j] + 
							  der2D_raw1_b[pm2][p][q][i]*a[0][j] + b[0][i]*der2D_raw1_a[pm2][p][q][j] +
							- der2D_raw1_a[pm2][p][q][i]*b[0][j] - a[0][i]*der2D_raw1_b[pm2][p][q][j] ;
        			  d2R(pm2,p,q) += 2*F[i]*F[j]* ( temp1*cos_del(i,j) + temp2*sin_del(i,j) ) ;				  
					}
				  }
				  else
    			  {
        			d2E(pm2,p,q) -= 2*F[i]*F[j]*temp1*cos_del(i,j) ;
					if ( !no_integ_last || N_meas >= 2 )	
					{
        			  temp1 = derD_raw1_a[i]*derDpq_raw1_a[pm2][p][q][j] + derDpq_raw1_a[pm2][p][q][i]*derD_raw1_a[j] + 
							  der2D_raw1_a[pm2][p][q][i]*a[0][j] + a[0][i]*der2D_raw1_a[pm2][p][q][j] ;
        			  d2R(pm2,p,q) += 2*F[i]*F[j]*temp1*cos_del(i,j) ;
					}
				  }
				}
			  }		// end of i,j looping
			}
		  }
		}		// end of p,q,pm2 looping
	  }		// end of do_der_D == 2 part

// end of preparations part; now calculate derivative wrt the Bef_term and then proceed
// with the parts requiring integration (in case of N_meas>=2)


      realnum Der_bef_term = 0;
	  for ( int nD = 0;  nD < NDact;  nD++ )
	  {			
		realnum diag_kl = 2; 	if ( kD[nD] == lD[nD] ) diag_kl = 1;
    	if ( kD[nD]>=N_meas ) {
		  Der_bef_term                    += diag_kl* cov->part[pm].deriv_D_r[ kD[nD] ][ lD[nD] ]*( c[kD[nD]][lD[nD]] - a[kD[nD]][lD[nD]] ); 
    	  if (!cov->no_imag) Der_bef_term += diag_kl* cov->part[pm].deriv_D_i[ kD[nD] ][ lD[nD] ]*( d[kD[nD]][lD[nD]] - b[kD[nD]][lD[nD]] );
		}
		else {
		  Der_bef_term 		    		  -= diag_kl* cov->part[pm].deriv_D_r[ kD[nD] ][ lD[nD] ]*a[kD[nD]][lD[nD]] ;
  		  if (!cov->no_imag) Der_bef_term -= diag_kl* cov->part[pm].deriv_D_i[ kD[nD] ][ lD[nD] ]*b[kD[nD]][lD[nD]] ;
		}
	  }



      if ( N_meas == 1 || no_integ_last ) 
      {
		if ( no_integ_last && N_meas == 1 )
      	  der_D[pm][k][l] = - Der_bef_term - dEdD ;
		else
      	  der_D[pm][k][l] = - Der_bef_term - dEdD - F0_isqrtR_sim*dRdD(pm,k,l) ;
		if ( cent )  der_D[pm][k][l] *= 0.5 ;
//cout << Der_bef_term << " " << dEdD <<" " << F0_isqrtR_sim*dRdD(pm,k,l) << endl;
//if (k==0&&l==1) cout << .5*dEdD << " " << .5*F[0]/Aux_Bessel*sim*dRdD(pm,k,l) << " " << .5*Der_bef_term << endl;
		if ( do_der_D == 2 ) 
		{
		  for (int  p = 0;  p <= k;  p++) 
		  {
			for (int  q = p+1;  (q < cov->DNum && p < k) || ( q <= l && p == k );  q++) 
			{
			  for (int pm2 = 0;  ( pm2 <= pm && p == k && q == l ) || ( pm2 < cov->N_part && (p < k || q < l) );  pm2++) 
			  {			
				realnum Der2_bef_term = 0; 
				for ( int nD = 0;  nD < cov->Dass[pm2][p][q].toNum;  nD++ )
				{
				  int pD = cov->Dass[pm2][p][q].to(nD,0);
				  int qD = cov->Dass[pm2][p][q].to(nD,1);
				  realnum diag_pp = 2; 	if ( pD == qD )	diag_pp *= .5;	
    			  if ( pD>=N_meas ) {
    				Der2_bef_term += diag_pp* (
					  cov->partpart[pm][pm2].deriv_D_r[p][q]*(c[pD][qD]-a[pD][qD]) + cov->part[pm2].deriv_D_r[pD][qD]*(derD_c[pD][qD]-derD_a[pD][qD]) );
					if (!cov->no_imag)  Der2_bef_term += diag_pp* (
					  cov->partpart[pm][pm2].deriv_D_i[p][q]*(d[pD][qD]-b[pD][qD]) + cov->part[pm2].deriv_D_i[pD][qD]*(derD_d[pD][qD]-derD_b[pD][qD]) );
				  }
				  else {
					Der2_bef_term -= diag_pp* ( 
					  cov->partpart[pm][pm2].deriv_D_r[p][q]*a[pD][qD] + cov->part[pm2].deriv_D_r[pD][qD]*derD_a[pD][qD] );
					if (!cov->no_imag)  Der2_bef_term -= diag_pp* ( 
  					  cov->partpart[pm][pm2].deriv_D_i[p][q]*b[pD][qD] + cov->part[pm2].deriv_D_i[pD][qD]*derD_b[pD][qD] );
  				  }
				}
				if ( no_integ_last && N_meas == 1 )
				  der2_D[pm][k][l][pm2][p][q] = - Der2_bef_term - d2E(pm2,p,q) ;
				else
				  der2_D[pm][k][l][pm2][p][q] = 
				  - Der2_bef_term - d2E(pm2,p,q) ;//- F0F0_iR_simder * dRdD(pm2,p,q)*dRdD(pm,k,l)
	    		  //- F0_isqrtR_sim * ( -dRdD(pm2,p,q)*dRdD(pm,k,l)*twoR_inv + d2R(pm2,p,q) );
				if ( cent )
				  der2_D[pm][k][l][pm2][p][q] *= .5;

				der2_D[pm2][p][q][pm][k][l] = der2_D[pm][k][l][pm2][p][q];
//	if (k==0&&l==1&&p==0&&q==1) cout << d2E(p,q) << " " << F0_isqrtR_sim *d2R(p,q) << " " << Der2_bef_term << endl;
			  }
			}
		  }
		}
//  if (k==0&&l==1)  cout << Der_bef_term << " " << Bef_term << " " << Aux2_der_D << " " << Aux1_der_D << " " << BessI1_term/BessI0_term << endl;
      }		// end of N_meas == 1 (or no_integ_last)  part

      if ( N_meas >= 2 && ! no_integ_last )
      {
        if ( improve_integ_now < 2 ) 
      	  dRdD(pm,k,l) += 2*F[1]*F[1]*( a[0][1]*derD_raw1_a[1] + b[0][1]*derD_raw1_b[1] ) ; 
		else 
		  if ( N_meas == 3 ) {
			dRdD(pm,k,l) += 2*F2sq_div_F1sq*F[2]*F[2]* a[1][2]*derD_raw2_a[2] ;
			dRdD_cos2 += 2*F2_div_F1*F[2]*F[2]* (a[0][2]*derD_raw2_a[2]+a[1][2]*derD_raw1_a[2]) ;
		  }
		if ( N_meas == 3 ) dRdD(pm,k,l) += 2*F[2]*F[2]* a[0][2]*derD_raw1_a[2] ;
		dEdD -= F[1]*F[1]*derD_a[1][1];
		if ( N_meas == 3 ) dEdD -= F[2]*F[2]*derD_a[2][2];
		dRdD_con = dRdD(pm,k,l);
		realnum tmp4dE12;	// only used if N_meas == 3
		realnum tmp4dR12;	// only used if N_meas == 3
		if ( N_meas == 3 && improve_integ_now < 2 )	// precomputing of terms needed for dE2dD, dRdD dependent on both 2. and 3. phase
		{
		  tmp4dE12 = 2*F[1]*F[2] * derD_a[1][2];
		  tmp4dR12 = 2*F[1]*F[2] * ( derD_raw1_a[1]*a[0][2] + a[0][1]*derD_raw1_a[2] ) ;
		}
		for (int  i = N_meas;  i < Num;  i++)				// precomputing of terms needed for dE2dD, dRdD
		{
		  realnum twoFlastFi = 2*F[last]*F[i];
		  tmp4dElast[i] = twoFlastFi * derD_a[last][i];
          tmp4dRlast[i] = twoFlastFi * ( 
			derD_raw1_a[i]*a[0][last] + a[0][i]*derD_raw1_a[last] + derD_raw1_b[i]*b[0][last] + b[0][i]*derD_raw1_b[last] ) ;
		  if ( !cov->no_imag )
		  {
			tmp4dElast2[i] = twoFlastFi * derD_b[last][i];
        	tmp4dRlast2[i] = twoFlastFi * ( 
			  derD_raw1_a[i]*b[0][last] + a[0][i]*derD_raw1_b[last] - derD_raw1_b[i]*a[0][last] - b[0][i]*derD_raw1_a[last] ) ;
		  }
		  if ( N_meas == 3 )
		  {
			realnum two_Flast2_Fi = 2*F[last2]*F[i];
			tmp4dElast2[i] = two_Flast2_Fi * derD_a[last2][i];
        	tmp4dRlast2[i] = two_Flast2_Fi * ( derD_raw1_a[i]*a[0][last2] + a[0][i]*derD_raw1_a[last2] ) ;
		  }
		}
		if ( do_der_D == 2 ) 
		{
		  for (int  p = 0;  p <= k;  p++) 
		  {
			for (int  q = p+1;  (q < cov->DNum && p < k) || ( q <= l && p == k );  q++) 
			{
			  for (int pm2 = 0;  ( pm2 <= pm && p == k && q == l ) || ( pm2 < cov->N_part && (p < k || q < l) );  pm2++) 
			  {			
    			der2_D[pm][k][l][pm2][p][q] = 0;
				d2R(pm2,p,q) += 2*F[1]*F[1]*( derDpq_raw1_a[pm2][p][q][1]*derD_raw1_a[1] + a[0][1]*der2D_raw1_a[pm2][p][q][1] +
										  derDpq_raw1_b[pm2][p][q][1]*derD_raw1_b[1] + b[0][1]*der2D_raw1_b[pm2][p][q][1] ) ;
				realnum derD_a22 = 0;
				for ( int nD = 0;  nD < cov->Dass[pm2][p][q].toNum;  nD++ )
				{
    			  int pD = cov->Dass[pm2][p][q].to(nD,0);
  				  int qD = cov->Dass[pm2][p][q].to(nD,1);
				  if ( pD == qD ) {
					comp_aux4inv( a, b, derD_a, derD_b, 1, 1, pD, qD, k, l );
					add_d2invD_r_pp( derD_a22, a, b, 1, 1, pm2, pD, pm, k, l );
				  }
				  else {
					if (!cov->no_imag)
					  { add_d2invdiagD_r( derD_a22, a, b, derD_a, derD_b, 1, pm2, pD, qD, pm, k, l ); }
					else
					  { add_d2invdiagD_r_no_imag( derD_a22, a, derD_a, 1, pm2, pD, qD, pm, k, l ); }
				  }
				}
				d2E(pm2,p,q) -= F[1]*F[1]*derD_a22 ;
        		for (int  i = N_meas;  i < Num;  i++)				// precomputing of terms needed for int_d2E2, int_d2R
				{
				  realnum twoF0Fi = 2*F[1]*F[i];
				  tmp4d2E1(pm2,i,p,q) = tmp4d2E2(pm2,i,p,q) = tmp4d2R1(pm2,i,p,q) = tmp4d2R2(pm2,i,p,q) = 0;
				  for ( int nD = 0;  nD < cov->Dass[pm2][p][q].toNum;  nD++ )
				  {
    				int pD = cov->Dass[pm2][p][q].to(nD,0);
  					int qD = cov->Dass[pm2][p][q].to(nD,1);
					if ( pD == qD ) {
					  comp_aux4inv( a, b, derD_a, derD_b, 1, i, pD, qD, k, l ) ;
  	  				  add_d2invD_r_pp( tmp4d2E1(pm2,i,p,q), a, b, 1, i, pm2, pD, pm, k, l ) ;
          			  if (!cov->no_imag) add_d2invD_i_pp( tmp4d2E2(pm2,i,p,q), a, b, 1, i, pm2, pD, pm, k, l ) ;
					}
					else 
    				  if (!cov->no_imag) {
						comp_aux4inv( a, b, derD_a, derD_b, 1, i, pD, qD, k, l ) ;
  	  					add_d2invD_r( tmp4d2E1(pm2,i,p,q), a, b, 1, i, pm2, pD, qD, pm, k, l ) ;
          				add_d2invD_i( tmp4d2E2(pm2,i,p,q), a, b, 1, i, pm2, pD, qD, pm, k, l ) ;
					  }
					  else
					  {
						comp_aux4inv_no_imag( a, b, derD_a, derD_b, 1, i, pD, qD, k, l ) ;
  	  					add_d2invD_r_no_imag( tmp4d2E1(pm2,i,p,q), a, b, 1, i, pm2, pD, qD, pm, k, l ) ;
					  }
				  }
				  tmp4d2E1(pm2,i,p,q) *= twoF0Fi;
    			  if (!cov->no_imag)
    			  {
					tmp4d2E2(pm2,i,p,q) *= twoF0Fi;
          			tmp4d2R1(pm2,i,p,q) = twoF0Fi * ( 
					  derD_raw1_a[i]*derDpq_raw1_a[pm2][p][q][1] + derDpq_raw1_a[pm2][p][q][i]*derD_raw1_a[1] + der2D_raw1_a[pm2][p][q][i]*a[0][1] + a[0][i]*der2D_raw1_a[pm2][p][q][1] + 
					  derD_raw1_b[i]*derDpq_raw1_b[pm2][p][q][1] + derDpq_raw1_b[pm2][p][q][i]*derD_raw1_b[1] + der2D_raw1_b[pm2][p][q][i]*b[0][1] + b[0][i]*der2D_raw1_b[pm2][p][q][1] ) ;
          			tmp4d2R2(pm2,i,p,q) = twoF0Fi * ( 
					  derD_raw1_a[i]*derDpq_raw1_b[pm2][p][q][1] + derDpq_raw1_a[pm2][p][q][i]*derD_raw1_b[1] + der2D_raw1_a[pm2][p][q][i]*b[0][1] + a[0][i]*der2D_raw1_b[pm2][p][q][1] + 
					- derD_raw1_b[i]*derDpq_raw1_a[pm2][p][q][1] - derDpq_raw1_b[pm2][p][q][i]*derD_raw1_a[1] - der2D_raw1_b[pm2][p][q][i]*a[0][1] - b[0][i]*der2D_raw1_a[pm2][p][q][1] ) ;
				  }
    			  else
          			tmp4d2R1(pm2,i,p,q) = twoF0Fi * ( 
					  derD_raw1_a[i]*derDpq_raw1_a[pm2][p][q][1] + derDpq_raw1_a[pm2][p][q][i]*derD_raw1_a[1] + der2D_raw1_a[pm2][p][q][i]*a[0][1] + a[0][i]*der2D_raw1_a[pm2][p][q][1] ) ; 
				}
			  }
			}
		  }
		}

// finally doing the integration here
		for ( int n = n_bottom+1;  n < n_top;  n++ )
//        for ( int n = 0;  n < Gauss_number_now_2;  n++ ) 
// temporary solution before the real time saving algorithm is implemented
//if ( improve_integ<2 || (n<Gauss_number_now_2*.75&&n>=Gauss_number_now_2*.25) || N_meas<3 ) 
//if (n<Gauss_number_now_2*.51&&n>=Gauss_number_now_2*.485) 
        for ( int m = m_bottom;  m < m_top;  (m==m_rej_bottom-1) ? m+=m_rej_num+1 : m++ ) 	
		if ( (N_meas==2&&for_derD___weight_E2_BessI0[m]!=0.) || 
			 (N_meas==3&&for_derD___weight_E2_BessI0_siras[n][m]!=0.) )
		{
		  realnum R_inv;
		  dRdD(pm,k,l) = dRdD_con;
		  dE2dD = 0;
		  if ( improve_integ_now >= 2 && N_meas == 2 ) 
		  {
			realnum cos_2last=cos_last[m]*cos_last[m]-sin_last[m]*sin_last[m];
			realnum sin_2last=2*sin_last[m]*cos_last[m];
			dE2dD -= cos_2last*F[0]*F[1]*2*derD_a[0][1];
			dRdD(pm,k,l) += cos_2last* dRdD_cos2 + sin_2last* dRdD_sin2 ;
		  }
		  else
		  {
    		if (!cov->no_imag)
        	  for (int  i = N_meas;  i < Num;  i++)
			  {
          		dE2dD -= tmp4dElast[i]*coslast_del[m][i] - tmp4dElast2[i]*sinlast_del[m][i] ;
	  			dRdD(pm,k,l) += tmp4dRlast[i]*coslast_del[m][i] + tmp4dRlast2[i]*sinlast_del[m][i] ;
			  }
			else
        	  for (int  i = N_meas;  i < Num;  i++)
			  {
          		dE2dD -= tmp4dElast[i]*coslast_del[m][i] ;
	  			dRdD(pm,k,l) += tmp4dRlast[i]*coslast_del[m][i] ;
			  }
		  }
		  if ( N_meas == 3 ) 
		  {
			if ( improve_integ_now >= 2 ) 
			{
			  realnum cos_2last2=cos_last2[n]*cos_last2[n]-sin_last2[n]*sin_last2[n];
			  realnum sin_2last2=2*sin_last2[n]*cos_last2[n];
			  dE2dD -= cos_2last2*F[0]*F[1]*2*derD_a[0][1];
			  dRdD(pm,k,l) += cos_2last2* dRdD_cos2 + sin_2last2* dRdD_sin2 ;
        	  for (int  i = N_meas;  i < Num;  i++)
			  {
				dRdD(pm,k,l) += 2 *F2sq_div_F1sq*F[2]*F[i]*(derD_raw2_a[i]*a[1][2]+a[1][i]*derD_raw2_a[2])*coslast_del[m][i] + 
			  	  cos_2last2*( 2*F[i]*F[2]*F2_div_F1*( derD_raw1_a[i]*a[1][2] + a[0][i]*derD_raw2_a[2] +
													   derD_raw2_a[i]*a[0][2] + a[1][i]*derD_raw1_a[2] ) * coslast_del[m][i] ) +
			  	  sin_2last2*( 2*F[i]*F[2]*F2_div_F1*( derD_raw1_a[i]*a[1][2] + a[0][i]*derD_raw2_a[2] -
													   derD_raw2_a[i]*a[0][2] - a[1][i]*derD_raw1_a[2] ) * sinlast_del[m][i] );
			  }
			}
			else 
			{
        	  for (int  i = N_meas;  i < Num;  i++)
			  {
          		dE2dD -= tmp4dElast2[i]*coslast2_del[n][i] ;
	  			dRdD(pm,k,l) += tmp4dRlast2[i]*coslast2_del[n][i] ;
			  }			
          	  dE2dD -= tmp4dE12*cos_last12_del[n][m] ;
	  		  dRdD(pm,k,l) += tmp4dR12*cos_last12_del[n][m] ;
			}
		  }
		  realnum dRdDkl = dRdD(pm,k,l);
		  if ( do_der_D == 2 ) 
		  {
			int_dRdD(pm,m,k,l) = dRdDkl;
			int_dE2dD(pm,m,k,l) = dE2dD;
			R_inv = 1./R[m];
			realnum dRdDkl_F1_F1_iR = F[0]*F[0]*dRdDkl*R_inv ;
			for (int  p = 0;  p <= k;  p++) 
			{
			  for (int  q = p+1;  (q < cov->DNum && p < k) || ( q <= l && p == k );  q++) 
			  {
				for (int pm2 = 0;  ( pm2 <= pm && p == k && q == l ) || ( pm2 < cov->N_part && (p < k || q < l) );  pm2++) 
				{
				  realnum int_d2E2pq = 0;
				  realnum int_d2Rpq = d2R(pm2,p,q);
    			  if (!cov->no_imag)
        			for (int  i = N_meas;  i < Num;  i++)
					{
					  int_d2E2pq -= tmp4d2E1(pm2,i,p,q)*coslast_del[m][i] - tmp4d2E2(pm2,i,p,q)*sinlast_del[m][i] ;
	  				  int_d2Rpq += tmp4d2R1(pm2,i,p,q)*coslast_del[m][i] + tmp4d2R2(pm2,i,p,q)*sinlast_del[m][i] ;
					}
				  else
        			for (int  i = N_meas;  i < Num;  i++)
					{
					  int_d2E2pq -= tmp4d2E1(pm2,i,p,q)*coslast_del[m][i] ;
	  				  int_d2Rpq += tmp4d2R1(pm2,i,p,q)*coslast_del[m][i] ;
					}
				  realnum int_dRdDpq = int_dRdD(pm2,m,p,q);
				  realnum int_dE2dDpq = int_dE2dD(pm2,m,p,q);
				  der2_D[pm][k][l][pm2][p][q] += 
					for_derD___F0_isqrtR_weight_E2_BessI1[m]*
					  ( int_dE2dDpq*dRdDkl + dE2dD*int_dRdDpq - int_dRdDpq*dRdDkl*R_inv + int_d2Rpq ) +
			        		  for_derD___weight_E2_BessI0[m]*
					  ( int_dE2dDpq*dE2dD + int_d2E2pq + int_dRdDpq*dRdDkl_F1_F1_iR ) ;
				}
			  }
			}
		  }
		  if ( N_meas == 2 ) 
        	der_D[pm][k][l] += for_derD___F0_isqrtR_weight_E2_BessI1[m]*dRdDkl +
  		                           for_derD___weight_E2_BessI0[m]*dE2dD ;
		  else	// N_meas == 3
        	der_D[pm][k][l] += for_derD___F0_isqrtR_weight_E2_BessI1_siras[n][m]*dRdDkl +
  		                           for_derD___weight_E2_BessI0_siras[n][m]*dE2dD ;
		}	// end of integration
//	if (k==0 && l==1) cout << Der_bef_term << " " << dEdD << " " << der_D[k][l]/integ << " " << integ << endl;
		if ( do_der_D == 2 ) 
		{
		  dIdD(pm,k,l) = der_D[pm][k][l];
		  for (int  p = 0;  p <= k;  p++) 
		  {
			for (int  q = p+1;  (q < cov->DNum && p < k) || ( q <= l && p == k );  q++) 
			{
			  for (int pm2 = 0;  ( pm2 <= pm && p == k && q == l ) || ( pm2 < cov->N_part && (p < k || q < l) );  pm2++) 
			  {			
				realnum Der2_bef_term = 0; 
				for ( int nD = 0;  nD < cov->Dass[pm2][p][q].toNum;  nD++ )
				{
				  int pD = cov->Dass[pm2][p][q].to(nD,0);
				  int qD = cov->Dass[pm2][p][q].to(nD,1);
				  realnum diag_pp = 2; 	if ( pD == qD )	diag_pp *= .5;	
    			  if ( pD>=N_meas ) { 
    				Der2_bef_term += diag_pp* (
					  cov->partpart[pm][pm2].deriv_D_r[pD][qD]*(c[pD][qD]-a[pD][qD]) + cov->part[pm2].deriv_D_r[pD][qD]*(derD_c[pD][qD]-derD_a[pD][qD]) );
					if (!cov->no_imag)  Der2_bef_term += diag_pp* (
					  cov->partpart[pm][pm2].deriv_D_i[pD][qD]*(d[pD][qD]-b[pD][qD]) + cov->part[pm2].deriv_D_i[pD][qD]*(derD_d[pD][qD]-derD_b[pD][qD]) );
				  }
				  else {
					Der2_bef_term -= diag_pp* ( 
					  cov->partpart[pm][pm2].deriv_D_r[pD][qD]*a[pD][qD] + cov->part[pm2].deriv_D_r[pD][qD]*derD_a[pD][qD] );
					if (!cov->no_imag)  Der2_bef_term -= diag_pp* ( 
  					  cov->partpart[pm][pm2].deriv_D_i[pD][qD]*b[pD][qD] + cov->part[pm2].deriv_D_i[pD][qD]*derD_b[pD][qD] );
  				  }
				}
//	if (k==0&& l==1&& p==0 && q==1 && pm==1 &&pm2==0) cout << -Der2_bef_term << " " << d2E(pm2,p,q) << " " << der2_D[pm][k][l][pm2][p][q]*integ_inv << " " << dIdD(pm2,p,q)*dIdD(pm,k,l)*integ_inv*integ_inv << endl;
//	if (k==0&& l==1&& p==0 && q==1 && pm==1 &&pm2==1) cout << -Der2_bef_term << " " << d2E(pm2,p,q) << " " << der2_D[pm][k][l][pm2][p][q]*integ_inv << " " << dIdD(pm2,p,q)*dIdD(pm,k,l)*integ_inv*integ_inv << endl;
	  			der2_D[pm][k][l][pm2][p][q] = 
				  - Der2_bef_term - d2E(pm2,p,q) - ( der2_D[pm][k][l][pm2][p][q]*integ - dIdD(pm2,p,q)*dIdD(pm,k,l) )*integ_inv*integ_inv ;
				der2_D[pm2][p][q][pm][k][l] = der2_D[pm][k][l][pm2][p][q];
			  }
			}
		  }
		}
//	if (k==0&& l==1&& pm==3)		cout << - Der_bef_term << " " << - dEdD << " " << - der_D[pm][k][l] << " " << integ_inv << endl; // don't forget at e2argmax when comparing the values
		der_D[pm][k][l] = - Der_bef_term - dEdD - der_D[pm][k][l]*integ_inv ;
	  }
	  }
    }
  }
  }
//  cout.precision(10);    
//  cout << "2. der.: " << der2_phF[Num-1][Num-1] << " " << der2_ph2[Num-1][Num-1] << " " << der2_F2[Num-1][Num-1] << endl;  
//  cout << "1. der.:  " << der_F[Num-1] << " " << der_ph[Num-1] << endl;
//  cout << "1. der. wrt D:" << der_D[0][2] + der_D[0][3] + der_D[1][2] + der_D[1][3] << " " << der_D[0][1] << endl;  
//  cout << "2. der. wrt D:" << der2_D[0][2][0][2] + der2_D[0][3][0][2] + der2_D[1][2][0][2] + der2_D[1][3][0][2] +
//	der2_D[0][2][0][3] + der2_D[0][3][0][3] + der2_D[1][2][0][3] + der2_D[1][3][0][3] + 
//	der2_D[0][2][1][2] + der2_D[0][3][1][2] + der2_D[1][2][1][2] + der2_D[1][3][1][2] +
//	der2_D[0][2][1][3] + der2_D[0][3][1][3] + der2_D[1][2][1][3] + der2_D[1][3][1][3] << " " << der2_D[0][1][0][1] << endl;  

//  fval = fval * F[0]/PI/cov->re[0][0]*tab->ExpM(F[0]*F[0]/cov->re[0][0]);


// changing back the variables switched in order to improve the integration
//if (N_meas>2)
  if (improve_integ_now||num4phib>=0)
	ChangeIntegOrder( F, phase, 1 );

  return fval;
}


#undef d2E
#undef d2R
#undef comp_dinvD_r
#undef comp_dinvD_i
#undef comp_dinvdiagD_r 
#undef comp_dinvdiagD_i
#undef comp_aux4inv
#undef comp_aux4inv_no_imag
#undef add_d2invD_r
#undef add_d2invD_i
#undef add_d2invdiagD_r
#undef add_d2invdiagD_i
#undef comp_dinvD_r_no_imag
#undef comp_dinvdiagD_r_no_imag
#undef add_d2invD_r_no_imag
#undef add_d2invdiagD_r_no_imag
#undef int_dE2dD
#undef int_dRdD
#undef tmp4d2E1
#undef tmp4d2E2
#undef tmp4d2R1
#undef tmp4d2R2
#undef sin_del
#undef cos_del
#undef dIdD
#undef dRdD

#undef to



// interface for calling lapack for eigenvalue and eigenvectors calculation
extern "C" void zheevd_(char*,char*,int*,complex<double>*,int*,double*,complex<double>*,int*,double*,int*,int*,int*,int*);
extern "C" void cheevd_(char*,char*,int*,complex<float>*,int*,float*,complex<float>*,int*,float*,int*,int*,int*,int*);
extern "C" void dsyev_(char*,char*,int*,double*,int*,double*,double*,int*,/*int*,int*,*/int*);
extern "C" void ssyev_(char*,char*,int*,float*,int*,float*,float*,int*,/*int*,int*,*/int*);
extern "C" void dsyevd_(char*,char*,int*,double*,int*,double*,double*,int*,int*,int*,int*);
extern "C" void ssyevd_(char*,char*,int*,float*,int*,float*,float*,int*,int*,int*,int*);

template<typename realnum>
void lapack_eigen_trait_complex(char*a1,char*a2,int*a3,complex<realnum>*a4,int*a5,realnum*a6,complex<realnum>*a7,int*a8,realnum*a9,int*a10,int*a11,int*a12,int*a13)
{ };
template<>
void lapack_eigen_trait_complex(char*a1,char*a2,int*a3,complex<double>*a4,int*a5,double*a6,complex<double>*a7,int*a8,double*a9,int*a10,int*a11,int*a12,int*a13)	
{	
  zheevd_(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13);  
};
template<>
void lapack_eigen_trait_complex(char*a1,char*a2,int*a3,complex<float>*a4,int*a5,float*a6,complex<float>*a7,int*a8,float*a9,int*a10,int*a11,int*a12,int*a13)	
{	
  cheevd_(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13);  
};


template<typename realnum>
void lapack_eigen_trait_real(char*a1,char*a2,int*a3,realnum*a4,int*a5,realnum*a6,realnum*a7,int*a8,int*a9,int*a10,int*a11)
{ };
template<>
inline void lapack_eigen_trait_real(char*a1,char*a2,int*a3,double*a4,int*a5,double*a6,double*a7,int*a8,int*a9,int*a10,int*a11)
{ 
  dsyevd_(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11);
//  dsyev_(a1,a2,a3,a4,a5,a6,a7,a8,/*a9,a10,*/a11);
};
template<>
inline void lapack_eigen_trait_real(char*a1,char*a2,int*a3,float*a4,int*a5,float*a6,float*a7,int*a8,int*a9,int*a10,int*a11)
{ 
  ssyevd_(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11);
//  ssyev_(a1,a2,a3,a4,a5,a6,a7,a8,/*a9,a10,*/a11);
};



// Finds the inverses to cov matrices and their eigenvalues

template<typename realnum>
void likelihood<realnum>::InverseAndEigen()
{
  int Num_mod = Num - N_meas;
  realnum** inp_r = work.space2;
  realnum** inp_i = work.space2 + Num;


  if (no_imag)
  {
	for (int i = 0;  i < Num_mod;  i++) 
  	  for (int j = i;  j < Num_mod;  j++) 
    	lpck.inout2_r[i+Num_mod*j] = cov->re[i+N_meas][j+N_meas];
//	lapack_eigen_trait_real<realnum>( &lpck.ch1, &lpck.ch2, &Num_mod, lpck.inout2_r, &Num_mod, eigenvalues2[Rice], 
//		  lpck.work_r, &lpck.dim_w, /*lpck.iwork, &lpck.dim_iw,*/ &lpck.info);
	lapack_eigen_trait_real<realnum>( &lpck.ch1, &lpck.ch2, &Num_mod, lpck.inout2_r, &Num_mod, eigenvalues2[Rice], 
		  lpck.work_r, &lpck.dim_w, lpck.iwork, &lpck.dim_iw, &lpck.info);
	for (int i = 0;  i < Num_mod;  i++) 
	{
	  inv_eigenval[i] = 1./eigenvalues2[Rice][i];
	  for (int j = 0;  j < Num_mod;  j++)  {
		inp_r[i][j] = lpck.inout2_r[j+Num_mod*i];
		inp_i[i][j] = 0.;
	  }
	}
  }
  else
  {
	for (int i = 0;  i < Num_mod;  i++) 
  	  for (int j = i;  j < Num_mod;  j++) 
    	lpck.inout2_c[i+Num_mod*j] = complex<realnum>( cov->re[i+N_meas][j+N_meas], cov->im[i+N_meas][j+N_meas] );
	lapack_eigen_trait_complex<realnum>( &lpck.ch1, &lpck.ch2, &Num_mod, lpck.inout2_c, &Num_mod, eigenvalues2[Rice],
		  lpck.work_c, &lpck.dim_w, lpck.rwork, &lpck.dim_rw, lpck.iwork, &lpck.dim_iw, &lpck.info);
	for (int i = 0;  i < Num_mod;  i++) 
	{
	  inv_eigenval[i] = 1./eigenvalues2[Rice][i];
	  for (int j = 0;  j < Num_mod;  j++)        
	  {
		inp_r[i][j] = real(lpck.inout2_c[j+Num_mod*i]);
		inp_i[i][j] = imag(lpck.inout2_c[j+Num_mod*i]);
	  }
	}
  }

  InverseMatrix( Num_mod, inp_r, inp_i, c, d, eigenvalues2[Rice], inv_eigenval );

  if (no_imag)
  {
	for (int i = 0;  i < Num;  i++) 
  	  for (int j = i;  j < Num;  j++) 
		lpck.inout1_r[i+Num*j] = cov->re[i][j];
//	lapack_eigen_trait_real<realnum>( &lpck.ch1, &lpck.ch2, &Num, lpck.inout1_r, &Num, eigenvalues1[Rice], 
//		  lpck.work_r, &lpck.dim_w, /*lpck.iwork, &lpck.dim_iw,*/ &lpck.info);
	lapack_eigen_trait_real<realnum>( &lpck.ch1, &lpck.ch2, &Num, lpck.inout1_r, &Num, eigenvalues1[Rice], 
		  lpck.work_r, &lpck.dim_w, lpck.iwork, &lpck.dim_iw, &lpck.info);
	for (int i = 0;  i < Num;  i++) 
	{
	  inv_eigenval[i] = 1./eigenvalues1[Rice][i];
	  for (int j = 0;  j < Num;  j++) {
		inp_r[i][j] = lpck.inout1_r[j+Num*i];
		inp_i[i][j] = 0.;
	  }
	}
  }
  else
  {
	for (int i = 0;  i < Num;  i++) 
  	  for (int j = i;  j < Num;  j++) 
		lpck.inout1_c[i+Num*j] = complex<realnum>( cov->re[i][j], cov->im[i][j] );
	lapack_eigen_trait_complex<realnum>( &lpck.ch1, &lpck.ch2, &Num, lpck.inout1_c, &Num, eigenvalues1[Rice], 
		  lpck.work_c, &lpck.dim_w, lpck.rwork, &lpck.dim_rw, lpck.iwork, &lpck.dim_iw, &lpck.info);
	for (int i = 0;  i < Num;  i++) 
	{
	  inv_eigenval[i] = 1./eigenvalues1[Rice][i];
	  for (int j = 0;  j < Num;  j++)        
	  {
		inp_r[i][j] = real(lpck.inout1_c[j+Num*i]);
		inp_i[i][j] = imag(lpck.inout1_c[j+Num*i]);
	  }
	}
  }

  InverseMatrix( Num, inp_r, inp_i, a, b, eigenvalues1[Rice], inv_eigenval );
  
}



template<typename realnum>
realnum likelihood<realnum>::TransfToProb()
{
  realnum prob; 
  if ( Tab_func ) 	prob = tab->ExpM(fval);
  else				prob = exp(-fval);
  realnum prob_inv = 0;
  if ( prob > 0. ) 	prob_inv = 1./prob;
  if ( do_der_Fph )
  {
    for (int i = N_meas;  i < Num;  i++) 
    {
	  der_F[i] *= -prob;
	  der_ph[i] *= -prob;
	  if ( do_der_Fph == 2 )
  	  for (int j = N_meas;  j < Num;  j++) 
  	  {
		der2_F2[i][j] = prob_inv*der_F[i]*der_F[j] - prob*der2_F2[i][j];
		der2_phF[i][j] = prob_inv*der_ph[i]*der_F[j] - prob*der2_phF[i][j];
		der2_ph2[i][j] = prob_inv*der_ph[i]*der_ph[j] - prob*der2_ph2[i][j];
	  }
	}
  }
  if ( do_der_D )
  {
	for (int  k = 0;  k < cov->DNum;  k++)
  	  for (int  l = k+1;  l < cov->DNum;  l++)
		for ( int pm = 0;  pm < cov->N_part;  pm++ )			// now looping over k,l,pm to get derivatives wrt D[pm,k,l]
		{
		  der_D[pm][k][l] *= -prob;
		  if ( do_der_D == 2 )
		  {
			for (int  p = 0;  p <= k;  p++) 
			  for (int  q = p+1;  (q < cov->DNum && p < k) || ( q <= l && p == k );  q++) 
				for (int pm2 = 0;  ( pm2 <= pm && p == k && q == l ) || ( pm2 < cov->N_part && (p < k || q < l) );  pm2++) 
				{
				  der2_D[pm][k][l][pm2][p][q] = 
					prob_inv*der_D[pm][k][l]*der_D[pm2][p][q] - prob*der2_D[pm][k][l][pm2][p][q];
				  der2_D[pm2][p][q][pm][k][l] = der2_D[pm][k][l][pm2][p][q];
				}
		  }
		}
  }
  return fval = prob;
}

// prints real part of actual inverse matrices (array a or c)
template<typename realnum>
void likelihood<realnum>::PrintInvMat( int modelmat )
{
  cout << "Inverted cov. matrix "; 
  if ( modelmat == 1 ) cout << "(model's matrix)"; 
  cout << " :" << endl;
  realnum **mat, **mat_im;
  int dim_mat;
  if ( modelmat == 1 ) {
	mat = c;
	mat_im = d;
	dim_mat = Num - N_meas;
  }
  else if ( modelmat == 0 ) {
	mat = a;
	mat_im = b;
	dim_mat = Num;
  }
  else cout << "Couldn't print requested inverted matrix - wrong matrix number requested.";

  for (int i=0; i<dim_mat; i++)
  {
	for (int j=0; j<dim_mat; j++)
	  if ( no_imag ) cout << setw(12) << mat[i][j] << " ";
	  else cout << setw(24) << "(" << mat[i][j] << "," << mat_im[i][j] << " ";
	cout << endl;
  }
}

// prints actual eigenvalues (of matrix a or c)
template<typename realnum>
void likelihood<realnum>::PrintEigenv( int modelmat )
{
  cout << "Eigenvalues of cov. matrix "; 
  if ( modelmat == 1 ) cout << "(model's matrix)"; 
  cout << " :" << endl;
  realnum *eigenv;
  int dim_mat;
  if ( modelmat == 1 ) {
	eigenv = eigenvalues2[0];
	dim_mat = Num - N_meas;
  }
  else if ( modelmat == 0 ) {
	eigenv = eigenvalues1[0];
	dim_mat = Num;
  }
  else cout << "Couldn't print requested eigenvalues - wrong matrix number requested.";

  for (int i=0; i<dim_mat; i++)
	cout << setw(12) << eigenv[i] << " ";
  cout << endl;
}


// deprecated
// the mixed 2.derivatives wrt D par. of func1 and func2 are not computed 
// (easy to do but would take a considerable amount of either time or memory and I don't use them at the moment anyway)
template<typename realnum>
realnum likelihood<realnum>::EvaluateSR_MLD( realnum *F1, realnum *phase1, realnum *F2, realnum *phase2 ) // temporary, change structure later !!!!!
{
  no_integ_last = 1;		// just to be sure here; should be ensured better in the future

  if (Rice)
  {
	no_integ_last = 0;
	cov->Make_matrix();		//	this stuff should not be here! move stuff belonging to covmat to covmat.h and let all this
	InverseAndEigen();		//  to be called by calling function!
	realnum fvalue = EvaluateSR(F1,phase1);	
	no_integ_last = 1;
	return fvalue;
  }
  else
  {

// set all auxiliary variables to 0
  SetCentRice(0,0);
  realnum fvalue = 0;
  if ( do_der_Fph )
  {
  	for (int i = 0;  i < Num;  i++) 
  	{
	  auxder_F[i] = auxder_ph[i] = 0;
	  if ( do_der_Fph == 2 )
  	  for (int j = 0;  j < Num;  j++) 
		auxder2_F2[i][j] = auxder2_phF[i][j] = auxder2_ph2[i][j] = 0;
	}
  }
  if ( do_der_D )
  {
	for (int  k = 0;  k < cov->DNum;  k++)
  	  for (int  l = k+1;  l < cov->DNum;  l++)
		for ( int pm = 0;  pm < cov->N_part;  pm++ )			// now looping over k,l,pm to get derivatives wrt D[pm,k,l]
		{
		  auxder_D[pm][k][l] = 0;
		  if ( do_der_D == 2 )
		  {
			for (int  p = 0;  p <= k;  p++) 
			  for (int  q = p+1;  (q < cov->DNum && p < k) || ( q <= l && p == k );  q++) 
				for (int pm2 = 0;  ( pm2 <= pm && p == k && q == l ) || ( pm2 < cov->N_part && (p < k || q < l) );  pm2++) 
				  auxder2_D[pm][k][l][pm2][p][q] = auxder2_D[pm2][p][q][pm][k][l] = 0.;
		  }
		}
  }

// backup of original setting
  int real_do_der_D = do_der_D;
  int real_do_der_Fph = do_der_Fph;
  int real_N_part = cov->N_part;

// precomputing func1 values
  realnum* func1 = new realnum[Gauss_number];	// allocate in Alloc later
  int NumDpar_func1 = 1;						// this variable should probably become part of llhood class for no_integ_part (and be set up at beginning so as N_part)
  SetCentRice(0,1);
  do_der_D = do_der_Fph = 0;
  cov->N_part = NumDpar_func1;
  cov->Make_matrix();
  InverseAndEigen();
  realnum P_F0;  			// the "redefinition constant" (1/P(F,al)) of MLDR (doesn't affect derivatives)
  if (Tab_func) P_F0 = F1[0]/PI/cov->re[0][0]*tab->ExpM(F1[0]*F1[0]/cov->re[0][0]);
  else			P_F0 = F1[0]/PI/cov->re[0][0]*exp(-F1[0]*F1[0]/cov->re[0][0]);
  for ( nil_m = 0;  nil_m < Gauss_number;  nil_m++ )
  {
	EvaluateSR(F1,phase1);
	func1[nil_m] = TransfToProb();
  }

// computing the derivatives wrt D from func2, precomputing func2 values and calculating FOM
// also coputing derivatives wrt F,ph of the second function (heavy atoms F,ph)
  realnum* func2 = new realnum[Gauss_number];	// allocate in Alloc later
  realnum FOM2sin = 0, FOM2cos = 0, FOM2 = 0;
  SetCentRice(0,0);
  do_der_D = real_do_der_D;
  do_der_Fph = real_do_der_Fph;		// this can be removed if don't want to compute heavy atom derivatives!
  cov->N_part = real_N_part - NumDpar_func1;
  cov->Make_matrix();
  InverseAndEigen();
  for ( nil_m = 0;  nil_m < Gauss_number;  nil_m++ )
  {
	EvaluateSR(F2,phase2);
	func2[nil_m] = TransfToProb();
	fvalue +=	Gauss_weight[nil_m]*func1[nil_m]*func2[nil_m];
	if ( do_der_Fph )
	{
  	  for (int i = N_meas;  i < Num;  i++) 
  	  {
		auxder_F[i] += Gauss_weight[nil_m]*der_F[i]*func1[nil_m];
		auxder_ph[i] += Gauss_weight[nil_m]*der_ph[i]*func1[nil_m];
		if ( do_der_Fph == 2 )
  		for (int j = N_meas;  j < Num;  j++) 
  		{
		  auxder2_F2[i][j] += Gauss_weight[nil_m]*der2_F2[i][j]*func1[nil_m];
		  auxder2_phF[i][j] += Gauss_weight[nil_m]*der2_phF[i][j]*func1[nil_m];
		  auxder2_ph2[i][j] += Gauss_weight[nil_m]*der2_ph2[i][j]*func1[nil_m];
		}
	  }
	}
	if ( do_der_D )
	{
	  for (int  k = 0;  k < cov->DNum;  k++)
  		for (int  l = k+1;  l < cov->DNum;  l++)
		  for ( int pm = 0;  pm < cov->N_part;  pm++ )		// now looping over k,l,pm to get derivatives wrt D[pm,k,l]				
		  {		// 			|		carefully - the NumDpar_func1 change is here!
			auxder_D[pm+NumDpar_func1][k][l] += Gauss_weight[nil_m]*der_D[pm][k][l]*func1[nil_m];
			if ( do_der_D == 2 )
			{
			  for (int  p = 0;  p <= k;  p++) 
				for (int  q = p+1;  (q < cov->DNum && p < k) || ( q <= l && p == k );  q++) 
				  for (int pm2 = 0;  ( pm2 <= pm && p == k && q == l ) || ( pm2 < cov->N_part && (p < k || q < l) );  pm2++) 
						// 				|	                   |   carefully - the NumDpar_func1 changes are here!
					auxder2_D[pm+NumDpar_func1][k][l][pm2+NumDpar_func1][p][q] += Gauss_weight[nil_m]*der2_D[pm][k][l][pm2][p][q]*func1[nil_m];
			}
		  }
	}
//    FOM2sin += Gauss_weight[nil_m]*func1[nil_m]*func2[nil_m] * sin_last[nil_m];
//    FOM2cos += Gauss_weight[nil_m]*func1[nil_m]*func2[nil_m] * cos_last[nil_m];
    realnum cos_delta;
    if (Tab_func) 	cos_delta = tab->Cos(phase1[1])*cos_last[nil_m] + tab->Sin_charged(phase1[1])*sin_last[nil_m] ;    
	else			cos_delta = cos(phase1[1])*cos_last[nil_m] + sin(phase1[1])*sin_last[nil_m] ;    
    FOM2 += Gauss_weight[nil_m]*func1[nil_m]*func2[nil_m] * cos_delta;
  }
  FOM = FOM2 / fvalue;  
//  FOM2 = sqrt( FOM2sin*FOM2sin + FOM2cos*FOM2cos ) /fvalue;
//  cout << "FOM = " << FOM << " " << FOM2 << endl;  


  cov->N_part = real_N_part; 	// end of work with 2. function (mixed 2.derivateves wrt D are not computed)

// computing the derivatives wrt D from func1 and wrt F,ph (those which are flagged to be computed)
  do_der_Fph = real_do_der_Fph;
  SetCentRice(0,1);
  cov->Make_matrix();
  InverseAndEigen();
  for ( nil_m = 0;  nil_m < Gauss_number;  nil_m++ )
  {
	EvaluateSR(F1,phase1);
	TransfToProb();
	if ( do_der_Fph )
	{
  	  for (int i = N_meas;  i < Num;  i++) 
  	  {
		auxder_F[i] += Gauss_weight[nil_m]*der_F[i]*func2[nil_m];
		auxder_ph[i] += Gauss_weight[nil_m]*der_ph[i]*func2[nil_m];
		if ( do_der_Fph == 2 )
  		for (int j = N_meas;  j < Num;  j++) 
  		{
		  auxder2_F2[i][j] += Gauss_weight[nil_m]*der2_F2[i][j]*func2[nil_m];
		  auxder2_phF[i][j] += Gauss_weight[nil_m]*der2_phF[i][j]*func2[nil_m];
		  auxder2_ph2[i][j] += Gauss_weight[nil_m]*der2_ph2[i][j]*func2[nil_m];
		}
	  }
	}
	if ( do_der_D )
	{
	  for (int  k = 0;  k < cov->DNum;  k++)
  		for (int  l = k+1;  l < cov->DNum;  l++)
		  for ( int pm = 0;  pm < cov->N_part;  pm++ )			// now looping over k,l,pm to get derivatives wrt D[pm,k,l]
		  {
			auxder_D[pm][k][l] += Gauss_weight[nil_m]*der_D[pm][k][l]*func2[nil_m];
			if ( do_der_D == 2 )
			{
			  for (int  p = 0;  p <= k;  p++) 
				for (int  q = p+1;  (q < cov->DNum && p < k) || ( q <= l && p == k );  q++) 
				  for (int pm2 = 0;  ( pm2 <= pm && p == k && q == l ) || ( pm2 < cov->N_part && (p < k || q < l) );  pm2++) 
					auxder2_D[pm][k][l][pm2][p][q] += Gauss_weight[nil_m]*der2_D[pm][k][l][pm2][p][q]*func2[nil_m];
			}
		  }
	}
  }
  cov->N_part = real_N_part; // end of work with 1. function
  SetCentRice(0,0);	// rice flag was originally set to 0 (otherwise Rice would have been called), centricity is assumed acentric at the moment

// assigning final values (of -log target) to der variables

  realnum Pi_fvalue_inv = 1 / fvalue;	// fvalue is divided by PI at the moment
  if ( do_der_Fph )
  {		//	carefully - the 2.function N_meas and Num has to be taken into account for der while the 1.function (N_meas=1) for auxder - this is no more true!
//	int Nmd = N_meas - 1;
  	for (int i = 1;  i < Num;  i++) 
  	{
	  der_F[i] = -auxder_F[i]*Pi_fvalue_inv;
	  der_ph[i] = -auxder_ph[i]*Pi_fvalue_inv;
	  if ( do_der_Fph == 2 )
  	  for (int j = 1;  j < Num;  j++) 
  	  {
		der2_F2[i][j] = -auxder2_F2[i][j]*Pi_fvalue_inv + auxder_F[i]*auxder_F[j]*Pi_fvalue_inv*Pi_fvalue_inv;
		der2_phF[i][j] = -auxder2_phF[i][j]*Pi_fvalue_inv + auxder_ph[i]*auxder_F[j]*Pi_fvalue_inv*Pi_fvalue_inv;
		der2_ph2[i][j] = -auxder2_ph2[i][j]*Pi_fvalue_inv + auxder_ph[i]*auxder_ph[j]*Pi_fvalue_inv*Pi_fvalue_inv;
	  }
	}
  }
  if ( do_der_D )
  {
	for (int  k = 0;  k < cov->DNum;  k++)
  	  for (int  l = k+1;  l < cov->DNum;  l++)
		for ( int pm = 0;  pm < cov->N_part;  pm++ )			// now looping over k,l,pm to get derivatives wrt D[pm,k,l]
		{
		  der_D[pm][k][l] = -auxder_D[pm][k][l]*Pi_fvalue_inv;
		  if ( do_der_D == 2 )
		  {
			for (int  p = 0;  p <= k;  p++) 
			  for (int  q = p+1;  (q < cov->DNum && p < k) || ( q <= l && p == k );  q++) 
			    for (int pm2 = 0;  ( pm2 <= pm && p == k && q == l ) || ( pm2 < cov->N_part && (p < k || q < l) );  pm2++) 
				{
				  der2_D[pm][k][l][pm2][p][q] = -auxder2_D[pm][k][l][pm2][p][q]*Pi_fvalue_inv +
												auxder_D[pm][k][l]*auxder_D[pm2][p][q]*Pi_fvalue_inv*Pi_fvalue_inv;
				  der2_D[pm2][p][q][pm][k][l] = der2_D[pm][k][l][pm2][p][q];
				}
		  }
		}
  }


  delete [] func1;
  delete [] func2;

  if ( P_F0 > 0. ) fvalue /= P_F0;
  return fval = -log(fvalue*PI) ;

  }
}

}
