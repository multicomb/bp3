#include <iomanip>

namespace multivar_llhood {
  
// handling with different D parameters schemes
struct Dassignment 			// assignes D parameters from YxY matrix of D parameters (Y<=Num) to the indices of the Num x Num matrix
							// i.e. determinates which term of cov. matrix is function of given D parameter
{
  int toNum;
  int *to_array;	// this will be 2-dimensional - Nx2 - 'to' will be used to access it (the original (*)[2] implementation was replaced by this due to problems with some compilers)
};
#define to(i,j) to_array[2*i+j]


template <typename realnum>
class likelihood;

template <typename realnum>
struct partial_models_data
{
  complex<realnum> *sigma_P;				// sigma_P values for all models
  realnum** D;								// matrix of D parameters 
  realnum **deriv_D_r, **deriv_D_i;			// real and imaginary  part of derivative of covariance matrix wrt D
};

template <typename realnum>
class covar_matrix
{
  friend class likelihood<realnum>;
  public:

	int N_part;								// number of partial models; at the moment 1 or 2 is supported
    int Num;								// total number of measurements + models
    int N_meas;								// number of measurements
	complex<realnum> *sigma_N;				// sigma_N values for all data
	realnum* sig_meas;						// measurement errors of data
	partial_models_data<realnum>* part;		// sigma_P, D param. and 1.deriv. wrt D[p,k,l] for all models and partial models
	partial_models_data<realnum>** partpart;// 'exchange' sigma_P and 2.deriv. wrt D[p1,k,l,p2,k,l] (all others are assumed 0) for all models and partial models
    int no_imag;							// if no_imag then only real parts will be taken into account (default 0)   (ovsem to neznamena zeby sa nikde imag. cleny v kode nepouzivali-tam kde je to casovo nevyhodne sa nepouzivaju, tam kde je to jedno sa nadalej pouzivaju)
//	Complex *sigma_P1;
//	realnum** D1;								// matrix of D parameters 
    covar_matrix( likelihood<realnum> *l, int N_part, int Num_p, int N_meas_p, string D_type, int do_der_D, int no_imagin );
	covar_matrix( covar_matrix<realnum>&, likelihood<realnum> *l );
	covar_matrix( likelihood<realnum> *l ) : llh(l) { re = 0; };
	covar_matrix& operator=( covar_matrix<realnum>& );
	~covar_matrix();
    void Make_matrix();						// computes re and im from sigmas and D's
	void Print();
	void SetZeroRows(int,int,int,int,int,int);	// specifies the order numbers of rows(=columns) for which the covariances should be 0.
	int CheckZeroRow(int row);				// returns 1 if row has been set to be covariance zeroed, 0 otherwise
	int GetNumZeroRows();					// returns the number of rows with covariaces zeroed

  private:

    likelihood<realnum> * llh;  
	realnum **re, **im;						// real and imaginary part of covariance matrix  
	int DNum;					// matrix of D's will be of dimenstion DNum x DNum; the diagonal terms are not assumed to be used
	int maxtoNum;
	string Dtype;				// one of { "SAD", "SIR", "SIRH", "single", ... }
	Dassignment ***Dass;		// first index is part. model number, second and third reffer to the indices of above mentioned YxY matrix (so, every partial model has different D-scheme)
	int *workspace;				// maxtoNum
	void SetAllocDtype( string &Dtype_toset );
	void Alloc( int do_der_D );
	void DeAlloc( int do_der_D );
	void DeAllocDtype( );
	inline void alloc_test(void* name, char* func_name=(char*)"") { if (!name) likelihood<realnum>::Error( 101, func_name ); }
	int zero_cov_row_num; 		// number of rows/columns for which cov terms will be set to 0.; default = 0
	int* zero_cov_row;			// the order numbers of rows to be covariance zeroed
};





//#define alloc_test(name)   if (!name) likelihood<realnum>::Error( 101, func_name )

template <typename realnum>
void covar_matrix<realnum>::DeAlloc( int do_der_D )
{
  char* func_name = (char*) "covar_matrix::DeAlloc";
  delete [] sigma_N;
  delete [] sig_meas;
  for ( int p = 0;  p < N_part;  p++ )
  {
	for ( int i = 0;  i < Num;  i++ ) 
	{
	  if ( do_der_D )
	  {
		for ( int p2 = 0;  p2 < N_part;  p2++ )
		{
		  delete [] partpart[p][p2].deriv_D_r[i];
		  delete [] partpart[p][p2].deriv_D_i[i];
		}
		delete [] part[p].deriv_D_r[i];
		delete [] part[p].deriv_D_i[i];
	  }
	  delete [] part[p].D[i];
	}
	for ( int p2 = 0;  p2 < N_part;  p2++ )
	{
	  delete [] partpart[p][p2].sigma_P;
	}
	if ( do_der_D )
	{
	  delete [] part[p].deriv_D_r;
	  delete [] part[p].deriv_D_i;
	  for ( int p2 = 0;  p2 < N_part;  p2++ )
	  {
		delete [] partpart[p][p2].deriv_D_r;
		delete [] partpart[p][p2].deriv_D_i;
	  }
	}
	delete [] part[p].D;
	delete [] part[p].sigma_P;
	delete [] partpart[p];
  }
  delete [] part;
  delete [] partpart;

  for ( int i = 0;  i < Num;  i++ )
  {
    delete [] re[i]; delete [] im[i];
  }
  delete [] re; 	delete [] im;
  delete [] zero_cov_row;
  if ( do_der_D )
  {
	delete [] workspace;
  }
};


#define alloc_test_init(var, type, dim, func) \
{ var = new type[dim]; alloc_test(var,func); memset(var,'\0',sizeof(type)*dim); }

template <typename realnum>
void covar_matrix<realnum>::Alloc( int do_der_D )
{
  char* func_name = (char*) "covar_matrix::Alloc";
  part = new partial_models_data<realnum>[N_part];					alloc_test(part);
  partpart = new partial_models_data<realnum>*[N_part];				alloc_test(partpart);
  alloc_test_init( sigma_N, complex<realnum>, Num, func_name )
//  sigma_N = new complex<realnum>[Num];								alloc_test(sigma_N);
    alloc_test_init( sig_meas, realnum, N_meas, func_name )
//  sig_meas = new realnum[N_meas];									alloc_test(sig_meas);
  for ( int p = 0;  p < N_part;  p++ )
  {
	part[p].D = new realnum*[Num];									alloc_test(part[p].D);
	
	part[p].sigma_P = new complex<realnum>[Num];					alloc_test(part[p].sigma_P);
	
	partpart[p] = new partial_models_data<realnum>[N_part];			alloc_test(partpart[p]);
	for ( int p2 = 0;  p2 < N_part;  p2++ )
	{
	  partpart[p][p2].sigma_P = new complex<realnum>[Num];			alloc_test(partpart[p][p2].sigma_P);
	}
	if ( do_der_D )
	{
	  part[p].deriv_D_r = new realnum*[Num]; 						alloc_test(part[p].deriv_D_r);
	  part[p].deriv_D_i = new realnum*[Num];						alloc_test(part[p].deriv_D_i);
	  for ( int p2 = 0;  p2 < N_part;  p2++ )
	  {
		partpart[p][p2].deriv_D_r = new realnum*[Num]; 				alloc_test(partpart[p][p2].deriv_D_r);
		partpart[p][p2].deriv_D_i = new realnum*[Num];				alloc_test(partpart[p][p2].deriv_D_i);
	  }
	}
	for ( int i = 0;  i < Num;  i++ ) 
	{
	  alloc_test_init( part[p].D[i], realnum, Num, func_name )
//	  part[p].D[i] = new realnum[Num];								alloc_test(part[p].D[i]);
	  if ( do_der_D )
	  {
	    alloc_test_init( part[p].deriv_D_r[i], realnum, Num, func_name )
//		part[p].deriv_D_r[i] = new realnum[Num];					alloc_test(part[p].deriv_D_r[i]);
	      alloc_test_init( part[p].deriv_D_i[i], realnum, Num, func_name )
//		part[p].deriv_D_i[i] = new realnum[Num];					alloc_test(part[p].deriv_D_i[i]);
		for ( int p2 = 0;  p2 < N_part;  p2++ )
		{
		  alloc_test_init( partpart[p][p2].deriv_D_r[i], realnum, Num, func_name )
//		  partpart[p][p2].deriv_D_r[i] = new realnum[Num];			alloc_test(partpart[p][p2].deriv_D_r[i]);
		    alloc_test_init( partpart[p][p2].deriv_D_i[i], realnum, Num, func_name )
//		  partpart[p][p2].deriv_D_i[i] = new realnum[Num];			alloc_test(partpart[p][p2].deriv_D_i[i]);
		}
	  }
	}
  }
  zero_cov_row = new int[Num];										alloc_test(zero_cov_row);
  re = new realnum*[Num]; 	im = new realnum*[Num];					alloc_test(re);	alloc_test(im);
  for ( int i = 0;  i < Num;  i++ )   
  {
    alloc_test_init( re[i], realnum, Num, func_name )
      alloc_test_init( im[i], realnum, Num, func_name )	
//    re[i] = new realnum[Num]; im[i] = new realnum[Num];				alloc_test(re[i]); alloc_test(im[i]);
  }
  if ( do_der_D )
  {
    alloc_test_init( workspace, int, 2*maxtoNum, func_name )
//	workspace = new int[2*maxtoNum];								alloc_test(workspace);
  }
}


template <typename realnum>
void covar_matrix<realnum>::DeAllocDtype( )
{
  if ( Dtype != "no" )
  {
	for ( int p = 0;  p < N_part;  p++ )
	{
	  for ( int i = 0;  i < DNum;  i++ )	
	  {
		for ( int j = i+1;  j < DNum;  j++ )	// working just with the upper triangle without diagonal terms
		{
		  delete [] Dass[p][i][j].to_array;
		}
		delete [] Dass[p][i];
	  }
	  delete [] Dass[p];
	}
	delete [] Dass;
  }
}


template <typename realnum>
void covar_matrix<realnum>::SetAllocDtype( string& Dtype_toset )
{
  char* func_name = (char*) "likelihood::SetAllocDtype";
//  Dtype = Dtype_toset;
//  workaround for Irix
  Dtype.insert(0,Dtype_toset);
  if ( Dtype == "no" ) Dass = 0; 
  else
	Dass = new Dassignment**[N_part];

  if ( Dtype == "no" )
  {
  }
// single is basically rice... (although can be used for other purposes)
  else if ( Dtype == "single" )
  {
	DNum = Num;
	maxtoNum = 1;
	if (N_part>2) maxtoNum = 2;
	for ( int p = 0;  p < N_part;  p++ ) 
	{
	  Dass[p] = new Dassignment*[DNum];
	  for ( int i = 0;  i < DNum;  i++ ) Dass[p][i] = new Dassignment[DNum];

	  for ( int i = 0;  i < DNum;  i++ )	
	  {
		for ( int j = i+1;  j < DNum;  j++ )	// working just with the upper triangle without diagonal terms
		{
		  Dass[p][i][j].toNum = 1;
		  if (p==2) Dass[p][i][j].toNum = 2;
		  if (p==1||p==3) Dass[p][i][j].toNum = 0;
		  Dass[p][i][j].to_array = new int[ 2*Dass[p][i][j].toNum ];
		  if (p!=1&&p!=3) Dass[p][i][j].to(0,0) = i;
		  if (p!=1&&p!=3) Dass[p][i][j].to(0,1) = j;
		  if (p==2) Dass[p][i][j].to(1,0) = 1;
		  if (p==2) Dass[p][i][j].to(1,1) = 1;
		}
	  }
	}
  }
  else if ( Dtype == "p+l" )
  {
	if ( Num == 4 )
	// D parameters not yet defined in the refmac way!!! This stands for all "p+l" num=4 cases.
	{
	  DNum = Num/2;
	  maxtoNum = 4;
	  for ( int p = 0; p < N_part; p++ )
	  {
		Dass[p] = new Dassignment*[DNum];
		for ( int i = 0;  i < DNum;  i++ ) Dass[p][i] = new Dassignment[DNum];
		if (p==0)       Dass[p][0][1].toNum = 1;
		if (p==1||p==2) Dass[p][0][1].toNum = 1;
		if (p==3)       Dass[p][0][1].toNum = 1;
		Dass[p][0][1].to_array = new int[ 2*Dass[p][0][1].toNum ];
		
		if ( p == 0 )
		{
		  Dass[p][0][1].to(0,0) = 0;		Dass[p][0][1].to(0,1) = 2;
		//  Dass[p][0][1].to(1,0) = 1;		Dass[p][0][1].to(1,1) = 2;
		//  Dass[p][0][1].to(2,0) = 0;		Dass[p][0][1].to(2,1) = 3;
		//  Dass[p][0][1].to(3,0) = 1;		Dass[p][0][1].to(3,1) = 3;
		}
		else if ( p == 1 )
		{
		  Dass[p][0][1].to(0,0) = 0;		Dass[p][0][1].to(0,1) = 1;
		//  Dass[p][0][1].to(1,0) = 1;		Dass[p][0][1].to(1,1) = 2;
		}
		else if ( p == 2 )
		{
//!!!!!!!1
//		  Dass[p][0][1].to(0,0) = 0;		Dass[p][0][1].to(0,1) = 3;
//		  Dass[p][0][1].to(1,0) = 2;		Dass[p][0][1].to(1,1) = 3;
		  Dass[p][0][1].to(0,0) = 1;		Dass[p][0][1].to(0,1) = 2;
		}
		else if ( p == 3 )
		{
		  Dass[p][0][1].to(0,0) = 1;		Dass[p][0][1].to(0,1) = 3;
		}
		else likelihood<realnum>::Error(106, func_name, (void*) &Num );
	  }
	}
    else if ( Num==3 )
    // D parameters defined in the refmac way (before D2 was D1 and vice versa) This stands for all "p+l" num=3 cases.
    {
  	  DNum = 2 ;
  	  maxtoNum = 2;
  	  for ( int p = 0;  p < N_part;  p++ )
  	  {
		Dass[p] = new Dassignment*[DNum];
		for ( int i = 0;  i < DNum;  i++ ) Dass[p][i] = new Dassignment[DNum];
		Dass[p][0][1].toNum = 1;
		Dass[p][0][1].to_array = new int[ 2*Dass[p][0][1].toNum ];
		if ( p == 0 )
		{
		//  Dass[p][0][1].to(0,0) = 0;		
		//  Dass[p][0][1].to(0,1) = 2;
		  Dass[p][0][1].to(0,0) = 1;		Dass[p][0][1].to(0,1) = 2;
		}
		else if ( p == 1 )
		{
		  Dass[p][0][1].to(0,0) = 0;		Dass[p][0][1].to(0,1) = 1;
//!!!!!
		//  Dass[p][0][1].to(1,0) = 0;
		//  Dass[p][0][1].to(1,1) = 0;
		}
//!!!!!!!!
		else if ( p == 2 )
		{
		  Dass[p][0][1].to(0,0) = 0;		Dass[p][0][1].to(0,1) = 2;
		//  Dass[p][0][1].to(1,0) = 0;
		//  Dass[p][0][1].to(1,1) = 0;
		}
  	  }  
  	}
	else likelihood<realnum>::Error(106, func_name, (void*) &Num ); 
  }	
  else if ( Dtype == "sad" || Dtype == "sadh" ) 
  {
	if ( Num%2 == 0 ) 
	{
	  DNum = Num/2;
	  maxtoNum = 4;
	  if (N_part>2) maxtoNum = 7;
	  for ( int p = 0;  p < N_part;  p++ ) 
	  {
		Dass[p] = new Dassignment*[DNum];
		for ( int i = 0;  i < DNum;  i++ ) Dass[p][i] = new Dassignment[DNum];

		for ( int i = 0;  i < DNum;  i++ )	
		{
		  for ( int j = i+1;  j < DNum;  j++ )	// working just with the upper triangle
		  {
			Dass[p][i][j].toNum = 4;
			if (p==2) Dass[p][i][j].toNum = 7;
			if (p==3) Dass[p][i][j].toNum = 1;
// these are for Raj's derivatives
			if (p==4) Dass[p][i][j].toNum = 3;
			if (p==5) Dass[p][i][j].toNum = 4;
			Dass[p][i][j].to_array = new int[ 2*Dass[p][i][j].toNum ];
		  
			for ( int k = 0;  k < Dass[p][i][j].toNum;  k++ )
			{
			  if (p!=3)	  Dass[p][i][j].to(k,0) = 2*i + k/2;
			  if (p!=3)	  Dass[p][i][j].to(k,1) = 2*j + k%2;
//			cout << i << " " << j << " " << k << "     " << Dass[i][j].to[k][0] << " " << Dass[i][j].to[k][1] << endl;
			}
			if (p==2) Dass[p][i][j].to(4,0) = 3;
			if (p==2) Dass[p][i][j].to(4,1) = 3;
			if (p==2) Dass[p][i][j].to(5,0) = 2;
			if (p==2) Dass[p][i][j].to(5,1) = 2;
			if (p==2) Dass[p][i][j].to(6,0) = 0;
			if (p==2) Dass[p][i][j].to(6,1) = 1;			
			if (p==3) Dass[p][i][j].to(0,0) = 0;
			if (p==3) Dass[p][i][j].to(0,1) = 1;
// these are for Raj's derivatives
			if (p==4) Dass[p][i][j].to(0,0) = 0;		if (p==4) Dass[p][i][j].to(0,1) = 3;
			if (p==4) Dass[p][i][j].to(1,0) = 1;		if (p==4) Dass[p][i][j].to(1,1) = 3;
			if (p==4) Dass[p][i][j].to(2,0) = 3;		if (p==4) Dass[p][i][j].to(2,1) = 3;
			if (p==5) Dass[p][i][j].to(0,0) = 0;		if (p==5) Dass[p][i][j].to(0,1) = 2;
			if (p==5) Dass[p][i][j].to(1,0) = 1;		if (p==5) Dass[p][i][j].to(1,1) = 2;
			if (p==5) Dass[p][i][j].to(2,0) = 2;		if (p==5) Dass[p][i][j].to(2,1) = 2;
			if (p==5) Dass[p][i][j].to(3,0) = 0;		if (p==5) Dass[p][i][j].to(3,1) = 1;
		  }
		}
	  }
	}
	else if ( Num == 3 )
	{
	  DNum = 2;
	  maxtoNum = 2;
	  for ( int p = 0;  p < N_part;  p++ ) 
	  {
		Dass[p] = new Dassignment*[DNum];
		for ( int i = 0;  i < DNum;  i++ ) Dass[p][i] = new Dassignment[DNum];
		Dass[p][0][1].toNum = 2;
		Dass[p][0][1].to_array = new int[ 2*Dass[p][0][1].toNum ];
		Dass[p][0][1].to(0,0) = 0;
		Dass[p][0][1].to(0,1) = 2;
		Dass[p][0][1].to(1,0) = 1;
		Dass[p][0][1].to(1,1) = 2;
	  }
	}
	else if ( Num == 5 ) // with an additional (eg density modified) real space Fc
	{
	  DNum = 2;
	  maxtoNum = 4;
//!!!!11	  if (N_part>2) maxtoNum = 7;
	  if (N_part>2) maxtoNum = 11;
	  for ( int p = 0;  p < N_part;  p++ ) 
	  {
		Dass[p] = new Dassignment*[DNum];
		for ( int i = 0;  i < DNum;  i++ ) Dass[p][i] = new Dassignment[DNum];

		for ( int i = 0;  i < DNum;  i++ )	
		{
		  for ( int j = i+1;  j < DNum;  j++ )	// working just with the upper triangle
		  {
			Dass[p][i][j].toNum = 4;
			if (p==0) Dass[p][i][j].toNum = 6; 
/*!!!!!!!*/			if (p==2) Dass[p][i][j].toNum = 7; 
			if (p==2) Dass[p][i][j].toNum = 11; 
			if (p==3) Dass[p][i][j].toNum = 1; 
			if (p==5) Dass[p][i][j].toNum = 2;
//!!!!!!!!!!!!
			if (p==1) Dass[p][i][j].toNum = 6;
// these are for Raj's derivatives
			if (p==6) Dass[p][i][j].toNum = 4;
			if (p==7) Dass[p][i][j].toNum = 4;
///*!!!!!*/			if (p==4/*||p==5|| p==3*/) Dass[p][i][j].toNum = 0;
			Dass[p][i][j].to_array = new int[ 2*Dass[p][i][j].toNum ];
		  
			for ( int k = 0;  k < Dass[p][i][j].toNum;  k++ )
			{
			  if (p<3)	  Dass[p][i][j].to(k,0) = 2*i + k/2;
			  if (p<3)	  Dass[p][i][j].to(k,1) = 2*j + k%2;
//			cout << i << " " << j << " " << k << "     " << Dass[i][j].to[k][0] << " " << Dass[i][j].to[k][1] << endl;
			}
			if (p==0) Dass[p][i][j].to(4,0) = 0;
			if (p==0) Dass[p][i][j].to(4,1) = 4;
			if (p==0) Dass[p][i][j].to(5,0) = 1;
			if (p==0) Dass[p][i][j].to(5,1) = 4;
		  // "standard" scaling
			if (p==2) Dass[p][i][j].to(4,0) = 3;
			if (p==2) Dass[p][i][j].to(4,1) = 3;
			if (p==2) Dass[p][i][j].to(5,0) = 2;
			if (p==2) Dass[p][i][j].to(5,1) = 2;
			if (p==2) Dass[p][i][j].to(6,0) = 0;
			if (p==2) Dass[p][i][j].to(6,1) = 1;			
//!!!!!!!!!!!!
			if (p==2) Dass[p][i][j].to(7,0) = 4;
			if (p==2) Dass[p][i][j].to(7,1) = 4;
			if (p==2) Dass[p][i][j].to(8,0) = 0;
			if (p==2) Dass[p][i][j].to(8,1) = 4;
			if (p==2) Dass[p][i][j].to(9,0) = 1;
			if (p==2) Dass[p][i][j].to(9,1) = 4;
			if (p==2) Dass[p][i][j].to(10,0) = 3;
			if (p==2) Dass[p][i][j].to(10,1) = 4;
		  // D between standard and additional Fc
			if (p==3) Dass[p][i][j].to(0,0) = 3;
			if (p==3) Dass[p][i][j].to(0,1) = 4;
		  // additional model scaling
			if (p==4) Dass[p][i][j].to(0,0) = 0; 
			if (p==4) Dass[p][i][j].to(0,1) = 3; 
			if (p==4) Dass[p][i][j].to(1,0) = 1; 
			if (p==4) Dass[p][i][j].to(1,1) = 3; 
			if (p==4) Dass[p][i][j].to(2,0) = 2; 
			if (p==4) Dass[p][i][j].to(2,1) = 3; 
			if (p==4) Dass[p][i][j].to(3,0) = 3; 
			if (p==4) Dass[p][i][j].to(3,1) = 3; 
		   // between observed and (not anymore additional) calculated
			if (p==5) Dass[p][i][j].to(0,0) = 0; 
			if (p==5) Dass[p][i][j].to(0,1) = 4; 
			if (p==5) Dass[p][i][j].to(1,0) = 1; 
			if (p==5) Dass[p][i][j].to(1,1) = 4; 
//!!!!!!!!!!!
			if (p==1) Dass[p][i][j].to(4,0) = 0; 
			if (p==1) Dass[p][i][j].to(4,1) = 4; 
			if (p==1) Dass[p][i][j].to(5,0) = 1; 
			if (p==1) Dass[p][i][j].to(5,1) = 4; 
// these are for Raj's derivatives
			if (p==6) Dass[p][i][j].to(0,0) = 0;		if (p==6) Dass[p][i][j].to(0,1) = 4;
			if (p==6) Dass[p][i][j].to(1,0) = 1;		if (p==6) Dass[p][i][j].to(1,1) = 4;
			if (p==6) Dass[p][i][j].to(2,0) = 3;		if (p==6) Dass[p][i][j].to(2,1) = 4;
			if (p==6) Dass[p][i][j].to(3,0) = 4;		if (p==6) Dass[p][i][j].to(3,1) = 4;
			if (p==7) Dass[p][i][j].to(0,0) = 0;		if (p==7) Dass[p][i][j].to(0,1) = 2;
			if (p==7) Dass[p][i][j].to(1,0) = 1;		if (p==7) Dass[p][i][j].to(1,1) = 2;
			if (p==7) Dass[p][i][j].to(2,0) = 2;		if (p==7) Dass[p][i][j].to(2,1) = 2;
			if (p==7) Dass[p][i][j].to(3,0) = 0;		if (p==7) Dass[p][i][j].to(3,1) = 1;
		  }
		}
	  }
	}
	else likelihood<realnum>::Error(106, func_name, (void*) &Num );
  }
  else if ( Dtype == "sadh2d" ) 
  {
	if ( Num == 4 ) 
	{
	  DNum = Num/2;
	  maxtoNum = 2;
	  for ( int p = 0;  p < N_part;  p++ ) 
	  {
		Dass[p] = new Dassignment*[DNum];
		for ( int i = 0;  i < DNum;  i++ ) Dass[p][i] = new Dassignment[DNum];

		for ( int i = 0;  i < DNum;  i++ )	
		{
		  for ( int j = i+1;  j < DNum;  j++ )	// working just with the upper triangle
		  {
			Dass[p][i][j].toNum = 2;
			Dass[p][i][j].to_array = new int[ 2*Dass[p][i][j].toNum ];
		  
			for ( int k = 0;  k < Dass[p][i][j].toNum;  k++ )
			{
			  if ( p == 0 ) {
				Dass[p][i][j].to(k,0) = k;
				Dass[p][i][j].to(k,1) = 3;
			  } 
			  else if ( p == 1 ) {
				Dass[p][i][j].to(k,0) = k;
				Dass[p][i][j].to(k,1) = 2;
			  }
//			cout << i << " " << j << " " << k << "     " << Dass[i][j].to[k][0] << " " << Dass[i][j].to[k][1] << endl;
			}
		  }
		}
	  }
	}
	else likelihood<realnum>::Error(106, func_name, (void*) &Num );
  }
  else if ( Dtype == "sir" ) 
  {
	if ( Num == 4 ) 
	{
	  DNum = 2;
	  maxtoNum = 7;
	  for ( int p = 0;  p < N_part;  p++ ) 
	  {
		Dass[p] = new Dassignment*[DNum];
		for ( int i = 0;  i < DNum;  i++ ) Dass[p][i] = new Dassignment[DNum];

		// working just with the upper triangle
		Dass[p][0][1].toNum = 4;
		if (p==2) Dass[p][0][1].toNum = 7;
		Dass[p][0][1].to_array = new int[ 2*Dass[p][0][1].toNum ] ;
		  
//		for ( int k = 0;  k < Dass[p][0][1].toNum;  k++ )
//		{
		  if ( p == 0 )
		  {
			Dass[p][0][1].to(0,0) = 0;  	Dass[p][0][1].to(0,1) = 2;
			Dass[p][0][1].to(1,0) = 0;  	Dass[p][0][1].to(1,1) = 3;
			Dass[p][0][1].to(2,0) = 1;  	Dass[p][0][1].to(2,1) = 2;
			Dass[p][0][1].to(3,0) = 1;  	Dass[p][0][1].to(3,1) = 3;
//			Dass[p][0][1].to(k,0) = k/2;
//			Dass[p][0][1].to(k,1) = 2 + k%2;
		  }
		  else if ( p == 1 )
		  {
			Dass[p][0][1].to(0,0) = 0;  	Dass[p][0][1].to(0,1) = 1;
			Dass[p][0][1].to(1,0) = 1;  	Dass[p][0][1].to(1,1) = 2;
			Dass[p][0][1].to(2,0) = 2;  	Dass[p][0][1].to(2,1) = 3;
			Dass[p][0][1].to(3,0) = 0;  	Dass[p][0][1].to(3,1) = 3;
//			Dass[p][0][1].to(k,0) = k%3;
//			Dass[p][0][1].to(k,1) = k%3 + 1 + (k/3)*2;
		  }
		  else if ( p == 2 )
		  {
			Dass[p][0][1].to(0,0) = 0;		Dass[p][0][1].to(0,1) = 2;
			Dass[p][0][1].to(1,0) = 0;		Dass[p][0][1].to(1,1) = 3;
			Dass[p][0][1].to(2,0) = 1;		Dass[p][0][1].to(2,1) = 2;
			Dass[p][0][1].to(3,0) = 1;		Dass[p][0][1].to(3,1) = 3;
			Dass[p][0][1].to(4,0) = 2;		Dass[p][0][1].to(4,1) = 2;
			Dass[p][0][1].to(5,0) = 2;		Dass[p][0][1].to(5,1) = 3;
			Dass[p][0][1].to(6,0) = 3;		Dass[p][0][1].to(6,1) = 3;
		  }
		  else likelihood<realnum>::Error( 111, func_name );
//			cout << i << " " << j << " " << k << "     " << Dass[i][j].to[k][0] << " " << Dass[i][j].to[k][1] << endl;
//		}
	  }
	}
  }
  else if ( Dtype == "sirh" ) 
  {
	if ( Num == 4 ) 
	{
	  DNum = 2;
	  maxtoNum = 7;
	  for ( int p = 0;  p < N_part;  p++ ) 
	  {
		Dass[p] = new Dassignment*[DNum];
		for ( int i = 0;  i < DNum;  i++ ) Dass[p][i] = new Dassignment[DNum];

		// working just with the upper triangle
		Dass[p][0][1].toNum = 4;
		if (p==2) Dass[p][0][1].toNum = 7;
		if (p==3||p==4) Dass[p][0][1].toNum = 2;
		Dass[p][0][1].to_array = new int[ 2*Dass[p][0][1].toNum ] ;
		  
//		for ( int k = 0;  k < Dass[p][0][1].toNum;  k++ )
//		{
		  if ( p == 0 )
		  {
			Dass[p][0][1].to(0,0) = 0;  	Dass[p][0][1].to(0,1) = 2;
			Dass[p][0][1].to(1,0) = 0;  	Dass[p][0][1].to(1,1) = 3;
			Dass[p][0][1].to(2,0) = 1;  	Dass[p][0][1].to(2,1) = 2;
			Dass[p][0][1].to(3,0) = 1;  	Dass[p][0][1].to(3,1) = 3;
//			Dass[p][0][1].to(k,0) = k/2;
//			Dass[p][0][1].to(k,1) = 2 + k%2;
		  }
		  else if ( p == 1 )
		  {
			Dass[p][0][1].to(0,0) = 0;  	Dass[p][0][1].to(0,1) = 1;
			Dass[p][0][1].to(1,0) = 1;  	Dass[p][0][1].to(1,1) = 2;
			Dass[p][0][1].to(2,0) = 2;  	Dass[p][0][1].to(2,1) = 3;
			Dass[p][0][1].to(3,0) = 0;  	Dass[p][0][1].to(3,1) = 3;
//			Dass[p][0][1].to(k,0) = k%3;
//			Dass[p][0][1].to(k,1) = k%3 + 1 + (k/3)*2;
		  }
		  else if ( p == 2 )
		  {
			Dass[p][0][1].to(0,0) = 0;		Dass[p][0][1].to(0,1) = 2;
			Dass[p][0][1].to(1,0) = 0;		Dass[p][0][1].to(1,1) = 3;
			Dass[p][0][1].to(2,0) = 1;		Dass[p][0][1].to(2,1) = 2;
			Dass[p][0][1].to(3,0) = 1;		Dass[p][0][1].to(3,1) = 3;
			Dass[p][0][1].to(4,0) = 2;		Dass[p][0][1].to(4,1) = 2;
			Dass[p][0][1].to(5,0) = 2;		Dass[p][0][1].to(5,1) = 3;
			Dass[p][0][1].to(6,0) = 3;		Dass[p][0][1].to(6,1) = 3;
		  }
		  else if ( p == 3 )
		  {
			Dass[p][0][1].to(0,0) = 1;		Dass[p][0][1].to(0,1) = 3;
//			Dass[p][0][1].to(0,0) = 3;		Dass[p][0][1].to(0,1) = 3;
//			Dass[p][0][1].to(0,0) = 2;		Dass[p][0][1].to(0,1) = 3;
		  }
		  else if ( p == 4 )
		  {
			Dass[p][0][1].to(0,0) = 1;		Dass[p][0][1].to(0,1) = 2;
		  }
		  else likelihood<realnum>::Error( 111, func_name );
//		}
	  }
	}
  }
  else if ( Dtype == "sirph" ) 
  {
	if ( Num == 3 ) 
	{
	  DNum = 2;
	  maxtoNum = 3;
	  for ( int p = 0;  p < N_part;  p++ ) 
	  {
		Dass[p] = new Dassignment*[DNum];
		for ( int i = 0;  i < DNum;  i++ ) Dass[p][i] = new Dassignment[DNum];

		// working just with the upper triangle
		if (p==0||p==1) Dass[p][0][1].toNum = 1;
		if (p==2) 		Dass[p][0][1].toNum = 3;
		Dass[p][0][1].to_array = new int[ 2*Dass[p][0][1].toNum ] ;
		  
		if ( p == 0 ) {
		  Dass[p][0][1].to(0,0) = 1;  	Dass[p][0][1].to(0,1) = 2;
		}
		else if ( p == 1 ) {
		  Dass[p][0][1].to(0,0) = 0;  	Dass[p][0][1].to(0,1) = 1;
		}
		else if ( p == 2 ) {
		  Dass[p][0][1].to(0,0) = 0;  	Dass[p][0][1].to(0,1) = 2;
		  Dass[p][0][1].to(1,0) = 1;  	Dass[p][0][1].to(1,1) = 2;
		  Dass[p][0][1].to(2,0) = 2;  	Dass[p][0][1].to(2,1) = 2;
		}
		else likelihood<realnum>::Error( 111, func_name );
	  }
	}
  }
  else if ( Dtype == "mldr" ) 
  {
	if ( Num == 3 )
	{
	  DNum = 2;
	  maxtoNum = 1;
	  for ( int p = 0;  p < N_part;  p++ )
	  {
		Dass[p] = new Dassignment*[DNum];
		for ( int i = 0;  i < DNum;  i++ ) Dass[p][i] = new Dassignment[DNum];

		// working just with the upper triangle
		Dass[p][0][1].toNum = 1;
		Dass[p][0][1].to_array = new int[ 2*Dass[p][0][1].toNum ] ;
		  
		for ( int k = 0;  k < Dass[p][0][1].toNum;  k++ )
		{
		  if ( p == 0 )
		  {
			Dass[p][0][1].to(k,0) = 0;
			Dass[p][0][1].to(k,1) = 1;
		  }
		  else if ( p == 1 )
		  {
			Dass[p][0][1].to(k,0) = 0;
			Dass[p][0][1].to(k,1) = 1;
		  }
		  else likelihood<realnum>::Error( 111, func_name );
//			cout << i << " " << j << " " << k << "     " << Dass[i][j].to[k][0] << " " << Dass[i][j].to[k][1] << endl;
		}
	  }
	}
  }
// SIRAS with 2 D parameters with models of F_N,F_D+,F_D-
  else if ( Dtype == "sras" ) 
  {
	if ( Num == 6 ) 
	{
	  DNum = 2; // however, this only means one D parameter - the second one is added as "partial" (the same as in SIR case)
	  maxtoNum = 9;
	  for ( int p = 0;  p < N_part;  p++ ) 
	  {
		Dass[p] = new Dassignment*[DNum];
		for ( int i = 0;  i < DNum;  i++ ) Dass[p][i] = new Dassignment[DNum];

		// working just with the upper triangle
		if ( p == 0 ) Dass[p][0][1].toNum = 9;
		else if ( p == 1 ) Dass[p][0][1].toNum = 8;
		else likelihood<realnum>::Error( 111, func_name );

		Dass[p][0][1].to_array = new int[ 2*Dass[p][0][1].toNum ] ;
		  
		if ( p == 0 )
		{
		  Dass[p][0][1].to(0,0) = 0;		Dass[p][0][1].to(0,1) = 3;
		  Dass[p][0][1].to(1,0) = 0;		Dass[p][0][1].to(1,1) = 4;
		  Dass[p][0][1].to(2,0) = 0;		Dass[p][0][1].to(2,1) = 5;
		  Dass[p][0][1].to(3,0) = 1;		Dass[p][0][1].to(3,1) = 3;
		  Dass[p][0][1].to(4,0) = 1;		Dass[p][0][1].to(4,1) = 4;
		  Dass[p][0][1].to(5,0) = 1;		Dass[p][0][1].to(5,1) = 5;
		  Dass[p][0][1].to(6,0) = 2;		Dass[p][0][1].to(6,1) = 3;
		  Dass[p][0][1].to(7,0) = 2;		Dass[p][0][1].to(7,1) = 4;
		  Dass[p][0][1].to(8,0) = 2;		Dass[p][0][1].to(8,1) = 5;
		}
		else if ( p == 1 )
		{
		  Dass[p][0][1].to(0,0) = 0;		Dass[p][0][1].to(0,1) = 1;
		  Dass[p][0][1].to(1,0) = 0;		Dass[p][0][1].to(1,1) = 2;
		  Dass[p][0][1].to(2,0) = 0;		Dass[p][0][1].to(2,1) = 4;
		  Dass[p][0][1].to(3,0) = 0;		Dass[p][0][1].to(3,1) = 5;
		  Dass[p][0][1].to(4,0) = 1;		Dass[p][0][1].to(4,1) = 3;
		  Dass[p][0][1].to(5,0) = 2;		Dass[p][0][1].to(5,1) = 3;
		  Dass[p][0][1].to(6,0) = 3;		Dass[p][0][1].to(6,1) = 4;
		  Dass[p][0][1].to(7,0) = 3;		Dass[p][0][1].to(7,1) = 5;
		}
	  }
	} else likelihood<realnum>::Error(106, func_name, (void*) &Num );
  }
// SIRAS with 3 D parameters with models of F_N,F_D+,F_D-
  else if ( Dtype == "sras3d" )
  {
	if ( Num == 6 ) 
	{
	  DNum = 2; // however, this only means one D parameter - the second one is added as "partial" (the same as in SIR case)
				// it could be changed to "right" scheme  for siras3d but left the partial model's representation for compatibility with siras (and easier mantainability)
	  maxtoNum = 5;
	  for ( int p = 0;  p < N_part;  p++ ) 
	  {
		Dass[p] = new Dassignment*[DNum];
		for ( int i = 0;  i < DNum;  i++ ) Dass[p][i] = new Dassignment[DNum];

		// working just with the upper triangle
		if ( p == 0 ) Dass[p][0][1].toNum = 5;
		else if ( p == 1 ) Dass[p][0][1].toNum = 4;
		else if ( p == 2 ) Dass[p][0][1].toNum = 4;
		else likelihood<realnum>::Error( 111, func_name );

		Dass[p][0][1].to_array = new int[ 2*Dass[p][0][1].toNum ] ;
		  
		if ( p == 0 )
		{
		  Dass[p][0][1].to(0,0) = 0;		Dass[p][0][1].to(0,1) = 3;
		  Dass[p][0][1].to(1,0) = 1;		Dass[p][0][1].to(1,1) = 4;
		  Dass[p][0][1].to(2,0) = 1;		Dass[p][0][1].to(2,1) = 5;
		  Dass[p][0][1].to(3,0) = 2;		Dass[p][0][1].to(3,1) = 4;
		  Dass[p][0][1].to(4,0) = 2;		Dass[p][0][1].to(4,1) = 5;
		}
		else if ( p == 1 )
		{
		  Dass[p][0][1].to(0,0) = 0;		Dass[p][0][1].to(0,1) = 1;
		  Dass[p][0][1].to(1,0) = 0;		Dass[p][0][1].to(1,1) = 2;
		  Dass[p][0][1].to(2,0) = 3;		Dass[p][0][1].to(2,1) = 4;
		  Dass[p][0][1].to(3,0) = 3;		Dass[p][0][1].to(3,1) = 5;
		}
		else if ( p == 2 )
		{
		  Dass[p][0][1].to(0,0) = 0;		Dass[p][0][1].to(0,1) = 4;
		  Dass[p][0][1].to(1,0) = 0;		Dass[p][0][1].to(1,1) = 5;
		  Dass[p][0][1].to(2,0) = 1;		Dass[p][0][1].to(2,1) = 3;
		  Dass[p][0][1].to(3,0) = 2;		Dass[p][0][1].to(3,1) = 3;
		}
	  }
	} else likelihood<realnum>::Error(106, func_name, (void*) &Num );
  }
// SIRAS with 4 D parameters with models of F_N, F_D (non-anom), F_Danom
  else if ( Dtype == "sras4d" )
  {
	if ( Num == 6 ) 
	{
	  DNum = 2; // however, this only means one D parameter - the second one is added as "partial" (the same as in SIR case)
				// it could be changed to "right" scheme  for siras4d but left the partial model's representation for compatibility with siras (and easier mantainability)
	  maxtoNum = 13;
	  for ( int p = 0;  p < N_part;  p++ ) 
	  {
		Dass[p] = new Dassignment*[DNum];
		for ( int i = 0;  i < DNum;  i++ ) Dass[p][i] = new Dassignment[DNum];

		// working just with the upper triangle
		if ( p == 0 ) Dass[p][0][1].toNum = 5;
		else if ( p == 1 ) Dass[p][0][1].toNum = 3;
		else if ( p == 4 ) Dass[p][0][1].toNum = 3;
//!!!!!!!!!!!!!!!
//		else if ( p == 3 ) Dass[p][0][1].toNum = 2;
		else if ( p == 3 ) Dass[p][0][1].toNum = 3;
		else if ( p == 2 ) Dass[p][0][1].toNum = 13;
		else if ( p == 5 ) Dass[p][0][1].toNum = 2;
		else likelihood<realnum>::Error( 111, func_name );

		Dass[p][0][1].to_array = new int[ 2*Dass[p][0][1].toNum ] ;
		  
		if ( p == 0 )
		{
		  Dass[p][0][1].to(0,0) = 0;		Dass[p][0][1].to(0,1) = 3;
		  Dass[p][0][1].to(1,0) = 1;		Dass[p][0][1].to(1,1) = 4;
		  Dass[p][0][1].to(2,0) = 2;		Dass[p][0][1].to(2,1) = 4;
		  Dass[p][0][1].to(3,0) = 1;		Dass[p][0][1].to(3,1) = 5;
		  Dass[p][0][1].to(4,0) = 2;		Dass[p][0][1].to(4,1) = 5;
		}
		else if ( p == 1 )
		{
		  Dass[p][0][1].to(0,0) = 0;		Dass[p][0][1].to(0,1) = 1;
		  Dass[p][0][1].to(1,0) = 0;		Dass[p][0][1].to(1,1) = 2;
		  Dass[p][0][1].to(2,0) = 3;		Dass[p][0][1].to(2,1) = 4;
		}
		else if ( p == 4 )
		{
		  Dass[p][0][1].to(0,0) = 0;		Dass[p][0][1].to(0,1) = 4;
		  Dass[p][0][1].to(1,0) = 1;		Dass[p][0][1].to(1,1) = 3;
		  Dass[p][0][1].to(2,0) = 2;		Dass[p][0][1].to(2,1) = 3;
		}
		else if ( p == 3 )
		{
		  Dass[p][0][1].to(0,0) = 1;		Dass[p][0][1].to(0,1) = 5;
		  Dass[p][0][1].to(1,0) = 2;		Dass[p][0][1].to(1,1) = 5;
//!!!!!!
		  Dass[p][0][1].to(2,0) = 0;		Dass[p][0][1].to(2,1) = 3;
        }
		else if ( p == 5 )
		{
		  Dass[p][0][1].to(0,0) = 1;		Dass[p][0][1].to(0,1) = 4;
		  Dass[p][0][1].to(0,0) = 2;		Dass[p][0][1].to(0,1) = 4;
		}
		else if ( p == 2 )
		{
		  Dass[p][0][1].to(0,0) = 0;		Dass[p][0][1].to(0,1) = 3;
		  Dass[p][0][1].to(1,0) = 0;		Dass[p][0][1].to(1,1) = 4;
		  Dass[p][0][1].to(2,0) = 1;		Dass[p][0][1].to(2,1) = 2;
		  Dass[p][0][1].to(3,0) = 1;		Dass[p][0][1].to(3,1) = 3;
		  Dass[p][0][1].to(4,0) = 1;		Dass[p][0][1].to(4,1) = 4;
		  Dass[p][0][1].to(5,0) = 1;		Dass[p][0][1].to(5,1) = 5;
		  Dass[p][0][1].to(6,0) = 2;		Dass[p][0][1].to(6,1) = 3;
		  Dass[p][0][1].to(7,0) = 2;		Dass[p][0][1].to(7,1) = 4;
		  Dass[p][0][1].to(8,0) = 2;		Dass[p][0][1].to(8,1) = 5;
		  Dass[p][0][1].to(9,0) = 3;		Dass[p][0][1].to(9,1) = 3;
		  Dass[p][0][1].to(10,0) = 3;		Dass[p][0][1].to(10,1) = 4;
		  Dass[p][0][1].to(11,0) = 4;		Dass[p][0][1].to(11,1) = 4;
		  Dass[p][0][1].to(12,0) = 5;		Dass[p][0][1].to(12,1) = 5;
		}
	  }
	} else likelihood<realnum>::Error(106, func_name, (void*) &Num );
  }
// the same as sras4d but without native model ie 5x5 for scaling purposes
  else if ( Dtype == "srasph" )
  {
	if ( Num == 5 ) 
	{
	  DNum = 2; // however, this only means one D parameter - the second one is added as "partial" (the same as in SIR case)
				// it could be changed to "right" scheme  for siras4d but left the partial model's representation for compatibility with siras (and easier mantainability)
	  maxtoNum = 8;
	  for ( int p = 0;  p < N_part;  p++ ) 
	  {
		Dass[p] = new Dassignment*[DNum];
		for ( int i = 0;  i < DNum;  i++ ) Dass[p][i] = new Dassignment[DNum];

		// working just with the upper triangle
		if ( p == 0 ) Dass[p][0][1].toNum = 5;
		else if ( p == 1 ) Dass[p][0][1].toNum = 3;
		else if ( p == 3 ) Dass[p][0][1].toNum = 2;
		else if ( p == 2 ) Dass[p][0][1].toNum = 8;
//		else if ( p == 4 ) Dass[p][0][1].toNum = 1;
// these are for Raj's derivatives
		else if ( p == 4 ) Dass[p][0][1].toNum = 3;
		else if ( p == 5 ) Dass[p][0][1].toNum = 4;
		else likelihood<realnum>::Error( 111, func_name );

		Dass[p][0][1].to_array = new int[ 2*Dass[p][0][1].toNum ] ;
		  
		if ( p == 0 )
		{
		  Dass[p][0][1].to(0,0) = 1;		Dass[p][0][1].to(0,1) = 3;
		  Dass[p][0][1].to(1,0) = 2;		Dass[p][0][1].to(1,1) = 3;
		  Dass[p][0][1].to(2,0) = 1;		Dass[p][0][1].to(2,1) = 4;
		  Dass[p][0][1].to(3,0) = 2;		Dass[p][0][1].to(3,1) = 4;
		  Dass[p][0][1].to(4,0) = 0;		Dass[p][0][1].to(4,1) = 3;
		}
		else if ( p == 1 )
		{
		  Dass[p][0][1].to(0,0) = 0;		Dass[p][0][1].to(0,1) = 1;
		  Dass[p][0][1].to(1,0) = 0;		Dass[p][0][1].to(1,1) = 2;
		  Dass[p][0][1].to(2,0) = 0;		Dass[p][0][1].to(2,1) = 3;
		}
		else if ( p == 3 )
		{
		  Dass[p][0][1].to(0,0) = 1;		Dass[p][0][1].to(0,1) = 4;
		  Dass[p][0][1].to(1,0) = 2;		Dass[p][0][1].to(1,1) = 4;
//		  Dass[p][0][1].to(1,0) = 4;		Dass[p][0][1].to(1,1) = 4;
        }
		else if ( p == 2 )
		{
		  Dass[p][0][1].to(0,0) = 1;		Dass[p][0][1].to(0,1) = 2;
		  Dass[p][0][1].to(1,0) = 0;		Dass[p][0][1].to(1,1) = 3;
		  Dass[p][0][1].to(2,0) = 1;		Dass[p][0][1].to(2,1) = 3;
		  Dass[p][0][1].to(3,0) = 1;		Dass[p][0][1].to(3,1) = 4;
		  Dass[p][0][1].to(4,0) = 2;		Dass[p][0][1].to(4,1) = 3;
		  Dass[p][0][1].to(5,0) = 2;		Dass[p][0][1].to(5,1) = 4;
		  Dass[p][0][1].to(6,0) = 3;		Dass[p][0][1].to(6,1) = 3;
		  Dass[p][0][1].to(7,0) = 4;		Dass[p][0][1].to(7,1) = 4;
		}
//		else if ( p == 4 )
//		{
//		  Dass[p][0][1].to(0,0) = 0;		Dass[p][0][1].to(0,1) = 3;
//		}
// these are for Raj's derivatives
		else if ( p == 4 )
		{
		  Dass[p][0][1].to(0,0) = 1;		Dass[p][0][1].to(0,1) = 3;
		  Dass[p][0][1].to(1,0) = 2;		Dass[p][0][1].to(1,1) = 3;
		  Dass[p][0][1].to(2,0) = 3;		Dass[p][0][1].to(2,1) = 3;
        }
		else if ( p == 5 )
		{
		  Dass[p][0][1].to(0,0) = 1;		Dass[p][0][1].to(0,1) = 4;
		  Dass[p][0][1].to(1,0) = 2;		Dass[p][0][1].to(1,1) = 4;
		  Dass[p][0][1].to(2,0) = 4;		Dass[p][0][1].to(2,1) = 4;
		  Dass[p][0][1].to(3,0) = 1;		Dass[p][0][1].to(3,1) = 2;
        }
	  }
	} else likelihood<realnum>::Error(106, func_name, (void*) &Num );
  }
  else if ( Dtype == "mad" ) 
  {
	if ( Num == 8 ) 
	{
	  DNum = 2; 
	  maxtoNum = 16;
	  for ( int p = 0;  p < N_part;  p++ ) 
	  {
		Dass[p] = new Dassignment*[DNum];
		for ( int i = 0;  i < DNum;  i++ ) Dass[p][i] = new Dassignment[DNum];

		// working just with the upper triangle
		if ( p == 0 ) Dass[p][0][1].toNum = 16;
		else if ( p == 1 ) Dass[p][0][1].toNum = 0;
		else likelihood<realnum>::Error( 111, func_name );

		Dass[p][0][1].to_array = new int[ 2*Dass[p][0][1].toNum ] ;
		  
		if ( p == 0 )
		{
          for (int i=0; i<Dass[p][0][1].toNum; i++) 
          {
		    Dass[p][0][1].to(i,0) = i/4;		Dass[p][0][1].to(i,1) = i%4;
          }
		}
	  }
	} else likelihood<realnum>::Error(106, func_name, (void*) &Num );
  }
  else likelihood<realnum>::Error( 108, func_name );
}


// Basic constructor

template <typename realnum>
covar_matrix<realnum>::covar_matrix( likelihood<realnum> *l, int Num_partial, int Number, int Num_meas, string Dtype_toset, 
	int do_der_D, int no_imaginary ) : llh(l) 
{
  char* func_name = (char*) "covar_matrix::covar_matrix";
  /*
  if ( Num_partial > 6 || Num_partial < 1 )		
	likelihood<realnum>::Error( 110, func_name, (void*) &Num_partial );
  */
  N_part = Num_partial;
  Num = Number;
  N_meas = Num_meas;
  if ( Dtype_toset == "sras" || Dtype_toset == "sras3d" || Dtype_toset == "sras4d" || Dtype_toset == "srasph" || 
	  Dtype_toset == "sadh" || Dtype_toset == "mad" || Dtype_toset == "p+l" || Dtype_toset == "sirph" )
  {
	no_imag = 1; // imaginary terms not supported for SIRAS or MAD
  }
  else no_imag = no_imaginary;
  if ( ( Dtype_toset == "sras" || Dtype_toset == "sirph" ) && N_part < 2 ) {
	N_part = 2;
	likelihood<realnum>::Error( 202, func_name, (void*) &N_part );
  }
  if ( ( Dtype_toset == "sras3d" ) && N_part < 3 ) {
	N_part = 3;
	likelihood<realnum>::Error( 202, func_name, (void*) &N_part );
  }
  if ( ( Dtype_toset == "sras4d" ) && N_part < 5 ) {
	N_part = 5;
	likelihood<realnum>::Error( 202, func_name, (void*) &N_part );
  }
  if ( Dtype_toset == "p+l" && Num == 4 && N_part < 4 ) {
	N_part = 4;
	likelihood<realnum>::Error( 202, func_name, (void*) &N_part );
  }
  zero_cov_row_num = 0;
// workaround for Irix
  Dtype.reserve(8);
  SetAllocDtype( Dtype_toset );
  Alloc( do_der_D ); 
}


// Copy constructor

template <typename realnum>
covar_matrix<realnum>::covar_matrix( covar_matrix<realnum>& cm_in, likelihood<realnum> *l )  : llh(l) 
{
  char* func_name = (char*) "covar_matrix::operator=";
  if ( cm_in.N_part > 6 || cm_in.N_part < 1 )	
	likelihood<realnum>::Error( 110, func_name, (void*) &(cm_in.N_part) );
  N_part = cm_in.N_part;
  Num = cm_in.Num;
  N_meas = cm_in.N_meas;
  no_imag = cm_in.no_imag;
  zero_cov_row_num = cm_in.zero_cov_row_num;
// workaround for Irix
  Dtype.reserve(8);
  SetAllocDtype( cm_in.Dtype );
  Alloc( cm_in.llh->do_der_D ); 
}


// Assignment opereator

template <typename realnum>
covar_matrix<realnum>& covar_matrix<realnum>::operator=( covar_matrix<realnum>& cm_in )
{
  char* func_name = (char*) "covar_matrix::operator=";
  if ( re != 0 )
  {
	DeAlloc( llh->do_der_D );
	DeAllocDtype( );
  }
  if ( cm_in.N_part > 6 || cm_in.N_part < 1 )	
	likelihood<realnum>::Error( 110, func_name, (void*) &(cm_in.N_part) );
  N_part = cm_in.N_part;
  Num = cm_in.Num;
  N_meas = cm_in.N_meas;
  no_imag = cm_in.no_imag;
  zero_cov_row_num = cm_in.zero_cov_row_num;
// workaround for Irix
  Dtype.reserve(8);
  SetAllocDtype( cm_in.Dtype );
  Alloc( cm_in.llh->do_der_D ); 
  return *this;
}


template <typename realnum>
covar_matrix<realnum>::~covar_matrix()
{
  DeAlloc( llh->do_der_D );
  DeAllocDtype( );
}

#undef alloc_test



template <typename realnum>
void covar_matrix<realnum>::Make_matrix( )
{
  char* func_name = (char*) "covar_matrix::Make_matrix";
  for (int i = 0;  i < llh->Num;  i++) {
    for (int j = i;  j < llh->Num;  j++) 
  	  if ( llh->do_der_D )
		for (int k=0; k<N_part; k++)  { 
		  part[k].deriv_D_r[i][j] = 0.;
		  if ( llh->do_der_D>1 )  
			for (int l=0; l<N_part; l++)   partpart[k][l].deriv_D_r[i][j] = 0.;
		}
  }
  realnum scale = 1.;
  if ( ( Dtype == "sad" || Dtype == "sadh" || Dtype == "sras4d" || Dtype == "srasph" || 
	     Dtype == "sirph" || Dtype == "sirh" || Dtype == "sir" || Dtype == "single" ) && N_part > 2 )   
	scale = part[2].D[0][1];
  realnum scale2 = 1.;
  if ( Dtype == "sadh" && N_part > 4 )
	scale2 = part[4].D[0][1];
  
  if ( llh->N_meas == 1 )
  for (int i = 0;  i < llh->Num;  i++) {
    for (int j = i;  j < llh->Num;  j++) {
      if ( i < llh->N_meas && j < llh->N_meas ) {	re[i][j] = real(sigma_N[0] + sig_meas[i]*sig_meas[i]) ;
    											im[i][j] = 0 ; 						
				    						  }
      else			      
	  {
//		int N_part_eff = N_part;
        int N_part_eff = 1;
//		if ( Dtype == "sir" || Dtype == "sirh" ) N_part_eff = N_part / 2;
		
        if ( i == j )            	{ re[i][j] = im[i][j] = 0 ; 
									  for  ( int p = 0;  p < N_part_eff;  p++ ) 
//										re[i][j] += real(part[p].sigma_P[1]) ;
//!!!!!!!!!!!!!!!!
										re[i][j] += /*scale**/real(part[p].sigma_P[0]) ;
										if ( llh->do_der_D && N_part>2 )	
										{		
//										  part[2].deriv_D_r[i][j] = real(part[0].sigma_P[0]);
										}
//  the following is an attempt to  improve the 2 part. models matrix
									  if ( N_part_eff == 2 ) re[i][j] += 2*real(partpart[0][1].sigma_P[0]);
									}	
        else  			    		{ re[i][j] = im[i][j] = 0;
//!!!!!!!!!!!
									  re[i][j] = part[0].D[i][j]* /*scale**/ real(part[0].sigma_P[0]) ;
									  if ( llh->do_der_D )	
									  {	
										part[0].deriv_D_r[i][j] = /*scale**/real(part[0].sigma_P[0]);
										part[1].deriv_D_r[i][j] = 0;
										if (N_part>2) {
//										  part[2].deriv_D_r[i][j] = part[0].D[i][j]*real(part[0].sigma_P[0]);
										  if (llh->do_der_D>1) partpart[0][2].deriv_D_r[i][j] = real(part[0].sigma_P[0]);
										}
									  }
//!!!!!!!!!!!!!!!!
//									  for  ( int p = 0;  p < N_part_eff;  p++ ) 
//									  {
//										re[i][j] += real(part[p].D[i][j]*part[p].sigma_P[0]) ;//* sqrt(real(sigma_N[0]/part[p].sigma_P[0]));
//										if ( llh->do_der_D )	
//										{	
//										  part[p].deriv_D_r[i][j] = real( part[p].sigma_P[0] ) ;//* sqrt(real(sigma_N[0]/part[p].sigma_P[0]));
//										  part[p].deriv_D_i[i][j] = 0;
//  the following is an attempt to  improve the 2 part. models matrix
//										  if ( N_part_eff == 2 ) 
//										  {
//											part[p].deriv_D_r[i][j] += sqrt( part[1-p].D[i][j]/part[p].D[i][j] )*real(partpart[0][1].sigma_P[0]);
//											partpart[p][p].deriv_D_r[i][j] = -.5*sqrt( part[1-p].D[i][j]/part[p].D[i][j] )/part[p].D[i][j]*real(partpart[0][1].sigma_P[0]);
//											partpart[p][1-p].deriv_D_r[i][j] = .5*sqrt( 1/part[1-p].D[i][j]/part[p].D[i][j] )*real(partpart[0][1].sigma_P[0]);
//										  }
//										}
//									  }  
//  the following is an attempt to  improve the 2 part. models matrix
//									  if ( N_part_eff == 2 ) re[i][j] += 2*sqrt(part[0].D[i][j]*part[1].D[i][j])*real(partpart[0][1].sigma_P[0]);
									}  
	  }
    }
  }


  if ( llh->N_meas == 2 )
  {
  // these definitions are just to make code bit more readable and slightly quicker for the most used case of 2 obs., 2 mod. and no imag. terms
  realnum partsigP_r[2][3], sigPex_r[2], sigN_r[2];
  if ( llh->Num >= 4 )
  {
    partsigP_r[0][0] = real(part[0].sigma_P[0]); 
	partsigP_r[0][1] = real(part[0].sigma_P[1]);
	partsigP_r[0][2] = real(part[0].sigma_P[2]);
	sigPex_r[0] = real(partpart[0][1].sigma_P[0]);
	if ( N_part == 2 )
	{
	  partsigP_r[1][0] = real(part[1].sigma_P[0]);
	  partsigP_r[1][1] = real(part[1].sigma_P[1]);
	  sigPex_r[1] = real(partpart[0][1].sigma_P[1]);
	}
	sigN_r[0] = real(sigma_N[0]);
	sigN_r[1] = real(sigma_N[1]);
  }
  for (int i = 0;  i < llh->Num;  i++) {
    int i_D = i/2 ;
    for (int j = i;  j < llh->Num;  j++) {
      int j_D = j/2 ;
      int primed = (i+j) % 2 ;
// making the "measured" part of cov matrix
      if ( i < llh->N_meas && j < llh->N_meas ) 
	  {	
		if ( Dtype == "mldr"  && llh->Num == 3 ) 	// inverted order of observed data supposed (derivative first)
		{
      	  if ( i == j ) 	{ int i_inv = 1-i;
							  re[i][j] = real(sigma_N[i_inv] + sig_meas[i]*sig_meas[i]) ;
    						  im[i][j] = 0 ; 
							}	
		  else				{ re[i][j] = part[0].D[i][j] * real(sigma_N[0]) ;
				    		  im[i][j] = 0 ;
							  if ( llh->do_der_D )
							  {	
								part[0].deriv_D_r[i][j] = real( sigma_N[0] );
								part[0].deriv_D_i[i][j] = 0;
							  }
							}
		}
		else if ( Dtype == "sir" || Dtype == "sirh" || Dtype == "sirph" || Dtype == "p+l" ) 
		{
      	  if ( i == j ) 	{ re[i][j] = real(sigma_N[i] + sig_meas[i]*sig_meas[i]) ; 	
    						  im[i][j] = 0 ; 
							}	
		  else				{ re[i][j] = part[1].D[i][j] * real(sigma_N[0]) ;
				    		  im[i][j] = 0 ;
							  if ( llh->do_der_D )
							  {	
								part[1].deriv_D_r[i][j] = real( sigma_N[0] );
								part[1].deriv_D_i[i][j] = 0;
						 	  }
							}
		}
		else //sad etc
		{
      	  if ( i == j ) 	{ re[i][j] = sigN_r[0] + sig_meas[i]*sig_meas[i] ; 	
    						  im[i][j] = 0 ; 
							}	
//!!!!!!!!
//		  else				{ re[i][j] = real(sigma_N[1]) ;
//				    		  if (!no_imag) im[i][j] = imag(sigma_N[1]) ;
//								else im[i][j] = 0. ;
		  else				{ re[i][j] = sigN_r[0] - (sigN_r[0]-sigN_r[1])*scale;
							  if (N_part>3&&N_part<6) re[i][j] = 
								sigN_r[0] - (sigN_r[0]-sigN_r[1])* scale* part[3].D[0][1] ;
				    		  if (!no_imag) im[i][j] = imag(sigma_N[1]) ;
								else im[i][j] = 0. ;
							  if ( llh->do_der_D  ) {	
							    if (N_part>2)
								  part[2].deriv_D_r[i][j] = - (sigN_r[0]-sigN_r[1]) ;
							    if (N_part>3&&N_part<6) {
							      part[2].deriv_D_r[i][j] *= part[3].D[0][1] ;
								  part[3].deriv_D_r[i][j] = - (sigN_r[0]-sigN_r[1])*scale ;
								  if (llh->do_der_D>1) partpart[2][3].deriv_D_r[i][j] = - (sigN_r[0]-sigN_r[1]);
								}
// for Raj
								if (N_part>5)
								  if (Num==4)      part[5].deriv_D_r[i][j] = - 2*scale*part[3].D[0][1] ;
								  else if (Num==5 && N_part>7) part[7].deriv_D_r[i][j] = - 2*scale*part[3].D[0][1] ;
                              }
							}
		}
	  }

// making the "model" part of cov matrix
      else		// if not in N_meas x N_meas upper left square
	  {
	    if ( Dtype == "single" )  	{ re[i][j] = im[i][j] = 0;
									  for  ( int p = 0;  p < N_part;  p++ ) 
									  {
										re[i][j] += real(part[p].D[i][j]*part[p].sigma_P[primed]) ; 
				    					if (!no_imag) im[i][j] += imag(part[p].D[i][j]*part[p].sigma_P[primed]) ;
										if ( llh->do_der_D )
										{
										  part[p].deriv_D_r[i][j] = real( part[p].sigma_P[primed] );
										  if (!no_imag) part[p].deriv_D_i[i][j] = imag( part[p].sigma_P[primed] );
										}
									  }
									}
		else if ( Dtype == "p+l" )
		{
		  if ( llh->Num == 4 ) 
		  {
			im[i][j] = 0 ;
		    if ( i == 3 && j==3 )
		    {  
			  re[i][j] = real(part[0].sigma_P[1]) ;
			}
			else if ( i==2 && j==2 )
			{
			  re[i][j] = real(part[0].sigma_P[0]) ;
			}
			else if ( i==2 && j==3 )
			{ // D4
//!!!!
              re[i][j] = 0.;
//			  re[i][j] = real( part[2].D[0][1] * part[0].sigma_P[1] );
//			  if ( llh->do_der_D ) 
//			    part[2].deriv_D_r[i][j] = real(part[0].sigma_P[1]) ;
			}
			else if ( i==1 && j==2) 
			{
//!!!!!!!!!!
//			  re[i][j] = real( part[0].D[0][1] * part[1].D[0][1] * part[0].sigma_P[0] ) ;
//			  if (llh->do_der_D ) {
//			    part[0].deriv_D_r[i][j] = real( part[1].D[0][1] * part[0].sigma_P[0] ) ;
//				part[1].deriv_D_r[i][j] = real( part[0].D[0][1] * part[0].sigma_P[0] ) ;
//			  }
			  re[i][j] = real( /*part[0].D[0][1] **/ part[2].D[0][1] * part[0].sigma_P[0] ) ;
			  if (llh->do_der_D ) {
			//    part[0].deriv_D_r[i][j] = real( part[2].D[0][1] * part[0].sigma_P[0] ) ;
				part[2].deriv_D_r[i][j] = real( /*part[0].D[0][1] **/ part[0].sigma_P[0] ) ;
			  }
			}
			else if ( i==1 && j==3 )
			{
			  re[i][j] = real(/* part[0].D[0][1] * */part[3].D[0][1] * part[0].sigma_P[1] ) ;
			  if ( llh->do_der_D ) {
			    //part[0].deriv_D_r[i][j] = real(part[3].D[0][1] * part[0].sigma_P[1]) ;
			    part[3].deriv_D_r[i][j] = real(/*part[0].D[0][1] **/ part[0].sigma_P[1]) ;
			  }
			}
			else if ( i==0 && j==2 )
			{
			  re[i][j] = real( part[0].D[0][1] * part[0].sigma_P[0] ) ;
			  if ( llh->do_der_D ) {
			    part[0].deriv_D_r[i][j] = real( part[0].sigma_P[0] ) ;
		  	  }
			}
			else if ( i==0 && j==3 )
			{
//!!!!!!!!
              re[i][j] = 0.;
//			  re[i][j] = real( part[0].D[0][1] * part[2].D[0][1] * part[0].sigma_P[1] ) ;
//			  if ( llh->do_der_D ) {
//			    part[0].deriv_D_r[i][j] = real( part[2].D[0][1] * part[0].sigma_P[1] ) ;
//			    part[2].deriv_D_r[i][j] = real( part[0].D[0][1] * part[0].sigma_P[1] ) ;
//			  }
			}
		  }
		  else if ( llh->Num == 3 ) 
		  {
		    if ( i==2 && j==2 )
		    {
			  re[i][j] = real(part[0].sigma_P[0]) ;
		    }
		    else if ( i==0 && j==2 )
		    {
			  re[i][j] = real(/*part[0].D[0][1] **/ part[2].D[0][1] * part[0].sigma_P[0]) ;
			  if ( llh->do_der_D ) { 
			//	part[0].deriv_D_r[i][j] = real( part[2].D[0][1] * part[0].sigma_P[0] ) ;
				part[2].deriv_D_r[i][j] = real(/* part[0].D[0][1] **/ part[0].sigma_P[0] ) ;
              }
		    }
		    else if ( i==1 && j==2 )
		    {
//!!!!
			  re[i][j] = real( /*part[1].D[0][1] **/ part[0].D[0][1] * part[0].sigma_P[0] ) ;
			  if ( llh->do_der_D ) {
//				part[1].deriv_D_r[i][j] = real( part[0].D[0][1] * part[0].sigma_P[0] ) ;
				part[0].deriv_D_r[i][j] = real( /*part[1].D[0][1] **/ part[0].sigma_P[0] ) ;
			  }
		    }
		  }
		  else likelihood<realnum>::Error(106,func_name,(void*) &llh->Num);	
		}
		else if ( Dtype == "sadh" )
		{
    	  if ( llh->Num == 4 )
		  {
			im[i][j] = 0.; 
      		if ( (i==0||i==1) && j==3 ) { re[i][j] = part[0].D[i_D][j_D] * scale * partsigP_r[0][0] ;//* sqrt(real(sigma_N[0])/partsigP_r[0][0]);
										  if ( llh->do_der_D ) /**/{
											part[0].deriv_D_r[i][j] = scale * partsigP_r[0][0] ;//* sqrt(real(sigma_N[0])/partsigP_r[0][0]);
										/**/part[1].deriv_D_r[i][j] = 0;//part[0].D[i_D][j_D] * partsigP_r[0][0] ;//* sqrt(real(sigma_N[0])/partsigP_r[0][0]);
											if ( N_part>2 ) {
											  part[2].deriv_D_r[i][j] = part[0].D[i_D][j_D] * partsigP_r[0][0];
											  if (llh->do_der_D>1) partpart[0][2].deriv_D_r[i][j] = partsigP_r[0][0];
											}
// for Raj
							  				if (N_part>5) 
							  				  part[4].deriv_D_r[i][j] = part[0].D[i_D][j_D] * scale ;
										  /**/}
										}  
      		else if ( (i==0||i==1) && j==2 ) 	
      									{ realnum sign = 1.;  if (i==1) sign = -1.;
//										  re[i][j] = part[1/*0*/].D[i_D][j_D] * partsigP_r[0][1] ;
										  re[i][j] = sign *part[0].D[i_D][j_D] * part[1].D[i_D][j_D] * scale * partsigP_r[0][1] ;//* sqrt(real(sigma_N[0])/partsigP_r[0][1]);
										  if ( llh->do_der_D ) /**/{
//											part[1/*0*/].deriv_D_r[i][j] = partsigP_r[0][1];
										/**/part[0].deriv_D_r[i][j] = sign*part[1].D[i_D][j_D] * scale * partsigP_r[0][1] ;//* sqrt(real(sigma_N[0])/partsigP_r[0][1]);
										/**/part[1].deriv_D_r[i][j] = sign*part[0].D[i_D][j_D] * scale * partsigP_r[0][1] ;//* sqrt(real(sigma_N[0])/partsigP_r[0][1]);
											if ( N_part>2 ) part[2].deriv_D_r[i][j] = sign*part[0].D[i_D][j_D] * part[1].D[i_D][j_D] * partsigP_r[0][1];
											//for mixed 2. derivs - not used as of now
											if ( llh->do_der_D>1 ) {
											  partpart[0][1].deriv_D_r[i][j] = sign* scale * partsigP_r[0][1] ;
											  if ( N_part>2 ) partpart[0][2].deriv_D_r[i][j] = sign*part[1].D[i_D][j_D] * partsigP_r[0][1] ;
											  if ( N_part>2 ) partpart[1][2].deriv_D_r[i][j] = sign*part[0].D[i_D][j_D] * partsigP_r[0][1] ;
											}
// for Raj
							  				if (N_part>5) 
							  				  part[5].deriv_D_r[i][j] = sign*part[0].D[i_D][j_D]*part[1].D[i_D][j_D] * scale ;
										  /**/}
										}  
      		else if ( i==3 && j==3 )    { re[i][j] = scale*partsigP_r[0][0] ;
										  if ( llh->do_der_D && N_part>2 ) 
										    part[2].deriv_D_r[i][j] = partsigP_r[0][0];
// for Raj
							  			  if (llh->do_der_D && N_part>5) 
							  				part[4].deriv_D_r[i][j] = scale ;
										}	
			else if ( i==2 && j==2 )	{ re[i][j] = scale*partsigP_r[0][1] ;
										  if ( llh->do_der_D && N_part>2 )
										    part[2].deriv_D_r[i][j] = partsigP_r[0][1];
// for Raj
							  			  if (llh->do_der_D && N_part>5) 
							  				part[5].deriv_D_r[i][j] = scale ;
										}
      		else if ( i==2 && j==3 ) 	  re[i][j] = 0.;  
		  }
    	  else if ( llh->Num == 5 )
		  {
			im[i][j] = 0.; 
      		if ( (i==0||i==1) && j==3 ) { re[i][j] = part[0].D[i_D][j_D] * scale2 * partsigP_r[0][0] ;//* sqrt(real(sigma_N[0])/partsigP_r[0][0]);
										  if ( llh->do_der_D ) /**/{
											part[0].deriv_D_r[i][j] = scale2 * partsigP_r[0][0] ;//* sqrt(real(sigma_N[0])/partsigP_r[0][0]);
										/**/part[1].deriv_D_r[i][j] = 0;//part[0].D[i_D][j_D] * partsigP_r[0][0] ;//* sqrt(real(sigma_N[0])/partsigP_r[0][0]);
											if ( N_part>2 ) {
											  part[4].deriv_D_r[i][j] = part[0].D[i_D][j_D] * partsigP_r[0][0];
											  if (llh->do_der_D>1) partpart[0][4].deriv_D_r[i][j] = partsigP_r[0][0];
											}
										  /**/}
										}  
      		else if ( (i==0||i==1) && j==2 ) 	
      									{ realnum sign = 1.;  if (i==1) sign = -1.;
//										  re[i][j] = part[1/*0*/].D[i_D][j_D] * partsigP_r[0][1] ;
										  re[i][j] = sign /**part[0].D[i_D][j_D]*/ * part[1].D[i_D][j_D] * scale * partsigP_r[0][1] ;//* sqrt(real(sigma_N[0])/partsigP_r[0][1]);
										  if ( llh->do_der_D ) /**/{
//											part[1/*0*/].deriv_D_r[i][j] = partsigP_r[0][1];
										//part[0].deriv_D_r[i][j] = sign*part[1].D[i_D][j_D] * scale * partsigP_r[0][1] ;//* sqrt(real(sigma_N[0])/partsigP_r[0][1]);
										/**/part[1].deriv_D_r[i][j] = sign/**part[0].D[i_D][j_D]*/ * scale * partsigP_r[0][1] ;//* sqrt(real(sigma_N[0])/partsigP_r[0][1]);
											if ( N_part>2 ) part[2].deriv_D_r[i][j] = sign/**part[0].D[i_D][j_D]*/ * part[1].D[i_D][j_D] * partsigP_r[0][1];
											//for mixed 2. derivs - not used as of now
											if ( llh->do_der_D>1 ) {
											  //partpart[0][1].deriv_D_r[i][j] = sign* scale * partsigP_r[0][1] ;
											  //if ( N_part>2 ) partpart[0][2].deriv_D_r[i][j] = sign*part[1].D[i_D][j_D] * partsigP_r[0][1] ;
											  if ( N_part>2 ) partpart[1][2].deriv_D_r[i][j] = sign/**part[0].D[i_D][j_D]*/ * partsigP_r[0][1] ;
											}
										  /**/}
// for Raj
							  				if (N_part>7) 
							  				  part[7].deriv_D_r[i][j] = sign*part[0].D[i_D][j_D]*part[1].D[i_D][j_D] * scale ;
										}  
      		else if ( i==3 && j==3 )    { re[i][j] = scale2*partsigP_r[0][0];
										  if ( llh->do_der_D && N_part>2 ) {
										    part[4].deriv_D_r[i][j] = partsigP_r[0][0] ;
										  }
										}	
			else if ( i==2 && j==2 )	{ re[i][j] = scale*partsigP_r[0][1] ;
										  if ( llh->do_der_D && N_part>2 ) {
										    part[2].deriv_D_r[i][j] = partsigP_r[0][1];
										  }
// for Raj
							  			  if (llh->do_der_D && N_part>7) 
							  				part[7].deriv_D_r[i][j] = scale ;
										}
      		else if ( i==2 && (j==3||j==4) ) 	  
      									  re[i][j] = 0.;  
      		else if ( (i==0||i==1) && j==4 ) 	
      									{ re[i][j] = /*part[0].D[0][1] **/ part[5].D[0][1] * scale * partsigP_r[0][2] ;//* sqrt(real(sigma_N[0])/partsigP_r[0][0]);
										  if ( llh->do_der_D ) /**/{
											//part[0].deriv_D_r[i][j] = part[5].D[0][1] * scale * partsigP_r[0][2] ;//* sqrt(real(sigma_N[0])/partsigP_r[0][0]);
											if (N_part>5) part[5].deriv_D_r[i][j] = /*part[0].D[0][1] **/ scale * partsigP_r[0][2] ;//* sqrt(real(sigma_N[0])/partsigP_r[0][0]);
											if (N_part>4) part[2].deriv_D_r[i][j] = /*part[0].D[0][1] **/ part[5].D[0][1] * partsigP_r[0][2] ;//* sqrt(real(sigma_N[0])/partsigP_r[0][0]);
											//if (llh->do_der_D>1&&N_part>4) partpart[0][2].deriv_D_r[i][j] = part[5].D[0][1] * partsigP_r[0][2];
											//if (llh->do_der_D>1&&N_part>5) partpart[0][5].deriv_D_r[i][j] = part[4].D[0][1] * partsigP_r[0][2];
											if (llh->do_der_D>1&&N_part>5) partpart[2][5].deriv_D_r[i][j] = /*part[0].D[0][1] **/ partsigP_r[0][2];
// for Raj
							  				if (N_part>7) 
							  				  part[6].deriv_D_r[i][j] = part[5].D[0][1] * scale ;
///      									{ re[i][j] = 0.1*part[0].D[0][1] * part[1].D[0][1] * scale * partsigP_r[0][2] ;//* sqrt(real(sigma_N[0])/partsigP_r[0][0]);
///										  if ( llh->do_der_D ) /**/{
///											part[0].deriv_D_r[i][j] = 0.1*part[1].D[0][1] * scale * partsigP_r[0][2] ;//* sqrt(real(sigma_N[0])/partsigP_r[0][0]);
///											if (N_part>5) part[1].deriv_D_r[i][j] = 0.1*part[0].D[0][1] * scale * partsigP_r[0][2] ;//* sqrt(real(sigma_N[0])/partsigP_r[0][0]);
///											if (N_part>4) part[2].deriv_D_r[i][j] = 0.1*part[0].D[0][1] * part[1].D[0][1] * partsigP_r[0][2] ;//* sqrt(real(sigma_N[0])/partsigP_r[0][0]);
///											if (llh->do_der_D>1&&N_part>4) partpart[0][2].deriv_D_r[i][j] = 0.1*part[1].D[0][1] * partsigP_r[0][2];
///											if (llh->do_der_D>1&&N_part>5) partpart[0][1].deriv_D_r[i][j] = 0.1*part[4].D[0][1] * partsigP_r[0][2];
///											if (llh->do_der_D>1&&N_part>5) partpart[1][2].deriv_D_r[i][j] = 0.1*part[0].D[0][1] * partsigP_r[0][2];
										  /**/}
										}  
      		else if ( i==3 && j==4 ) 	{ re[i][j] = part[3].D[0][1] * scale * partsigP_r[0][2] ;//* sqrt(real(sigma_N[0])/partsigP_r[0][0]);
										  if ( llh->do_der_D ) /**/{
											if (N_part>3) part[3].deriv_D_r[i][j] = scale * partsigP_r[0][2] ;//* sqrt(real(sigma_N[0])/partsigP_r[0][0]);
											if (N_part>4) part[2].deriv_D_r[i][j] = part[3].D[0][1] * partsigP_r[0][2] ;//* sqrt(real(sigma_N[0])/partsigP_r[0][0]);
											if (llh->do_der_D>1&&N_part>4) partpart[2][3].deriv_D_r[i][j] = partsigP_r[0][2];
// for Raj
							  				if (N_part>7) 
							  				  part[6].deriv_D_r[i][j] = part[3].D[0][1] * scale ;
										  /**/}
										}  
      		else if ( i==4 && j==4 ) 	{ re[i][j] = scale * partsigP_r[0][2] ;//* sqrt(real(sigma_N[0])/partsigP_r[0][0]);
										  if ( llh->do_der_D ) /**/{
											if (N_part>4) part[2].deriv_D_r[i][j] = partsigP_r[0][2] ;//* sqrt(real(sigma_N[0])/partsigP_r[0][0]);
// for Raj
							  				if (N_part>7) 
							  				  part[6].deriv_D_r[i][j] = scale ;
										  /**/}
										}  
		  }
		  else likelihood<realnum>::Error(106,func_name,(void*) &llh->Num);	
        }
		else if ( Dtype == "sad" )
		{
    	  if ( llh->Num == 3 )
		  {
      		if ( i == j )            	{ re[i][j] = im[i][j] = 0;
										  for  ( int p = 0;  p < N_part;  p++ ) 
											re[i][j] += real(part[p].sigma_P[0]) ;
//  the following is an attempt to  improve the 2 part. models matrix
										  if ( N_part == 2 )  re[i][j] += 2*real(partpart[0][1].sigma_P[0]);
										}	
			else if ( i == 0 )			{ re[i][j] = im[i][j] = 0;
										  for  ( int p = 0;  p < N_part;  p++ )
										  {
											re[i][j] += part[p].D[i_D][j_D] * real( part[p].sigma_P[1] ) ;
											if (!no_imag) im[i][j] += - part[p].D[i_D][j_D] * imag( part[p].sigma_P[1] ) ; 
											if ( llh->do_der_D )
											{	
											  part[p].deriv_D_r[i][j] = real( part[p].sigma_P[1] );
											  if (!no_imag) part[p].deriv_D_i[i][j] = - imag( part[p].sigma_P[1] );
//  the following is an attempt to  improve the 2 part. models matrix
											  if ( N_part == 2 )
											  {
												part[p].deriv_D_r[i][j] += sqrt( part[1-p].D[i_D][j_D]/part[p].D[i_D][j_D] )*real(partpart[0][1].sigma_P[1]);
												partpart[p][p].deriv_D_r[i][j] = -.5*sqrt( part[1-p].D[i_D][j_D]/part[p].D[i_D][j_D] )/part[p].D[i_D][j_D]*real(partpart[0][1].sigma_P[1]);
												partpart[p][1-p].deriv_D_r[i][j] = .5*sqrt( 1/part[1-p].D[i_D][j_D]/part[p].D[i_D][j_D] )*real(partpart[0][1].sigma_P[1]);
											  }
											}
										  }
//  the following is an attempt to  improve the 2 part. models matrix
										  if ( N_part == 2 ) re[i][j] += 2*sqrt(part[0].D[i_D][j_D]*part[1].D[i_D][j_D])*real(partpart[0][1].sigma_P[1]);
										}
			else if ( i == 1 )			{ re[i][j] = im[i][j] = 0;
										  for  ( int p = 0;  p < N_part;  p++ )
										  {
											re[i][j] += part[p].D[i_D][j_D] * real( part[p].sigma_P[1] ) ;
											if (!no_imag) im[i][j] += part[p].D[i_D][j_D] * imag( part[p].sigma_P[1] ) ;
											if ( llh->do_der_D )
											{	
											  part[p].deriv_D_r[i][j] = real( part[p].sigma_P[1] );
											  if (!no_imag) part[p].deriv_D_i[i][j] = imag( part[p].sigma_P[1] );
//  the following is an attempt to  improve the 2 part. models matrix
											  if ( N_part == 2 )
											  {
												part[p].deriv_D_r[i][j] += sqrt( part[1-p].D[i_D][j_D]/part[p].D[i_D][j_D] )*real(partpart[0][1].sigma_P[1]);
												partpart[p][p].deriv_D_r[i][j] = -.5*sqrt( part[1-p].D[i_D][j_D]/part[p].D[i_D][j_D] )/part[p].D[i_D][j_D]*real(partpart[0][1].sigma_P[1]);
												partpart[p][1-p].deriv_D_r[i][j] = .5*sqrt( 1/part[1-p].D[i_D][j_D]/part[p].D[i_D][j_D] )*real(partpart[0][1].sigma_P[1]);
											  }
											}
										  }
//  the following is an attempt to  improve the 2 part. models matrix
										  if ( N_part == 2 ) re[i][j] += 2*sqrt(part[0].D[i_D][j_D]*part[1].D[i_D][j_D])*real(partpart[0][1].sigma_P[1]);
										}
		  }
    	  else if ( llh->Num == 4 )
		  {
      		if ( i == j )            	{ re[i][j] = im[i][j] = 0;    
										  for  ( int p = 0;  p < N_part;  p++ )
											re[i][j] += partsigP_r[p][0] ;
//  the following is an attempt to  improve the 2 part. models matrix
//										  if ( N_part == 2 ) re[i][j] += 2 * sigPex_r[0] ;
										}	
			else if ( j-i==1 && i%2==0 ){ re[i][j] = im[i][j] = 0;
										  for  ( int p = 0;  p < N_part;  p++ )
										  {
											re[i][j] += partsigP_r[p][1] ;
    										if (!no_imag) im[i][j] += imag( part[p].sigma_P[1] ) ; 
										  }
//  the following is an attempt to  improve the 2 part. models matrix
//										  if ( N_part == 2 ) re[i][j] += 2 * sigPex_r[0] ;
										}
      		else  			    		{ re[i][j] = im[i][j] = 0;   
										  realnum im_sign = 1;	if ( i%2 == 1 ) im_sign = -1;
										  for  ( int p = 0;  p < N_part;  p++ )
										  {
											re[i][j] += part[p].D[i_D][j_D] * partsigP_r[p][primed] ;
				    						if (!no_imag) im[i][j] += im_sign*imag( part[p].D[i_D][j_D] * part[p].sigma_P[primed] ) ;
											if ( llh->do_der_D )
											{	
											  part[p].deriv_D_r[i][j] = partsigP_r[p][primed];
											  if (!no_imag) part[p].deriv_D_i[i][j] = im_sign*imag( part[p].sigma_P[primed] );
//  the following is an attempt to  improve the 2 part. models matrix
//											  if ( N_part == 2 )
//											  {
//												part[p].deriv_D_r[i][j] += sqrt( part[1-p].D[i_D][j_D]/part[p].D[i_D][j_D] )*sigPex_r[1];
//												partpart[p][p].deriv_D_r[i][j] = -.5*sqrt( part[1-p].D[i_D][j_D]/part[p].D[i_D][j_D] )/part[p].D[i_D][j_D]*sigPex_r[1];
//												partpart[p][1-p].deriv_D_r[i][j] = .5*sqrt( 1/part[1-p].D[i_D][j_D]/part[p].D[i_D][j_D] )*sigPex_r[1];
//											  }
											}
										  }
//  the following is an attempt to  improve the 2 part. models matrix
//										  if ( N_part == 2 ) re[i][j] += 2*sqrt(part[0].D[i_D][j_D]*part[1].D[i_D][j_D])*sigPex_r[1];
										}  
		  }
		  else likelihood<realnum>::Error(106,func_name,(void*) &llh->Num);	
		}
		else if ( Dtype == "sir" )		// 2 part. models not set up so far
		{
		  if ( llh->Num == 4 )
		  {
      		if ( i == j )            	{ re[i][j] = im[i][j] = 0;
										  re[i][j] += partsigP_r[0][i-llh->N_meas] *scale;
										  if ( llh->do_der_D ) {	
											part[2].deriv_D_r[i][j] = partsigP_r[0][i-llh->N_meas];
											part[2].deriv_D_i[i][j] = 0;
										  }
										}	
			else if ( j-i==1 && i%2==0 ){ re[i][j] = im[i][j] = 0;
										  re[i][j] += part[1].D[0][1]*partsigP_r[0][0] *scale;
										  if ( llh->do_der_D ) {	
											part[0].deriv_D_r[i][j] = part[0].deriv_D_i[i][j] = 0;
											part[1].deriv_D_r[i][j] = partsigP_r[0][0] *scale;
											part[1].deriv_D_i[i][j] = 0;
											part[2].deriv_D_r[i][j] = part[1].D[0][1]*partsigP_r[0][0] ;
											part[2].deriv_D_i[i][j] = 0;
										  }
										}
      		else if ( (i+j)%2==0 ) 		{ re[i][j] = im[i][j] = 0;
										  re[i][j] +=  part[0].D[0][1] * partsigP_r[0][i] *scale;
										  if ( llh->do_der_D )
										  {	
											part[0].deriv_D_r[i][j] = partsigP_r[0][i] *scale;
											part[0].deriv_D_i[i][j] = 0;
											part[1].deriv_D_r[i][j] = part[1].deriv_D_i[i][j] = 0;
											part[2].deriv_D_r[i][j] = part[0].D[0][1] * partsigP_r[0][i] ;
											part[2].deriv_D_i[i][j] = 0;
										  }
										}  
      		else 				 		{ re[i][j] = im[i][j] = 0;
										  re[i][j] +=  part[0].D[0][1] * part[1].D[0][1] * partsigP_r[0][0] *scale;
										  if ( llh->do_der_D )
										  {	
											part[0].deriv_D_r[i][j] = part[1].D[0][1] * partsigP_r[0][0] *scale;
											part[0].deriv_D_i[i][j] = 0;
											part[1].deriv_D_r[i][j] = part[0].D[0][1] * partsigP_r[0][0] *scale;
											part[1].deriv_D_i[i][j] = 0;
											partpart[0][1].deriv_D_r[i][j] = partsigP_r[0][0] *scale;
											partpart[1][0].deriv_D_r[i][j] = partsigP_r[0][0] *scale;
											part[2].deriv_D_r[i][j] = part[0].D[0][1] * part[1].D[0][1] * partsigP_r[0][0] ;
											part[2].deriv_D_i[i][j] = 0;
										  }
										}  
		  }
		  else likelihood<realnum>::Error(106,func_name,(void*) &llh->Num);			  
		}
		else if ( Dtype == "sirph" )
		{
		  if ( llh->Num == 3 )
		  {
      		if ( i == j )           { re[i][j] = im[i][j] = 0;
									  re[i][j] += scale* real(part[0].sigma_P[1]) ;
									  if ( llh->do_der_D ) {	
									    part[2].deriv_D_r[i][j] = real(part[0].sigma_P[1]);
									  }
									}	
			else if ( i==0 )		{ re[i][j] = im[i][j] = 0;
									}
      		else if ( i==1 ) 		{ re[i][j] = im[i][j] = 0;
									  re[i][j] +=  scale* part[0].D[0][1] * real(part[0].sigma_P[1]) ;
									  if ( llh->do_der_D ) {	
										part[0].deriv_D_r[i][j] = scale * real(part[0].sigma_P[1]) ;
										part[2].deriv_D_r[i][j] = part[0].D[0][1]* real(part[0].sigma_P[1]) ;
									  }
									}  
		  }
		  else likelihood<realnum>::Error(106,func_name,(void*) &llh->Num);			  
		}
		else if ( Dtype == "mldr" )		// 2 part. models not set up so far
		{
		  if ( llh->Num == 3 )
		  {
      		if ( i == 2 )            	{ re[i][j] = im[i][j] = 0;
										  re[i][j] += real(part[0].sigma_P[0]) ;
										}	
      		if ( i == 1 )            	{ re[i][j] = im[i][j] = 0;
										  re[i][j] += real(part[0].sigma_P[1]) ;
										}	
      		if ( i == 0 )            	{ re[i][j] = im[i][j] = 0;
										  re[i][j] += real(part[0].sigma_P[0]) ;
										}	
		  }
		}
		else if ( Dtype == "sirh" )		// 2 part. models not set up so far
		{
		  if ( llh->Num == 4 )
		  {
      		if ( i == j )            	{ re[i][j] = im[i][j] = 0;
										  re[i][j] += partsigP_r[0][i-llh->N_meas]*scale ;
										  if ( llh->do_der_D ) {
											part[2].deriv_D_r[i][j] = partsigP_r[0][i-llh->N_meas];
										  }
										}	
			else if ( i==2 && j==3 )	{ re[i][j] = im[i][j] = 0;
										}
      		else if ( i==0 && j==2 ) 	{ re[i][j] = im[i][j] = 0;
										  re[i][j] +=  part[0].D[0][1] * partsigP_r[0][0] * scale;
										  if ( llh->do_der_D )
										  {	
											part[0].deriv_D_r[i][j] = partsigP_r[0][0] *scale;
											part[2].deriv_D_r[i][j] = part[0].D[0][1] * partsigP_r[0][0];
										  }
										}  
			else if ( i==0 && j==3 )	{ re[i][j] = im[i][j] = 0;
										}
      		else if ( i==1 && j==2 ) 	{ re[i][j] = im[i][j] = 0;
//										  re[i][j] +=  part[0].D[0][1] * part[1].D[0][1] * partsigP_r[0][0] * scale;
										  re[i][j] +=  part[4].D[0][1] * partsigP_r[0][0] * scale;
										  if ( llh->do_der_D )
										  {	
//											part[0].deriv_D_r[i][j] = part[1].D[0][1] * partsigP_r[0][0] *scale;
//											part[1].deriv_D_r[i][j] = part[0].D[0][1] * partsigP_r[0][0] *scale;
//											partpart[0][1].deriv_D_r[i][j] = partsigP_r[0][0] *scale;
//											partpart[1][0].deriv_D_r[i][j] = partsigP_r[0][0] *scale;
//											part[2].deriv_D_r[i][j] = part[0].D[0][1] * part[1].D[0][1] * partsigP_r[0][0];
											part[4].deriv_D_r[i][j] = partsigP_r[0][0] * scale;
										  }
										}  
      		else if ( i==1 && j==3 ) 	{ re[i][j] = im[i][j] = 0;
										  re[i][j] +=  part[3].D[0][1]* part[0].D[0][1]*  partsigP_r[0][1] *scale;
										  if ( llh->do_der_D )
										  {	
											part[0].deriv_D_r[i][j] = scale* part[3].D[0][1]* partsigP_r[0][1];
											part[2].deriv_D_r[i][j] = part[3].D[0][1]* part[0].D[0][1]* partsigP_r[0][1];
											part[3].deriv_D_r[i][j] = scale* part[0].D[0][1]* partsigP_r[0][1];
										  }
										}  
		  }
		  else likelihood<realnum>::Error(106,func_name,(void*) &llh->Num);			  
		}
	  }
    }
  }
  }

  if ( llh->N_meas == 3 ) 
  {
	if ( Dtype == "sras" || Dtype == "sras3d" || Dtype == "sras4d" || Dtype == "srasph" )
	{
  // just to make code bit more readable and slightly quicker for the most used (and currently the only supported) 
  // case of 3 obs., 3 mod. and no imag. terms and 1 part. model (although 2 might be requested for SIRAS,
  // but just to enable 2 different D parameters in the same matrix term)
	  realnum sigP[3];
	  if ( llh->Num == 6 || (llh->Num == 5 && Dtype == "srasph") )
	  {
  		sigP[0] = real(part[0].sigma_P[0]); 
		sigP[1] = real(part[0].sigma_P[1]);
		sigP[2] = real(part[0].sigma_P[2]);
	  }
	  else likelihood<realnum>::Error(106,func_name,(void*) &llh->Num);			  
	  for (int i=0; i<llh->Num; i++) 
		for (int j=i; j<llh->Num; j++) 
		  im[i][j]=0;
	  re[0][0] = real(sigma_N[0]) + sig_meas[0]*sig_meas[0]; 
	  re[1][1] = real(sigma_N[1]) + sig_meas[1]*sig_meas[1];
	  re[2][2] = real(sigma_N[1]) + sig_meas[2]*sig_meas[2];
	  re[0][1] = part[1].D[0][1]*real(sigma_N[0]); 
	  re[0][2] = part[1].D[0][1]*real(sigma_N[0]); 

	  re[0][3] = part[0].D[0][1]*sigP[0];
	  
	  if ( Dtype == "sras" || Dtype == "sras3d" ) 
	  {
		re[1][2] =  real(sigma_N[2]);
	    if ( Dtype == "sras" )
		  re[0][4] = re[0][5] = re[1][3] = re[2][3] = part[0].D[0][1]*part[1].D[0][1]*sigP[0];
	    else if ( Dtype == "sras3d" )
		  re[0][4] = re[0][5] = re[1][3] = re[2][3] = part[2].D[0][1]*sigP[0];
		re[1][4] = re[2][5] = part[0].D[0][1]*sigP[1];
		re[1][5] = re[2][4] = part[0].D[0][1]*sigP[2];
		re[3][3] = sigP[0];
		re[3][4] = re[3][5] = part[1].D[0][1]*sigP[0];
		re[4][4] = re[5][5] = sigP[1];
		re[4][5] = sigP[2];
	  } else if ( Dtype == "sras4d" )
	  { 
		re[1][2] = real(sigma_N[1]) - real(sigma_N[1]-sigma_N[2])*scale;
      //  re[0][3] *= scale;
/*!!!!!*/	    re[0][4] = re[1][3] = re[2][3] = /*scale**/ part[4].D[0][1]*sigP[0];
//*!!!!!*/	    re[0][4] = 0.; re[1][3] = re[2][3] = /*scale**/ part[4].D[0][1]*sigP[0];
//*!!!!!*/ re[0][4] = /*scale**/ part[4].D[0][1]*sigP[0];
	    re[0][5] = 0.;
		re[1][4] = re[2][4] = /*scale**/ part[0].D[0][1]*sigP[1] /*!!!! *part[5].D[0][1]*/ ;
		re[1][5] = part[0].D[0][1]*  /*scale**/ part[3].D[0][1] *  sigP[2];
		re[2][5] = - part[0].D[0][1]*  /*scale**/ part[3].D[0][1] *  sigP[2]; 
/*!!!!!!!!!*/		re[3][4] = /*scale**/ part[1].D[0][1]* sigP[0];
//*!!!!!!*/		re[3][4] = 0.;
		re[3][3] = /*scale**/ sigP[0];
		re[4][4] = /*scale**/ sigP[1];
		re[5][5] = /*scale**/ sigP[2];
		re[4][5] = re[3][5] = 0.;
	  } else if ( Dtype == "srasph" )
	  { 
		re[1][2] = real(sigma_N[1]) - real(sigma_N[1]-sigma_N[2])*scale  /**part[3].D[0][1]*/;
		re[0][3] = scale* part[1].D[0][1]*part[0].D[0][1]*sigP[0];
////		re[0][3] = scale* part[4].D[0][1]*part[0].D[0][1]*sigP[0];
		re[0][4] = 0.;
		re[1][3] = scale* part[0].D[0][1]*sigP[1];
		re[2][3] = scale* part[0].D[0][1]*sigP[1];
/**/		re[1][4] = /*part[0].D[0][1]* */ scale* part[3].D[0][1] * sigP[2];
/**/		re[2][4] = - /*part[0].D[0][1]* */ scale* part[3].D[0][1] * sigP[2]; 
		re[3][3] = scale* sigP[1];
		re[4][4] = scale* sigP[2];
		re[3][4] = 0.;
	  }	  
	  
	  if ( llh->do_der_D )
	  {	
		for ( int i=0; i<llh->Num; i++ )  		
		  for ( int j=i; j<llh->Num; j++ )
		  {
			part[0].deriv_D_r[i][j] = part[0].deriv_D_i[i][j] = part[1].deriv_D_r[i][j] = part[1].deriv_D_i[i][j] = 0;
			if ( Dtype == "sras3d" || Dtype == "sras4d" )   
			  part[2].deriv_D_r[i][j] = part[2].deriv_D_i[i][j] = 0;
			if ( Dtype == "sras4d" )   {
			  part[3].deriv_D_r[i][j] = part[3].deriv_D_i[i][j] = 0;
			  part[4].deriv_D_r[i][j] = part[4].deriv_D_i[i][j] = 0;
			}
		  }
		part[1].deriv_D_r[0][1] = real(sigma_N[0]); 
		part[1].deriv_D_r[0][2] = real(sigma_N[0]); 
		if ( Dtype != "srasph" ) part[0].deriv_D_r[0][3] = /*scale**/sigP[0];
	    if ( Dtype == "sras" || Dtype == "sras3d" ) 
	    {
		  if ( Dtype == "sras" ) {
		    part[0].deriv_D_r[0][4] = part[0].deriv_D_r[0][5] = 
		    part[0].deriv_D_r[1][3] = part[0].deriv_D_r[2][3] = part[1].D[0][1]*sigP[0] ;
		    part[1].deriv_D_r[0][4] = part[1].deriv_D_r[0][5] = 
		    part[1].deriv_D_r[1][3] = part[1].deriv_D_r[2][3] = part[0].D[0][1]*sigP[0] ;
		  }
          else if ( Dtype == "sras3d" )
		    part[2].deriv_D_r[0][4] = part[2].deriv_D_r[0][5] = 
		    part[2].deriv_D_r[1][3] = part[2].deriv_D_r[2][3] = sigP[0] ;
		  part[0].deriv_D_r[1][4] = part[0].deriv_D_r[2][5] = sigP[1];
		  part[0].deriv_D_r[1][5] = part[0].deriv_D_r[2][4] = sigP[2];
		  part[1].deriv_D_r[3][4] = part[1].deriv_D_r[3][5] = sigP[0];
		} else if ( Dtype == "sras4d" )
		{ 
		  part[2].deriv_D_r[1][2] = - real(sigma_N[1]-sigma_N[2]);
//*!!!!!!!*/		  part[4].deriv_D_r[0][4] = part[4].deriv_D_r[1][3] = part[4].deriv_D_r[2][3] = /*scale**/sigP[0] ;
/*!!!!!!!!*/		  part[4].deriv_D_r[0][4] = 0; part[4].deriv_D_r[1][3] = part[4].deriv_D_r[2][3] = /*scale**/sigP[0] ;
//*!!!!!!!!*/ part[4].deriv_D_r[0][4] = /*scale**/sigP[0] ;
		  part[0].deriv_D_r[1][4] = part[0].deriv_D_r[2][4] = /*scale**/sigP[1] /*!!!!! *part[5].D[0][1]*/;
//*!!!*/		  part[5].deriv_D_r[1][4] = part[5].deriv_D_r[2][4] = /*scale**/sigP[1] /*!!!!!*/ *part[0].D[0][1];
		  part[3].deriv_D_r[1][5] = part[0].D[0][1]* /*scale**/sigP[2];
		  part[3].deriv_D_r[2][5] = - part[0].D[0][1]* /*scale**/sigP[2];
		  part[0].deriv_D_r[1][5] = part[3].D[0][1]* /*scale**/ sigP[2];
		  part[0].deriv_D_r[2][5] = - part[3].D[0][1]* /*scale**/ sigP[2];
/*!!!!!*/		  part[1].deriv_D_r[3][4] = /*scale**/sigP[0];
		//  part[2].deriv_D_r[0][4] = part[2].deriv_D_r[1][3] = part[2].deriv_D_r[2][3] = part[4].D[0][1]*sigP[0] ;
        //  part[2].deriv_D_r[1][4] = part[2].deriv_D_r[2][4] = part[0].D[0][1]*sigP[1];
//		  part[2].deriv_D_r[1][5] = /*part[0].D[0][1]**/  part[3].D[0][1]* sigP[2];
//		  part[2].deriv_D_r[2][5] = - /*part[0].D[0][1]**/  part[3].D[0][1]* sigP[2];
		//  part[2].deriv_D_r[1][5] = part[0].D[0][1]*  part[3].D[0][1]* sigP[2];
		//  part[2].deriv_D_r[2][5] = - part[0].D[0][1]*  part[3].D[0][1]* sigP[2];
		//  part[2].deriv_D_r[3][4] = part[1].D[0][1]* sigP[0];
		//  part[2].deriv_D_r[3][3] = sigP[0];
		//  part[2].deriv_D_r[4][4] = sigP[1];
		//  part[2].deriv_D_r[5][5] = sigP[2];
		} else if ( Dtype == "srasph" )
		{ 
		  part[2].deriv_D_r[1][2] = - real(sigma_N[1]-sigma_N[2]) /** *part[3].D[0][1]*/;
		  part[0].deriv_D_r[0][3] = scale* part[1].D[0][1]*sigP[0];
		  part[1].deriv_D_r[0][3] = scale* part[0].D[0][1]*sigP[0];
////		  part[0].deriv_D_r[0][3] = scale* part[4].D[0][1]*sigP[0];
////		  part[4].deriv_D_r[0][3] = scale* part[0].D[0][1]*sigP[0];
//		  part[3].deriv_D_r[1][2] = - real(sigma_N[1]-sigma_N[2])*scale;
		  part[0].deriv_D_r[1][3] = part[0].deriv_D_r[2][3] = scale*sigP[1];
//		  part[0].deriv_D_r[1][4] = scale* part[3].D[0][1]*sigP[2];
//		  part[0].deriv_D_r[2][4] = - scale* part[3].D[0][1]*sigP[2];
		  part[3].deriv_D_r[1][4] = /*part[0].D[0][1]* */ scale*sigP[2];
		  part[3].deriv_D_r[2][4] = - /*part[0].D[0][1]* */ scale*sigP[2];
          part[2].deriv_D_r[1][3] = part[0].D[0][1]*sigP[1];
          part[2].deriv_D_r[2][3] = part[0].D[0][1]*sigP[1];
		  part[2].deriv_D_r[1][4] = /*part[0].D[0][1]**/  part[3].D[0][1]* sigP[2];
	  	  part[2].deriv_D_r[2][4] = - /*part[0].D[0][1]**/  part[3].D[0][1]* sigP[2];
		  part[2].deriv_D_r[3][3] = sigP[1];
		  part[2].deriv_D_r[4][4] = sigP[2];
          part[2].deriv_D_r[0][3] = part[0].D[0][1]*part[1].D[0][1]*sigP[0];
////          part[2].deriv_D_r[0][3] = part[0].D[0][1]*part[4].D[0][1]*sigP[0];
// these are for Raj's derivatives
                  if ( N_part > 4) {
		    part[4].deriv_D_r[1][3] = scale*part[0].D[0][1];
		    part[4].deriv_D_r[2][3] = scale*part[0].D[0][1];
		    part[4].deriv_D_r[3][3] = scale;
		    part[5].deriv_D_r[1][4] = scale*part[3].D[0][1]/**part[0].D[0][1]*/;
		    part[5].deriv_D_r[2][4] = - scale*part[3].D[0][1]/**part[0].D[0][1]*/;
		    part[5].deriv_D_r[4][4] = scale;
		    part[5].deriv_D_r[1][2] = - 2*scale;
		  }
		}
		
	  }
// copy all deriv terms wrt D terms to lower triangle
  	  if ( llh->do_der_D ) {
		for (int i = 0;  i < llh->Num;  i++) 
  		  for (int j = i;  j < llh->Num;  j++) 
			for (int k=0; k<N_part; k++)  { 
			  part[k].deriv_D_r[j][i] = part[k].deriv_D_r[i][j];
			  if ( llh->do_der_D>1 )  
				for (int l=0; l<N_part; l++)   partpart[k][l].deriv_D_r[j][i] = partpart[k][l].deriv_D_r[i][j];
	  }
  }
	}
	else likelihood<realnum>::Error(106,func_name,(void*) &llh->Num);
  }

  if ( llh->N_meas == 4 ) 
  {
	if ( Dtype == "mad" )
	{
  // just to make code bit more readable and slightly quicker for the most used (and currently the only supported) 
  // case of 4 obs., 4 mod. and no imag. terms
	  realnum sigP[6];
	  if ( llh->Num == 8 )
	  {
  		sigP[0] = real(part[0].sigma_P[0]); 
		sigP[1] = real(part[0].sigma_P[1]);
		sigP[2] = real(part[0].sigma_P[2]);
		sigP[3] = real(part[0].sigma_P[3]);
		sigP[4] = real(part[0].sigma_P[4]);
		sigP[5] = real(part[0].sigma_P[5]);
	  }
	  else likelihood<realnum>::Error(106,func_name,(void*) &llh->Num);			  
	  for (int i=0; i<llh->Num; i++) 
		for (int j=i; j<llh->Num; j++) 
		  im[i][j]=0;
      for (int i=0; i<llh->N_meas; i++) 
		re[i][i] = real(sigma_N[i-i%2]) + 2*sig_meas[i]*sig_meas[i];
	  re[0][1] = real(sigma_N[1]);
	  re[2][3] = real(sigma_N[3]);
	  re[0][2] = re[1][3] = real(sigma_N[4]);
	  re[0][3] = re[1][2] = real(sigma_N[5]);

      for (int i=0; i<llh->Num-llh->N_meas; i++)
		re[i+4][i+4] = sigP[i-i%2];
	  re[4][5] = sigP[1];
	  re[6][7] = sigP[3];
	  re[4][6] = re[5][7] = sigP[4];
	  re[4][7] = re[5][6] = sigP[5];

      for (int i=0; i<llh->Num-llh->N_meas; i++)
		re[i][i+4] = part[0].D[0][1]*sigP[i-i%2];
	  re[0][5] = re[1][4] = part[0].D[0][1]*sigP[1];
	  re[2][7] = re[3][6] = part[0].D[0][1]*sigP[3];
	  re[0][6] = re[1][7] = re[2][4] = re[3][5] = part[0].D[0][1]*sigP[4];
	  re[0][7] = re[1][6] = re[2][5] = re[3][4] = part[0].D[0][1]*sigP[5];

	  if ( llh->do_der_D )
	  {	
		for ( int i=0; i<llh->Num; i++ )  		
		  for ( int j=i; j<llh->Num; j++ )
			part[0].deriv_D_r[i][j] = part[0].deriv_D_i[i][j] = 0;
    	for (int i=0; i<llh->Num-llh->N_meas; i++)
		  part[0].deriv_D_r[i][i+4] = sigP[i-i%2];
		part[0].deriv_D_r[0][5] = part[1].deriv_D_r[0][4] = sigP[1];
		part[0].deriv_D_r[2][7] = part[0].deriv_D_r[3][6] = sigP[3];
		part[0].deriv_D_r[0][6] = part[0].deriv_D_r[1][7] = 
		part[0].deriv_D_r[2][4] = part[0].deriv_D_r[3][5] = sigP[4];
		part[0].deriv_D_r[0][7] = part[0].deriv_D_r[1][6] = 
		part[0].deriv_D_r[2][5] = part[0].deriv_D_r[3][4] = sigP[5];
	  }
	}
	else likelihood<realnum>::Error(106,func_name,(void*) &llh->Num);
  }
  
// set selected covariance terms to 0
  for ( int i=0; i<zero_cov_row_num; i++ )
  {
	int row2zero = zero_cov_row[i];
	if (row2zero<llh->Num) {
	  for ( int j=0; j<llh->Num; j++ )  {
		if ( j != row2zero ) {
		  re[row2zero][j] = re[j][row2zero] = 0.;
		  im[row2zero][j] = im[j][row2zero] = 0.;
		}
		else re[j][j] = 1.;
		if ( llh->do_der_D ) {
		  for (int k=0; k<N_part; k++)  {
			part[k].deriv_D_r[row2zero][j] = part[k].deriv_D_r[j][row2zero] = 0.;
			part[k].deriv_D_i[row2zero][j] = part[k].deriv_D_i[j][row2zero] = 0.;
		  }
		}
	  }
	}
	else likelihood<realnum>::Error(114,func_name);
  }
  
//copy the mixed 2. derivatives wrt terms to the lower triangle
  int k,l;
  if ( llh->do_der_D>1 )
	for (int i=0; i<N_part; i++) 
	  for (int j=i+1; j<N_part; j++) 
		for (int q=0; q<Dass[i][0][1].toNum; q++) 
		  for (int r=0; r<Dass[j][0][1].toNum; r++) 
			if ( (k=Dass[i][0][1].to(q,0))==Dass[j][0][1].to(r,0) && (l=Dass[i][0][1].to(q,1))==Dass[j][0][1].to(r,1) )
			  partpart[j][i].deriv_D_r[k][l] = partpart[i][j].deriv_D_r[k][l];
  
//  for (int i = 0;  i < llh->Num;  i++) 		for (int j = i;  j < llh->Num;  j++) 
//	if ( re[i][j] < 0. ) 
//	{ 
//	  int cov_terms[2] = { i, j }; 
//	  likelihood<realnum>::Error(201,func_name,(void*) &cov_terms);
//	  re[i][j] = 1.;
//	}
}


template <typename realnum>
void covar_matrix<realnum>::SetZeroRows( int r1=-1, int r2=-1, int r3=-1, int r4=-1, int r5=-1, int r6=-1 )
{
  int r[6] = {r1,r2,r3,r4,r5,r6};
  zero_cov_row_num=0;
  int i=0;
  while ( r[i] >= 0 ) {
	zero_cov_row_num++;
	zero_cov_row[i] = r[i];
	i++;
  }
}

template <typename realnum>
int covar_matrix<realnum>::CheckZeroRow( int row ) {
  for (int i=0; i<zero_cov_row_num; i++)
	if (zero_cov_row[i]==row) return 1;
  return 0;
}

template <typename realnum>
int covar_matrix<realnum>::GetNumZeroRows() {
  return zero_cov_row_num;
}


template <typename realnum>
void covar_matrix<realnum>::Print( )
{
  cout << "Actual cov. matrix (before inverting):" << endl;
  int wid=8; if ( !no_imag ) wid=16;
  for (int i=0; i<Num; i++)
  {
	for (int j=0; j<Num; j++) 
	  if ( j<i ) cout << setw(wid) << 0. << " ";
	  else if ( no_imag ) cout << setw(wid) << re[i][j] << " ";
	  else cout << setw(wid) << "(" << re[i][j] << "," << im[i][j] << ")";
	cout << endl;
  }
}

#undef to

}

