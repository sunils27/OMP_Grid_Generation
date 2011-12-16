#include "Matrix.h"
#include "Coefficients.h"
#include <iostream>
#include <fstream>
#include <windows.h>
using namespace std;
#include <omp.h>
#include <cmath>
#include <windows.h>
#include <process.h>
#include <gsl/gsl_vector.h>
 extern "C" {
   // Get some declarations
   #include "utils.h"
 }
const double eps = 1.0E-12;

//GLOBALS
typedef struct 
{
	Matrix* xgr;
	Matrix* ygr;
	double* eps;
	double xl, yl;
	string fName;
	//omp lock for eps updating
	omp_lock_t eps_lock;
	//run-time, wall clock
	double* times; //times[0] is st and times[1] is et
	int nthreads;
}thread_info, *ptr_tinfo;

Matrix* yvals;
double basethk;
double spacing;
int degree;
int nt; //num threads
//fin grid
Matrix* xfin;
Matrix* yfin;
//curved fluid grid
Matrix* xfluid;
Matrix* yfluid;
//lower base grid
Matrix* xlbase;
Matrix* ylbase;
Matrix* xubase;
Matrix* yubase;
Matrix* xlchan;
Matrix* ylchan;
Matrix* xuchan;
Matrix* yuchan;
thread_info t1info;
thread_info t2info;
thread_info t3info;
thread_info t4info;
thread_info t5info;
thread_info t6info;

void SetupFluidGrid( )
{
}

void readargs( int argc, char* argv[] )
{
	//order of args : numthreads basethk, spacing, y0, y1, y2, ...y(n-1)
	//for ( int i=0;i<argc;i++ )
		//cout<<i<<" :"<<argv[i]<<endl;
	{
		nt = atoi( argv[1] );
		basethk = atof( argv[2] );
		spacing = atof( argv[3] );
		int numyvals = argc - 4;
		degree = numyvals - 1;
		//assign the yvals matrix
		yvals = new Matrix( 1, numyvals );
		for ( int i=0;i<numyvals;i++ )
			yvals->SetAt( 0, i, atof(argv[i+4]) );
	}
}


//Poisson solver for elliptic grid generation
unsigned int __stdcall  PoissonGrid( void* pVoid )//( Matrix* x, Matrix* y, double* epsxy, double xlen, double ylen )
{
	double epsx, epsy;
	int iter, nt;
	//////////////////////////////unpack info
	ptr_tinfo tinfo = (ptr_tinfo)pVoid;
	Matrix* x = tinfo->xgr;
	Matrix* y = tinfo->ygr;
	double* epsxy = tinfo->eps;
	double xlen = tinfo->xl;
	double ylen = tinfo->yl;
	string fName = tinfo->fName;
	omp_lock_t eps_lock = tinfo->eps_lock;
	double* times = tinfo->times;
	nt = tinfo->nthreads;
	omp_set_num_threads( nt ); //set num threads

	//////////////////////////////end unpack
	//initialize matrices for x and y
	Matrix* xnew = new Matrix( x->GetNumRows(), x->GetNumCols() );
	Matrix* ynew = new Matrix( y->GetNumRows(), y->GetNumCols() );
	int i,j;
	iter = 0;
	//metrics
	double x_temp, y_temp, g11, g12, g22, den1, den2,den3,den4,den;
	double xSubXi,xSubEta,ySubXi,ySubEta;
	double deltaXi = xlen/(double)(x->GetNumRows()-1);
	double deltaEta = ylen/(double)(x->GetNumCols()-1);
	double J, Q;
	g11 = g12 = g22 = (double)(0);
	//derivatives
	xSubXi = xSubEta = ySubXi = ySubEta = (double)(0);
	double localx, localy = 0.0;
	times[0] = omp_get_wtime();
	//start the Jacobi algorithm
	do
	{

#pragma omp parallel for shared( x, y, iter,deltaXi,deltaEta,localx,localy ) private( i,x_temp,y_temp,g11,g12,g22,den1,den2,den3,den,J,Q,xSubXi,xSubEta,ySubXi,ySubEta,epsx,epsy )
		for ( i=1;i<x->GetNumRows()-1;i++ )
		{
			localx = localy = 0.0;
			for ( j=1;j<x->GetNumCols()-1;j++ )
			{
				localx = localy = 0.0;
				xSubXi = ( x->GetAt(i+1,j) - x->GetAt(i-1, j) )/(2.0*deltaXi);
				xSubEta = ( x->GetAt(i, j+1) - x->GetAt(i, j-1) )/(2.0*deltaEta);
				ySubXi = ( y->GetAt(i+1, j) - y->GetAt(i-1, j) )/(2.0*deltaXi);
				ySubEta = ( y->GetAt(i, j+1) - y->GetAt(i, j-1) )/(2.0*deltaEta);

				g11 = xSubEta*xSubEta + ySubEta*ySubEta;   
				g12 = xSubXi*xSubEta + ySubXi*ySubEta;
				g22 = xSubXi*xSubXi + ySubXi*ySubXi;
				J = xSubXi*ySubEta - xSubEta*ySubXi;
				Q = 0.0E-1;//*exp(-1.0*fabs(j*deltaXi-1.0));
				//-----------------denominators----------------------
				den1 = g11/(deltaXi*deltaXi);
				den2 = g12/(2.0*deltaXi*deltaEta);
				den3 = g22/(deltaEta*deltaEta);
				den4 = J*J*Q;
				den = 2.0*g11/(deltaXi*deltaXi)+2.0*g22/(deltaEta*deltaEta);
				x_temp = ( den1*(x->GetAt(i+1, j ) + 
					x->GetAt(i-1, j ) )
					-den2*( x->GetAt(i+1, j+1 ) - 
					x->GetAt(i+1, j-1) 
					- x->GetAt(i-1, j+1) + 
					x->GetAt(i-1, j-1) )
					+ den3*( x->GetAt(i, j+1 ) + 
					x->GetAt(i, j-1 ) ) 
					+ den4*ySubEta )/den;
				xnew->SetAt( i,j,x_temp );

				y_temp = ( den1*( y->GetAt(i+1, j ) + 
					y->GetAt(i-1, j ) )
					-den2*( y->GetAt(i+1, j+1 ) - 
					y->GetAt(i+1, j-1) 
					- y->GetAt(i-1, j+1) + 
					y->GetAt(i-1, j-1) )
					+ den3*( y->GetAt(i, j+1 ) + 
					y->GetAt(i, j-1 ) ) 
					+ den4*ySubEta )/den;

				ynew->SetAt( i,j,y_temp );
				epsx = fabs( x_temp - (x->GetAt(i,j)) );
				epsy = fabs( y_temp - (y->GetAt(i,j)) );
			}
			//update the residual under a lock, after each j station is computed
			//by a omp thread
			omp_set_lock( &eps_lock );
			{
				if ( epsx > localx )
					localx = epsx;
				if ( epsy > localy )
					localy = epsy;
			}
			omp_unset_lock( &eps_lock );
		}
		//update 
		for ( int p=1;p<x->GetNumRows()-1;p++ )
		{
			for ( int q=1;q<x->GetNumCols()-1;q++ )
			{
				x->SetAt(p,q, xnew->GetAt(p,q) );
				y->SetAt(p,q, ynew->GetAt(p,q) );
			}
		}
#pragma omp atomic
		iter++;
	}
	while( (localx>=eps) || ( localy>=eps) );
	epsxy[0] = localx;
	epsxy[1]= localy;
	epsxy[2] = iter;
	times[1] = omp_get_wtime();

	delete xnew;
	delete ynew;

	return 0;
}

int main ( int argc, char* argv[] )
{
	//some definitions before we start
	double xlen, ylen, xlen1, ylen1;
	xlen = 0.75;
	ylen = 0.5;
	xlen1 = 0.75;
	ylen1 = 0.5;
	double xmin, ymin, xmin1, ymin1 = 0.0;
	xmin = 0;
	ymin = 0;
	xmin1 = 0.0;
	//read cmd line args
	readargs( argc, argv );
	int nx, ny;
	nx = 101;
	ny = 41;
	//allocate 0th
	xfin = new Matrix( nx, ny );
	yfin = new Matrix( nx, ny );
	double times0[2];
	double times1[2];
	double times2[2];
	double times3[2];
	double times4[2];
	double times5[2];

	double eps0[3]; //residuals and num iterations
	double eps1[3]; 
	double eps2[3];
	double eps3[3];
	double eps4[3];
	double eps5[3];
	//

	double dx, dy, dx1, dy1;
	/*gsl_vector* yvalsvect = gsl_vector_calloc( 3 );
	for ( size_t i=0;i<yvalsvect->size;i++ )
	{
		gsl_vector_set( yvalsvect, i, yvals->GetAt(0,i) );
	}*/

	//generate coefficients
	//yvals->DebugMatrix( );
	Coefficients c( degree, yvals, xlen );
	Matrix* coeffs;
	c.CalculateCoeffs( );
	//gsl_vector* coeffvect = c.getCoeffs( );
	coeffs = c.getCoeffs2( );
	dx = xlen/(double)(nx-1);
	Matrix* xstn = new Matrix( nx, 1 );
	Matrix* ystn;
	for ( int i=0;i<nx;i++ )
		xstn->SetAt( i, 0, dx*(double)(i) );
	//xstn->DebugMatrix( );
	ystn = ComputeFunction( coeffs, xstn );
	//ystn->DebugMatrix( );
	//dy = ystn->GetAt(0,0)/(double)(ny-1);
	//cout<<ystn->GetAt(0,0)<<endl;
	//cout<<"coeffs :"<<endl;
	//coeffs->DebugMatrix();
	//cout<<"coeffs vector :"<<endl;
	//PrintVector( coeffvect );

	///////////////////////////////////////////////////////////////////////////////////////////
	//                                FIN GRID          ///////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////
	//lower boundary
	for ( int i=0;i<nx;i++ )
	{
		xfin->SetAt( i, 0, xmin + dx*(double)(i) );
		yfin->SetAt( i, 0, ymin + 0 );
	}
	//left boundary
	dy = ystn->GetAt(0,0)/(double)(ny-1);
	for ( int j=0;j<ny;j++ )
	{
		xfin->SetAt( 0, j, xmin + 0 );
		yfin->SetAt( 0, j, ymin + (ystn->GetAt(0,0)/(double)(ny-1))*(double)j );
	}
	//right boundary
	for ( int j=0;j<ny;j++ )
	{
		xfin->SetAt( nx-1, j, xmin + xlen );
		yfin->SetAt( nx-1, j, ymin + ( ystn->GetAt(nx-1,0)/(double)(ny-1) )*(double)j );
	}
	//upper boundary
	for ( int i=0;i<nx;i++ )
	{
		xfin->SetAt( i, ny-1, xmin + dx*(double)(i) );
		yfin->SetAt( i, ny-1, ymin + ystn->GetAt(i,0));//
	}
	//WriteGrid( x, y, "grid.dat" ); //use this for debugging Dirichlet boundaries

	//lock 0 //omp stuff
	omp_lock_t lock0;
	omp_init_lock(&lock0); //initialize lock 
	t1info.eps = eps0;
	t1info.xgr = xfin;
	t1info.ygr = yfin;
	t1info.xl = xlen;
	t1info.yl = ylen;
	t1info.fName = "gridfin.dat";
	t1info.eps_lock = lock0;
	t1info.times = times0;
	t1info.nthreads = 4;
	///////////////////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////////////////////
	//                                FLUID GRID                ///////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////
	dx1 = xlen1/(double)(nx-1);
	dy1 = ylen1/(double)(ny-1);
	//allocate 1st
	xfluid = new Matrix( nx, ny );
	yfluid = new Matrix( nx, ny );
	ymin1 = ystn->GetAt(0,0);
	//lower boundary
	for ( int i=0;i<nx;i++ )
	{
		xfluid->SetAt( i, 0, dx1*(double)(i) );
		yfluid->SetAt( i, 0, ystn->GetAt(i,0));//
	}
	//left boundary
	for ( int j=0;j<ny;j++ )
	{
		xfluid->SetAt( 0, j, 0 );
		yfluid->SetAt( 0, j, ymin1 + ((spacing-ystn->GetAt(0,0))/(double)(ny-1)) * (double)(j) );
	}
	//right boundary
	for ( int j=0;j<ny;j++ )
	{
		xfluid->SetAt( nx-1, j, xlen1 );
		yfluid->SetAt( nx-1, j, ystn->GetAt(nx-1,0) + ((spacing-ystn->GetAt(nx-1,0))/(double)(ny-1)) * (double)(j) );
	}
	//upper boundary
	for ( int i=0;i<nx;i++ )
	{
		xfluid->SetAt( i, ny-1, dx1*(double)(i) );
		yfluid->SetAt( i, ny-1, spacing );//
	}
	//WriteGrid( xfluid, yfluid, "grid2.dat" );

	omp_lock_t lock1;
	omp_init_lock(&lock1); //initialize lock 
	t2info.eps = eps1;
	t2info.xgr = xfluid;
	t2info.ygr = yfluid;
	t2info.xl = xlen1;
	t2info.yl = ylen1;
	t2info.fName = "gridfluid.dat";
	t2info.eps_lock = lock1;
	t2info.times = times1;
	t2info.nthreads = 4;
	///////////////////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////////////////////
	//                                LOWER BASE GRID          ////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////
	//allocate memory
	int nxsmall, nysmall;
	nxsmall = 21;
	nysmall = 41;
	xlbase = new Matrix( nxsmall, nysmall );
	ylbase = new Matrix( nxsmall, nysmall );
	//lower boundary
	for ( int i=0;i<nx;i++ )
	{
		xlbase->SetAt( i, 0, -( basethk/(double)(nxsmall-1) )*(double)(nxsmall-1-i) );
		ylbase->SetAt( i, 0, 0 );
	}
	//left boundary
	//dy = ystn->GetAt(0,0)/(double)(ny-1);
	for ( int j=0;j<ny;j++ )
	{
		xlbase->SetAt( 0, j, -basethk );
		ylbase->SetAt( 0, j, (ystn->GetAt(0,0)/(double)(nysmall-1))*(double)(j) );//ymin + (ystn->GetAt(0,0)/(double)(ny-1))*(double)j );
	}
	//right boundary
	for ( int j=0;j<ny;j++ )
	{
		xlbase->SetAt( nxsmall-1, j, 0 );
		ylbase->SetAt( nxsmall-1, j, ( ystn->GetAt(0,0)/(double)(nysmall-1) )*(double)j );
	}
	//upper boundary
	for ( int i=0;i<nx;i++ )
	{
		xlbase->SetAt( i, nysmall-1, -( basethk/(double)(nxsmall-1) )*(double)(nxsmall-1-i) );
		ylbase->SetAt( i, nysmall-1, ystn->GetAt(0,0));//
	}

	omp_lock_t lock2;
	omp_init_lock(&lock2); //initialize lock 
	t3info.eps = eps2;
	t3info.xgr = xlbase;
	t3info.ygr = ylbase;
	t3info.xl = basethk; //need to change these
	t3info.yl = 0.5;
	t3info.fName = "gridlbase.dat";
	t3info.eps_lock = lock2;
	t3info.times = times2;
	t3info.nthreads = 2;
	//PoissonGrid( (void*)&t3info );
	//WriteGrid( xlbase, ylbase, "gridlbase.dat" );
	///////////////////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////////////////////
	//                                UPPER BASE GRID          ////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////
	xubase = new Matrix( nxsmall, nysmall );
	yubase = new Matrix( nxsmall, nysmall );
	//lower boundary
	for ( int i=0;i<nx;i++ )
	{
		xubase->SetAt( i, 0, -( basethk/(double)(nxsmall-1) )*(double)(nxsmall-1-i) );
		yubase->SetAt( i, 0, ystn->GetAt(0,0) );
	}
	//left boundary
	//dy = ystn->GetAt(0,0)/(double)(ny-1);
	for ( int j=0;j<ny;j++ )
	{
		xubase->SetAt( 0, j, -basethk );
		yubase->SetAt( 0, j, ystn->GetAt(0,0) + ((spacing-ystn->GetAt(0,0))/(double)(nysmall-1))*(double)(j) );//ymin + (ystn->GetAt(0,0)/(double)(ny-1))*(double)j );
	}
	//right boundary
	for ( int j=0;j<ny;j++ )
	{
		xubase->SetAt( nxsmall-1, j, 0 );
		yubase->SetAt( nxsmall-1, j, ystn->GetAt(0,0) + ( (spacing-ystn->GetAt(0,0))/(double)(nysmall-1) )*(double)j );
	}
	//upper boundary
	for ( int i=0;i<nx;i++ )
	{
		xubase->SetAt( i, nysmall-1, -( basethk/(double)(nxsmall-1) )*(double)(nxsmall-1-i) );
		yubase->SetAt( i, nysmall-1, spacing);//
	}

	//WriteGrid( xubase, yubase, "gridubase.dat" );
	omp_lock_t lock3;
	omp_init_lock(&lock3); //initialize lock 
	t4info.eps = eps3;
	t4info.xgr = xubase;
	t4info.ygr = yubase;
	t4info.xl = basethk; //need to change these
	t4info.yl = 0.5;
	t4info.fName = "gridubase.dat";
	t4info.eps_lock = lock3;
	t4info.times = times3;
	t4info.nthreads = 2;
	///////////////////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////////////////////
	//                                LOWER CHANNEL FLUID      ////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////
	xlchan = new Matrix( nxsmall, nysmall );
	ylchan = new Matrix( nxsmall, nysmall );
	//ystn->DebugMatrix();
	//cout<<ystn->GetAt(nx-1,0)<<endl;
	//lower boundary
	for ( int i=0;i<nx;i++ )
	{
		xlchan->SetAt( i, 0, xlen + ( (1.0-xlen)/(double)(nxsmall-1) )*(double)(i) );
		ylchan->SetAt( i, 0, 0 );
	}
	//left boundary
	//dy = ystn->GetAt(0,0)/(double)(ny-1);
	for ( int j=0;j<ny;j++ )
	{
		xlchan->SetAt( 0, j, xlen );
		ylchan->SetAt( 0, j, ((ystn->GetAt(nx-1,0))/(double)(nysmall-1))*(double)(j) );//ymin + (ystn->GetAt(0,0)/(double)(ny-1))*(double)j );
	}
	//right boundary
	for ( int j=0;j<ny;j++ )
	{
		xlchan->SetAt( nxsmall-1, j, 1.0 );
		ylchan->SetAt( nxsmall-1, j, ( (ystn->GetAt(nx-1,0))/(double)(nysmall-1) )*(double)j );
	}
	//upper boundary
	for ( int i=0;i<nx;i++ )
	{
		xlchan->SetAt( i, nysmall-1, xlen + ( (1.0-xlen)/(double)(nxsmall-1) )*(double)(i) );
		ylchan->SetAt( i, nysmall-1, ystn->GetAt(nx-1,0));//
	}
	//WriteGrid( xlchan, ylchan, "gridlchan.dat" );
	omp_lock_t lock4;
	omp_init_lock(&lock4); //initialize lock 
	t5info.eps = eps4;
	t5info.xgr = xlchan;
	t5info.ygr = ylchan;
	t5info.xl = basethk; //need to change these
	t5info.yl = 0.5;
	t5info.fName = "gridlchan.dat";
	t5info.eps_lock = lock4;
	t5info.times = times4;
	t5info.nthreads = 2;
	///////////////////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////////////////////
	//                                UPPER CHANNEL FLUID      ////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////
	xuchan = new Matrix( nxsmall, nysmall );
	yuchan = new Matrix( nxsmall, nysmall );
	//ystn->DebugMatrix();
	//cout<<ystn->GetAt(nx-1,0)<<endl;
	//lower boundary
	for ( int i=0;i<nx;i++ )
	{
		xuchan->SetAt( i, 0, xlen + ( (1.0-xlen)/(double)(nxsmall-1) )*(double)(i) );
		yuchan->SetAt( i, 0, ystn->GetAt(nx-1,0) );
	}
	//left boundary
	//dy = ystn->GetAt(0,0)/(double)(ny-1);
	for ( int j=0;j<ny;j++ )
	{
		xuchan->SetAt( 0, j, xlen );
		yuchan->SetAt( 0, j, ystn->GetAt(nx-1,0) + ((spacing-ystn->GetAt(nx-1,0))/(double)(nysmall-1))*(double)(j) );//ymin + (ystn->GetAt(0,0)/(double)(ny-1))*(double)j );
	}
	//right boundary
	for ( int j=0;j<ny;j++ )
	{
		xuchan->SetAt( nxsmall-1, j, 1.0 );
		yuchan->SetAt( nxsmall-1, j, ystn->GetAt(nx-1,0) + ( (spacing-ystn->GetAt(nx-1,0))/(double)(nysmall-1) )*(double)j );
	}
	//upper boundary
	for ( int i=0;i<nx;i++ )
	{
		xuchan->SetAt( i, nysmall-1, xlen + ( (1.0-xlen)/(double)(nxsmall-1) )*(double)(i) );
		yuchan->SetAt( i, nysmall-1, spacing);//
	}
	//WriteGrid( xuchan, yuchan, "griduchan.dat" );
	omp_lock_t lock5;
	omp_init_lock(&lock5); //initialize lock 
	t6info.eps = eps5;
	t6info.xgr = xuchan;
	t6info.ygr = yuchan;
	t6info.xl = basethk; //need to change these
	t6info.yl = 0.5;
	t6info.fName = "griduchan.dat";
	t6info.eps_lock = lock5;
	t6info.times = times5;
	t6info.nthreads = 2;


	HANDLE t[6];
	t[0] = (HANDLE)::_beginthreadex( NULL, 0, PoissonGrid, (void*)&t1info, 0, NULL );
	t[1] = (HANDLE)::_beginthreadex( NULL, 0, PoissonGrid, (void*)&t2info, 0, NULL );
	t[2] = (HANDLE)::_beginthreadex( NULL, 0, PoissonGrid, (void*)&t3info, 0, NULL );
	t[3] = (HANDLE)::_beginthreadex( NULL, 0, PoissonGrid, (void*)&t4info, 0, NULL );
	t[4] = (HANDLE)::_beginthreadex( NULL, 0, PoissonGrid, (void*)&t5info, 0, NULL );
	t[5] = (HANDLE)::_beginthreadex( NULL, 0, PoissonGrid, (void*)&t6info, 0, NULL );
	::WaitForMultipleObjects( 6, t, true, INFINITE );

	//PoissonGrid( (void*)&t2info );
	cout<<"Time :"<<times0[1]-times0[0]<<" Eps x: "<<eps0[0]<<"  Eps y: "<<eps0[1]<<" Iterations: "<<eps0[2]<<endl;
	cout<<"Time :"<<times1[1]-times1[0]<<" Eps x: "<<eps1[0]<<"  Eps y: "<<eps1[1]<<" Iterations: "<<eps1[2]<<endl;
	cout<<"Time :"<<times2[1]-times2[0]<<" Eps x: "<<eps2[0]<<"  Eps y: "<<eps2[1]<<" Iterations: "<<eps2[2]<<endl;
	cout<<"Time :"<<times3[1]-times3[0]<<" Eps x: "<<eps3[0]<<"  Eps y: "<<eps3[1]<<" Iterations: "<<eps3[2]<<endl;
	cout<<"Time :"<<times4[1]-times4[0]<<" Eps x: "<<eps4[0]<<"  Eps y: "<<eps4[1]<<" Iterations: "<<eps4[2]<<endl;
	cout<<"Time :"<<times5[1]-times5[0]<<" Eps x: "<<eps5[0]<<"  Eps y: "<<eps5[1]<<" Iterations: "<<eps5[2]<<endl;

	WriteGrid( t1info.xgr, t1info.ygr, t1info.fName );
	WriteGrid( t2info.xgr, t2info.ygr, t2info.fName );
	WriteGrid( t3info.xgr, t3info.ygr, t3info.fName );
	WriteGrid( t4info.xgr, t4info.ygr, t4info.fName );
	WriteGrid( t5info.xgr, t5info.ygr, t5info.fName );
	WriteGrid( t6info.xgr, t6info.ygr, t6info.fName );

	//clean up
	delete xfin;
	delete yfin;
	delete xfluid;
	delete yfluid;
	delete xlbase;
	delete ylbase;
	delete xubase;
	delete yubase;
	delete xlchan;
	delete ylchan;
	delete xuchan;
	delete yuchan;

	delete yvals;

	//delete p;
	return 0;
}
