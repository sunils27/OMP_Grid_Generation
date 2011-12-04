#include "Matrix.h"
#include <iostream>
#include <fstream>
#include <windows.h>
using namespace std;
#include <omp.h>
#include <cmath>
#include <windows.h>
#include <process.h>
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
}thread_info, *ptr_tinfo;

void WriteGrid( Matrix* xgrid, Matrix* ygrid, string fName )
{
	cout<<"Writing grid to "<<fName<<" ...";
	ofstream fout;

	fout.open( fName.c_str() );
	fout<<"TITLE = "<<'"'<<"Solid Grid"<<'"'<<endl;
	fout<<"VARIABLES = "<<'"'<<'x'<<'"'<<"\t"<<'"'<<'y'<<'"'<<"\t"<<endl;
	fout<<"ZONE "<<"I= "<<xgrid->GetNumCols()<<"\t J= "<<xgrid->GetNumRows()<<endl;
	fout<<"DT = ( SINGLE SINGLE )"<<endl;
	for ( int j=0;j<xgrid->GetNumRows();j++ )
	{
		for ( int i=0;i<xgrid->GetNumCols();i++ )
		{

			fout<<xgrid->GetAt( j, i )<<" "<<
				ygrid->GetAt( j, i )<<endl;
		}
	}
	cout<<" Done "<<endl;
	cout<<"============================================================="<<endl;
	fout.close();
}
//Poisson solver for elliptic grid generation
unsigned int __stdcall  PoissonGrid( void* pVoid )//( Matrix* x, Matrix* y, double* epsxy, double xlen, double ylen )
{
	double epsx, epsy;
	int iter;
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
	thread_info t1info;
	thread_info t2info;
	int nt = atoi(argv[1]);
	int nx;
	int ny;
	nx = 101;
	ny = 101;
	//0th grid
	Matrix* x;
	Matrix* y;
	//allocate 0th
	x = new Matrix( nx, ny );
	y = new Matrix( nx, ny );
	double times0[2];
	//1st grid
	Matrix* x1;
	Matrix* y1;
	//allocate 1st
	x1 = new Matrix( nx, ny );
	y1 = new Matrix( nx, ny );
	double times1[2];
	double xlen, ylen, xlen1, ylen1;
	xlen = 0.2;
	ylen = 0.1;
	xlen1 = 0.1;
	ylen1 = 0.1;
	double xmin, ymin, xmin1, ymin1 = 0.0;
	xmin = 0;
	ymin = 0;
	xmin1 = 0.2;



	double dx, dy, dx1, dy1;
	dx = xlen/(double)(nx-1);
	dy = ylen/(double)(ny-1);
	//lower boundary
	for ( int i=0;i<nx;i++ )
	{
		x->SetAt( i, 0, xmin + dx*(double)(i) );
		y->SetAt( i, 0, ymin + 0 );
	}
	//left boundary
	for ( int j=0;j<ny;j++ )
	{
		x->SetAt( 0, j, xmin + 0 );
		y->SetAt( 0, j, ymin + dy*(double)j );
	}
	//right boundary
	for ( int j=0;j<ny;j++ )
	{
		x->SetAt( nx-1, j, xmin + xlen );
		y->SetAt( nx-1, j, ymin + dy*(double)j );
	}
	//upper boundary
	for ( int i=0;i<nx;i++ )
	{
		x->SetAt( i, ny-1, xmin + dx*(double)(i) );
		y->SetAt( i, ny-1, ymin + (ylen) +(double)(i)*dx/2.0 - dx*dx*2.0);//
	}
	//WriteGrid( x, y, "grid.dat" ); //use this for debugging Dirichlet boundaries

	double eps0[3]; //residuals and num iterations
	double eps1[3]; 
	omp_set_num_threads( nt ); //set num threads
	//lock 0 //omp stuff
	omp_lock_t lock0;
	omp_init_lock(&lock0); //initialize lock 
	t1info.eps = eps0;
	t1info.xgr = x;
	t1info.ygr = y;
	t1info.xl = xlen;
	t1info.yl = ylen;
	t1info.fName = "grid1.dat";
	t1info.eps_lock = lock0;
	t1info.times = times0;

	dx1 = xlen1/(double)(nx-1);
	dy1 = ylen1/(double)(ny-1);
	//lower boundary
	for ( int i=0;i<nx;i++ )
	{
		x1->SetAt( i, 0, xmin1 + dx1*(double)(i) );
		y1->SetAt( i, 0, ymin1 + 0 );
	}
	//left boundary
	for ( int j=0;j<ny;j++ )
	{
		x1->SetAt( 0, j, xmin1 + 0 );
		y1->SetAt( 0, j, ymin1 + dy1*(double)j );
	}
	//right boundary
	for ( int j=0;j<ny;j++ )
	{
		x1->SetAt( nx-1, j, xmin1 + xlen1 );
		y1->SetAt( nx-1, j, ymin1 + dy1*(double)j );
	}
	//upper boundary
	for ( int i=0;i<nx;i++ )
	{
		x1->SetAt( i, ny-1, xmin1 + dx1*(double)(i) );
		y1->SetAt( i, ny-1, ymin1 + (ylen1) +(double)(i)*dx1/2.0 - dx1*dx1*2.0);//
	}
	//WriteGrid( x1, y1, "grid2.dat" );

	omp_lock_t lock1;
	omp_init_lock(&lock1); //initialize lock 
	t2info.eps = eps1;
	t2info.xgr = x1;
	t2info.ygr = y1;
	t2info.xl = xlen1;
	t2info.yl = ylen1;
	t2info.fName = "grid2.dat";
	t2info.eps_lock = lock1;
	t2info.times = times1;

	HANDLE t[2];
	t[0] = (HANDLE)::_beginthreadex( NULL, 0, PoissonGrid, (void*)&t1info, 0, NULL );
	t[1] = (HANDLE)::_beginthreadex( NULL, 0, PoissonGrid, (void*)&t2info, 0, NULL );
	::WaitForMultipleObjects( 2, t, true, INFINITE );

	//PoissonGrid( (void*)&t2info );

	//PoissonGrid( x, y, eps, xlen, ylen );

	cout<<"Time :"<<times0[1]-times0[0]<<" Eps x: "<<eps0[0]<<"  Eps y: "<<eps0[1]<<" Iterations: "<<eps0[2]<<endl;
	cout<<"Time :"<<times1[1]-times1[0]<<" Eps x: "<<eps1[0]<<"  Eps y: "<<eps1[1]<<" Iterations: "<<eps1[2]<<endl;
	WriteGrid( t1info.xgr, t1info.ygr, t1info.fName );
	WriteGrid( t2info.xgr, t2info.ygr, t2info.fName );


	//clean up
	delete x;
	delete y;
	delete x1;
	delete y1;
	//delete p;
	return 0;
}

