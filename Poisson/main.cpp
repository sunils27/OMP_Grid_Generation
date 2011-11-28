#include "Matrix.h"
#include <iostream>
#include <fstream>
#include <windows.h>
using namespace std;
#include <omp.h>
#include <cmath>
#include "Poisson.h"
double epsx, epsy;
const double eps = 1.0E-10;

//omp lock for eps updating
omp_lock_t eps_lock;
int iter;

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
void PoissonGrid( Matrix* x, Matrix* y, double* epsxy, double xlen, double ylen )
{
	Matrix* xnew = new Matrix( x->GetNumRows(), x->GetNumCols() );
	Matrix* ynew = new Matrix( y->GetNumRows(), y->GetNumCols() );
	int i,j;
	iter = 0;
	double x_temp, y_temp, g11, g12, g22, den1, den2,den3,den4,den;
	double xSubXi,xSubEta,ySubXi,ySubEta;
	double deltaXi = xlen/(double)(x->GetNumRows()-1);
	double deltaEta = ylen/(double)(x->GetNumCols()-1);
	double J, Q;
	g11 = g12 = g22 = (double)(0);
	xSubXi = xSubEta = ySubXi = ySubEta = (double)(0);
	double localx, localy = 0.0;
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
				Q = 0.0;//E-1*exp(-0.0*fabs(j*deltaXi-1.0));
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

				//cout<<x_temp<<" "<<y_temp<<" ";
			}
			//cout<<endl;
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
	delete xnew;
	delete ynew;
}


int main ( int argc, char* argv[] )
{
	int nt = atoi(argv[1]);
	Matrix* x;
	Matrix* y;
	int nx;
	int ny;
	nx = 11;
	ny = 11;
	x = new Matrix( nx, ny );
	y = new Matrix( nx, ny );
	double xlen, ylen;
	xlen = 1.0;
	ylen = 1.0;


	double dx, dy;
	dx = xlen/(double)(nx-1);
	dy = ylen/(double)(ny-1);
	//lower boundary
	for ( int i=0;i<nx;i++ )
	{
		x->SetAt( i, 0, dx*(double)(i) );
		y->SetAt( i, 0, 0 );
	}
	//left boundary
	for ( int j=0;j<ny;j++ )
	{
		x->SetAt( 0, j, 0 );
		y->SetAt( 0, j, dy*(double)j );
	}
	//right boundary
	for ( int j=0;j<ny;j++ )
	{
		x->SetAt( nx-1, j, xlen );
		y->SetAt( nx-1, j, dy*(double)j );
	}
	//upper boundary
	for ( int i=0;i<nx;i++ )
	{
		x->SetAt( i, ny-1, dx*(double)(i) );
		y->SetAt( i, ny-1, (ylen) +(double)(i)*dx/2.0 - dx*dx*2.0);//
	}
	//WriteGrid( x, y, "grid.dat" ); //use this for debugging Dirichlet boundaries

	double eps[3]; //residuals and num iterations
	//omp stuff
	omp_init_lock(&eps_lock); //initialize lock 
	omp_set_num_threads( nt ); //set num threads
	double st = omp_get_wtime(); //keep start time

	PoissonGrid( x, y, eps, xlen, ylen );
	cout<<"Time :"<<omp_get_wtime()-st<<" Eps x: "<<eps[0]<<"  Eps y: "<<eps[1]<<" Iterations: "<<iter<<endl;
	WriteGrid( x, y, "grid.dat" );


	/*Poisson* p = new Poisson( nt );
	p->SetInput( x, y );
	//p->Add( x, y, sum );
	p->IterateAndSolve( 0.001, 10000 );
	//x->DebugMatrix();
	//y->DebugMatrix();
	//sum->DebugMatrix();
	p->WriteGrid( "grid.dat" );*/



	//clean up
	delete x;
	delete y;
	//delete p;
	return 0;
}

