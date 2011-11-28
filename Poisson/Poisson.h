#pragma once
#include <cmath>
#include <omp.h>
#include "Matrix.h"
#include <iostream>
using namespace std;
class Poisson
{
private:
	int numprocs;
	int numthreads;
	Matrix* xpvt;
	Matrix* ypvt;
	double tolx, toly, tol;

	void ComputeMetrics( int, int );
      //variables in computational domain
      //epsilon
      double eps_x;
      double eps_y;
      double maxX;
	  double maxY;

public:
	Poisson(int nt);
	void SetInput( Matrix*, Matrix* );
	int Solver( );
	void IterateAndSolve( double tol, int max_iter );
	void Add( Matrix*, Matrix*, Matrix* );
	void WriteGrid( string fName );
	void SolverSimple( Matrix* );
	~Poisson(void);
};
