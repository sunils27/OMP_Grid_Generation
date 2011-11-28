#ifndef MATRIX_H
#define MATRIX_H
#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
using namespace std;
class Matrix
{
private:
	double* matrix;
	int size;
	int nR;
	int nC;
public:
	Matrix(int nRows, int nCols);
	int GetAt( int r, int c, double& val );
	double GetAt( int r, int c );

	int SetAt( int r, int c, double val );
	double* GetMatrix( );
	int SetMatrix( double* mat, int numR, int numC );
	void DebugMatrix( );
	void WriteMatrix( string fileName );
	int GetNumRows( );
	int GetNumCols( );
	~Matrix(void);
};

#endif //MATRIX_H