#include "Matrix.h"
////////////////////////////////////////////////////////////////////////////////
Matrix::Matrix(int nRows, int nCols)
{
	srand( 8678764 );
	nR = nRows;
	nC = nCols;
	size = nR*nC;
	matrix = new double[size];
	//zero everything out
	for ( int i=0;i<size;i++ )
		matrix[i] = (double)rand()/(double)RAND_MAX;
}
////////////////////////////////////////////////////////////////////////////////
Matrix::~Matrix(void)
{
	delete [] matrix;
}
////////////////////////////////////////////////////////////////////////////////
int Matrix::GetAt(int r, int c, double& val)
{
	//TODO
	//check limits
	int index = -1;
	index = c+r*nC;
	if ( (index<0) || (index>=size) )
		return -1;
	else
	{
		val = matrix[index];
		return 0;		
	}
}
////////////////////////////////////////////////////////////////////////////////
double Matrix::GetAt( int r, int c )
{
	int index = -1;
	index = c+r*nC;
	return matrix[index];
}
////////////////////////////////////////////////////////////////////////////////
int Matrix::SetAt( int r, int c, double val )
{
	// caller should check for return code
	int index = -1;
	index = c+r*nC;
	if ( (index <0) || (index >= size) )
		return -1;
	else
		matrix[index] = val;
	return 0;
}
////////////////////////////////////////////////////////////////////////////////
double* Matrix::GetMatrix( )
{
	return matrix;
}
////////////////////////////////////////////////////////////////////////////////
void Matrix::DebugMatrix( )
{
	double val;
	for ( int i=0;i<nR;i++ )
	{
		for ( int j=0;j<nC;j++ )
		{
			this->GetAt(i,j,val);
			cout<<val<<" ";
		}
		cout<<endl;
	}
}
////////////////////////////////////////////////////////////////////////////////
int Matrix::SetMatrix( double* mat, int numR, int numC )
{

	//TODO: add error checking for bounds
	int counter = 0;
	for ( int r=0;r<numR;r++ )
	{
		for ( int c=0;c<numC;c++ )
		{
			this->SetAt( r,c,mat[counter] );
			counter++;
		}
	}
	return 0;
}
////////////////////////////////////////////////////////////////////////////////
void Matrix::WriteMatrix( string fileName )
{
	ofstream fout;
	fout.open( fileName.c_str() );
	double val;
	for ( int i=0;i<nR;i++ )
	{
		for ( int j=0;j<nC;j++ )
		{
			this->GetAt(i,j,val);
			fout<<val<<" ";
		}
		fout<<endl;
	}
	fout.close();
}
////////////////////////////////////////////////////////////////////////////////
int Matrix::GetNumCols( )
{
	return nC;
}
////////////////////////////////////////////////////////////////////////////////
int Matrix::GetNumRows( )
{
	return nR;
}
////////////////////////////////////////////////////////////////////////////////
