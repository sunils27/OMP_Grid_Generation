#include "Coefficients.h"
////////////////////////////////////////////////////////////////////////////////
Coefficients::Coefficients( int degreePoly, gsl_vector* yvals, double xl )
{
	polyDegree = degreePoly;
	xLength = xl;
	//cout<<"Degree of polynomial function :"<<polyDegree<<endl;
	//cout<<"Length of domain in x-direction :"<<xLength<<endl;
	//allocate vector for y values
	yValues = gsl_vector_calloc( polyDegree + 1 );
	for ( int i=0;i<polyDegree+1;i++ )
	{
		gsl_vector_set( yValues, i, gsl_vector_get( yvals, i ) );
		//cout<<gsl_vector_get( yValues, i )<<endl;
	}
	//allocate vector for coefficients to be calculated
	polyCoeffs = gsl_vector_calloc( polyDegree + 1 );
	//allocate vector for x values
	xValues = gsl_vector_calloc( polyDegree + 1 );
	//allocate memory for matrix
	matrix = gsl_matrix_calloc( polyDegree + 1, polyDegree + 1 );
}
////////////////////////////////////////////////////////////////////////////////
Coefficients::Coefficients( int degreePoly, Matrix* yvals, double xl )
{
	polyDegree = degreePoly;
	xLength = xl;
	//allocate vector for y values
	yValues = gsl_vector_calloc( polyDegree + 1 );
	for ( int i=0;i<polyDegree+1;i++ )
	{
		gsl_vector_set( yValues, i, yvals->GetAt( 0, i ) );
		//cout<<gsl_vector_get( yValues, i )<<endl;
	}
	//allocate vector for coefficients to be calculated
	polyCoeffs = gsl_vector_calloc( polyDegree + 1 );
	//allocate vector for x values
	xValues = gsl_vector_calloc( polyDegree + 1 );
	//allocate memory for matrix
	matrix = gsl_matrix_calloc( polyDegree + 1, polyDegree + 1 );
}
////////////////////////////////////////////////////////////////////////////////
Coefficients::~Coefficients( )
{
	gsl_vector_free( xValues );
	gsl_vector_free( yValues );
	gsl_vector_free( polyCoeffs );
	gsl_matrix_free( matrix );
}
////////////////////////////////////////////////////////////////////////////////
//RETURN THE DEGREE OF THE POLYNOMIAL
int Coefficients::getDegree( )
{
	return polyDegree;
}
////////////////////////////////////////////////////////////////////////////////
//RETURN THE COEFFICIENTS OF THE POLYNOMIAL
gsl_vector* Coefficients::getCoeffs( )
{
	return polyCoeffs;
}
////////////////////////////////////////////////////////////////////////////////
//RETURN THE Y VALUES INPUT TO THE CLASS
gsl_vector* Coefficients::get_Y_Vals( )
{
	return yValues;
}

////////////////////////////////////////////////////////////////////////////////
//CALCULATE THE COEFFICIENTS OF THE POLYNOMIAL
void Coefficients::CalculateCoeffs( )
{
	double x;
	double val;
	gsl_permutation* p = gsl_permutation_alloc ( polyDegree+1 );
	int s;
	//set the elements of the matrix
	for ( int i=0;i<polyDegree+1;i++ )
	{
		x = (double)(i)*(xLength)/(double)(polyDegree);
		for ( int j=0;j<polyDegree+1;j++ )
		{
			val = pow( x, (double)(polyDegree-j) );
			gsl_matrix_set( matrix, i, j, val );
		}
	}
	gsl_linalg_LU_decomp ( matrix, p, &s );
	gsl_linalg_LU_solve ( matrix, p, yValues, polyCoeffs );
	gsl_permutation_free ( p );
}
////////////////////////////////////////////////////////////////////////////////
void Coefficients::ShowMatrix( )	//use this for de-bugging only
{
	cout<<"Matrix is :"<<endl;
	for ( int i=0;i<polyDegree+1;i++ )
	{
		for ( int j=0;j<polyDegree+1;j++ )
		{
			cout<<gsl_matrix_get( matrix, i, j )<<" ";
		}
		cout<<endl;
	}
}
////////////////////////////////////////////////////////////////////////////////
Matrix* Coefficients::getCoeffs2( )
{
	Matrix* matcoeffs = new Matrix( polyCoeffs->size, 1 );
	for ( size_t i=0;i<polyCoeffs->size;i++ )
		matcoeffs->SetAt( i, 0, gsl_vector_get( polyCoeffs, i ) );
	return matcoeffs;
}