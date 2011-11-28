#include "Poisson.h"
////////////////////////////////////////////////////////////////////////////////
int Poisson::Solver( )
{
	double x_temp;
	double y_temp;
	double xVal;
	double yVal;

	eps_x, eps_y = (double)(1);
	maxX = maxY = (double)1;
	int nX;
	int nY;
	nX = xpvt->GetNumRows();
	nY = xpvt->GetNumCols();
	int i;
	int j;
	int limit;
	double deltaXi;
	double deltaEta;

	double g11;
	double g12;
	double g22;
	double xSubXi;
	double xSubEta;
	double ySubXi;
	double ySubEta;
	double J;   //Jacobian
	double Q;   //grid clustering function in the eta direction
	//variables for denominators
	double den1;
	double den2;
	double den3;
	double den4;
	double den;
	Matrix* x = xpvt;
	Matrix* y = ypvt;
	double st = omp_get_wtime();
	double eps = 1.0E-5;
	double epsx, epsy = 0;
	xSubXi,xSubEta,ySubXi,ySubEta,g11,g12,g22,J,Q,den1,den2,den3,den4,den,xVal,yVal = 0;
	deltaXi = (double)(1)/(xpvt->GetNumRows()-1);
	deltaEta = (double)(1)/(xpvt->GetNumCols()-1);
	int iter = 0;

	for ( iter=0;iter<10;iter++ )
	{
		eps = 0.0;
//#pragma omp parallel for shared(x,y,eps) private(i,x_temp,y_temp,deltaXi,deltaEta,xSubXi,xSubEta,ySubXi,ySubEta,g11,g12,g22,J,Q,den1,den2,den3,den4,den,xVal,yVal)
		for ( i=1;i<nX-1;i++ )
		{
			for ( j=1;j<nY-1;j++ )
			{
				x->GetAt( i, j, x_temp );
				y->GetAt( i, j, y_temp );
				g11 = g12 = g22 = (double)(0);
				xSubXi = xSubEta = ySubXi = ySubEta = (double)(0);

				xSubXi = ( x->GetAt(i+1,j) - 
					x->GetAt(i-1, j) )/(2.0*deltaXi);
				xSubEta = ( x->GetAt(i, j+1) - 
					x->GetAt(i, j-1) )/(2.0*deltaEta);
				ySubXi = ( y->GetAt(i+1, j) - 
					y->GetAt(i-1, j) )/(2.0*deltaXi);
				ySubEta = ( y->GetAt(i, j+1) - 
					y->GetAt(i, j-1) )/(2.0*deltaEta);

				g11 = xSubEta*xSubEta + ySubEta*ySubEta;   
				g12 = xSubXi*xSubEta + ySubXi*ySubEta;
				g22 = xSubXi*xSubXi + ySubXi*ySubXi;
				J = xSubXi*ySubEta - xSubEta*ySubXi;
				Q = -0.0*exp(-10.0*fabs(j*deltaEta-0.0));
				//-----------------denominators----------------------
				den1 = g11/(deltaXi*deltaXi);
				den2 = g12/(2.0*deltaXi*deltaEta);
				den3 = g22/(deltaEta*deltaEta);
				den4 = J*J*Q;
				den = 2.0*g11/(deltaXi*deltaXi)+2.0*g22/(deltaEta*deltaEta);


				xVal = ( den1*( x->GetAt(i+1, j ) + 
					x->GetAt(i-1, j ) )
					-den2*( x->GetAt(i+1, j+1 ) - 
					x->GetAt(i+1, j-1) 
					- x->GetAt(i-1, j+1) + 
					x->GetAt(i-1, j-1) )
					+ den3*( x->GetAt(i, j+1 ) + 
					x->GetAt(i, j-1 ) ) 
					+ den4*xSubEta )/den;

				x->SetAt( i, j, xVal );

				yVal = ( den1*( y->GetAt(i+1, j ) + 
					y->GetAt(i-1, j ) )
					-den2*( y->GetAt(i+1, j+1 ) - 
					y->GetAt(i+1, j-1) 
					- y->GetAt(i-1, j+1) + 
					y->GetAt(i-1, j-1) )
					+ den3*( y->GetAt(i, j+1 ) + 
					y->GetAt(i, j-1 ) ) 
					+ den4*ySubEta )/den;

				y->SetAt( i, j, yVal );

				//if ( fabs( xVal - x_temp ) > maxX ) maxX = fabs( xVal - x_temp );
				//if ( fabs( yVal - y_temp ) > maxY ) maxY = fabs( xVal - x_temp );
			}
		}
		//cout<<"Limit :"<<limit<<endl;
		//iter++;
	}
	while( eps > 1.0E-7 );
	cout<<"Time :"<<omp_get_wtime()-st<<" "<<iter<<endl;
	return 0;
}
////////////////////////////////////////////////////////////////////////////////
void Poisson::ComputeMetrics( int ix, int jy )
{
	/*	g11 = g12 = g22 = (double)(0);
	xSubXi = xSubEta = ySubXi = ySubEta = (double)(0);

	xSubXi = ( gsl_matrix_get(x, ix+1, jy) - 
	gsl_matrix_get(x, ix-1, jy) )/(2.0*deltaXi);
	xSubEta = ( gsl_matrix_get(x, ix, jy+1) - 
	gsl_matrix_get(x, ix, jy-1) )/(2.0*deltaEta);
	ySubXi = ( gsl_matrix_get(y, ix+1, jy) - 
	gsl_matrix_get(y, ix-1, jy) )/(2.0*deltaXi);
	ySubEta = ( gsl_matrix_get(y, ix, jy+1) - 
	gsl_matrix_get(y, ix, jy-1) )/(2.0*deltaEta);

	g11 = xSubEta*xSubEta + ySubEta*ySubEta;   
	g12 = xSubXi*xSubEta + ySubXi*ySubEta;
	g22 = xSubXi*xSubXi + ySubXi*ySubXi;
	J = xSubXi*ySubEta - xSubEta*ySubXi;
	Q = 0.0;
	//-----------------denominators----------------------
	den1 = g11/(deltaXi*deltaXi);
	den2 = g12/(2.0*deltaXi*deltaEta);
	den3 = g22/(deltaEta*deltaEta);
	den4 = J*J*Q;
	den = 2.0*g11/(deltaXi*deltaXi)+2.0*g22/(deltaEta*deltaEta);*/
}
////////////////////////////////////////////////////////////////////////////////
Poisson::Poisson(int nt=0)
{
	if ( nt == 0 )
	{
		numprocs = omp_get_num_procs( );
		numthreads = numprocs;
	}
	else
		numthreads = nt;
	omp_set_num_threads( numthreads );
}
////////////////////////////////////////////////////////////////////////////////
Poisson::~Poisson(void)
{
}
////////////////////////////////////////////////////////////////////////////////
void Poisson::SetInput( Matrix* xinp, Matrix* yinp )
{
	xpvt = xinp;
	ypvt = yinp;
}
////////////////////////////////////////////////////////////////////////////////
void Poisson::IterateAndSolve(double tol, int max_iter)
{
	int k;   //number of iterations
	k = 0;
	//while ( k< max_iter )
	{
		//#pragma omp parallel
		Solver( );

		//if ( ( maxX <=tol) && ( maxY<=tol ) ) break;
		k = k + 1;
		if ( k%100 == 0 ) cout<<k<<" Max x :"<<maxX<<" Max y :"<<maxY<<endl;
	}
	//cout<<"Number of iterations :"<<k<<endl;
	//cout<<"Eps x :"<<maxX<<" Eps y :"<<maxY<<endl;

}
////////////////////////////////////////////////////////////////////////////////
void Poisson::Add( Matrix* a, Matrix* b, Matrix* sum )
{
	double val;
	int i, j;
	double st = omp_get_wtime();
#pragma omp parallel for shared(sum,a,b) private(i,j, val)
	for ( i=0;i<a->GetNumRows();i++ )
	{
		for ( j=0;j<a->GetNumCols();j++ )
		{
			val = a->GetAt(i,j) + b->GetAt(i,j);
			sum->SetAt( i, j, val );
		}
	}
	cout<<"Time :"<<omp_get_wtime()-st<<endl;
}
////////////////////////////////////////////////////////////////////////////////
void Poisson::SolverSimple( Matrix* inp )
{


}
////////////////////////////////////////////////////////////////////////////////
void Poisson::WriteGrid( string fName )
{
	cout<<"Writing grid to "<<fName<<" ...";
	ofstream fout;

	fout.open( fName.c_str() );
	fout<<"TITLE = "<<'"'<<"Solid Grid"<<'"'<<endl;
	fout<<"VARIABLES = "<<'"'<<'x'<<'"'<<"\t"<<'"'<<'y'<<'"'<<"\t"<<endl;
	fout<<"ZONE "<<"I= "<<xpvt->GetNumRows()<<"\t J= "<<xpvt->GetNumCols()<<endl;
	fout<<"DT = ( SINGLE SINGLE )"<<endl;
	for ( int j=0;j<xpvt->GetNumCols();j++ )
	{
		for ( int i=0;i<xpvt->GetNumRows();i++ )
		{

			fout<<xpvt->GetAt( i, j )<<" "<<
				ypvt->GetAt( i, j )<<endl;
		}
		fout<<endl;      
	}
	cout<<" Done "<<endl;
	cout<<"============================================================="<<endl;
	fout.close();
}
////////////////////////////////////////////////////////////////////////////////