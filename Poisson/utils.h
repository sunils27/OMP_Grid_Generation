Matrix* ComputeFunction( Matrix* a, Matrix* xvals )
{
	//computes a function of the form y = a0 + a1*x + a2*x^2 + ... + a(n-1)*x^(n-1)
	//and returns the values for a vector of xvals, at each station
	Matrix* ret = new Matrix( xvals->GetNumRows(), 1 ); //there should be only one col
	double ytmp, xtmp, ctmp;
	//loop over each of xvals
	for ( int j=0;j<xvals->GetNumRows();j++ )
	{
		ytmp = 0.0;
		xtmp = xvals->GetAt(j,0);
		for ( int i=0;i<a->GetNumRows();i++ )
		{
			ctmp = a->GetAt(i,0); //coefficient
			ytmp = ytmp + ctmp*pow( xtmp,(double)(a->GetNumRows()-(i+1)) ); //y = sigma( ai * x^i )
		}
		ret->SetAt( j, 0, ytmp );
	}
	return ret;
}

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

void PrintVector( gsl_vector* vect )
{
	for ( size_t i=0;i<vect->size;i++ )
		cout<<i<<" : "<<gsl_vector_get( vect, i )<<endl;
}