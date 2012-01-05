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

void WriteGrid( Matrix* xgrid, Matrix* ygrid, Matrix* u, Matrix* T, string fName )
{
	cout<<"Writing grid to "<<fName<<" ...";
	ofstream fout;

	fout.open( fName.c_str() );
	fout<<"TITLE = "<<'"'<<"Solid Grid"<<'"'<<endl;
	fout<<"VARIABLES = "<<'"'<<'x'<<'"'<<"\t"<<'"'<<'y'<<'"'<<"\t";
	if ( u != NULL )
		fout<<'"'<<'u'<<'"'<<"\t";
	if ( T != NULL )
		fout<<'"'<<'T'<<'"'<<"\t";
	fout<<endl;
	fout<<"ZONE "<<"I= "<<xgrid->GetNumCols()<<"\t J= "<<xgrid->GetNumRows()<<endl;
	fout<<"DT = ( SINGLE SINGLE ";
	if ( u != NULL )
		fout<<"SINGLE ";
	if ( T != NULL )
		fout<<" SINGLE ";
	fout<<" )"<<endl;
	for ( int j=0;j<xgrid->GetNumRows();j++ )
	{
		for ( int i=0;i<xgrid->GetNumCols();i++ )
		{

			fout<<xgrid->GetAt( j, i )<<" "<<
				ygrid->GetAt( j, i )<<" ";
			if ( u != NULL )
				fout<<u->GetAt( j, i ) <<" ";
			if ( T != NULL )
				fout<<T->GetAt( j, i );
			fout<<endl;
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

void ApplyFieldDirichlet( Matrix* field, int i, int j )
{
	//apply Dirichlet BC on a field Matrix

}

void UpdateNeumannBoundary( Matrix& f, int alongi, bool iupdate, int alongj, bool jupdate )
{
	double tmp;
	if ( iupdate ) //update for a given i, along column
	{
		if ( alongi != 0 )
		{
			for ( int j=0;j<f.GetNumCols();j++ )
			{
				tmp = ( 4.0*f.GetAt(alongi-1,j) - f.GetAt(alongi-2,j) )/3.0;
				f.SetAt( alongi, j, tmp );
			}
		}
		else if ( alongi == 0 )
		{
			for ( int j=0;j<f.GetNumCols();j++ )
			{
				tmp = ( 4.0*f.GetAt(alongi+1,j) - f.GetAt(alongi+2,j) )/3.0;
				f.SetAt( alongi, j, tmp );
			}
		}
	}
	if ( jupdate ) //update for a given j, along row
	{
		if ( alongj != 0 )
		{
			for ( int i=0;i<f.GetNumRows();i++ )
			{
				tmp = ( 4.0*f.GetAt(i,alongj-1) - f.GetAt(i,alongj-2) )/3.0;
				f.SetAt( i, alongj, tmp );
			}
		}
		else if ( alongj == 0 )
		{
			for ( int i=0;i<f.GetNumRows();i++ )
			{
				tmp = ( 4.0*f.GetAt(i,alongj+1) - f.GetAt(i,alongj+2) )/3.0;
				f.SetAt( i, alongj, tmp );
			}
		}
	}
}

void testfunc( )
{
	cout<<"test func "<<endl;
}

typedef struct
{
	int a;
} params, *ptr_params;
void testfunc2( int a )
{
	//ptr_params ptr = (ptr_params)(params);
	//int a = ptr->a ;
	cout<<"test function 2 "<<a<<endl;
}

void readsyncfunc( Matrix* ds_from, Matrix* ds_to, omp_lock_t& lck, int synx, int syny )
{
	//acquire lock, copy from to to and release lock
	//lock ds_from, read ds_from, unlock ds_from, write to f (ds_to)
	omp_set_lock( &lck );
	{
		for ( int j=0;j<ds_to->GetNumCols();j++ )
		{
			ds_to->SetAt( synx, j, ds_from->GetAt(synx,j) );
		}
	}
	omp_unset_lock( &lck );
}

void writesyncfunc( Matrix* ds_from, Matrix* ds_to, omp_lock_t&  lck, int synx, int syny )
{
	//read ds_from, lock ds_to, write to ds_to, unlock ds_to
	omp_set_lock( &lck );
	{
		for ( int j=0;j<ds_to->GetNumCols();j++ )
		{
			ds_to->SetAt( synx, j, ds_from->GetAt(synx,j) );
		}
	}
	omp_unset_lock( &lck );
}