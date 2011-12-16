#ifndef COEFFICIENTS_H
#define COEFFICIENTS_H
//class to compute coefficients of the curve based on the input y-values
#include <cmath>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <iostream>
using namespace std;
#include <stdio.h>
#include "Matrix.h"


class Coefficients
{
   private:
      int polyDegree;
      gsl_vector* polyCoeffs; //values calculated, coefficients fo the polynomial
      gsl_vector* yValues; //values set from outside the y-coordinates

      gsl_vector* xValues; //values of X coordinates
      double xLength;   //length of domain in x-direction
      gsl_matrix* matrix;

   public:
      Coefficients( int degreePoly, gsl_vector* yvals, double xl );
	  Coefficients( int degreePoly, Matrix* yvals, double xl );
      ~Coefficients( );
      int getDegree( );
      gsl_vector* getCoeffs( );
	  Matrix* getCoeffs2( );
      gsl_vector* get_Y_Vals( );
      void CalculateCoeffs( );
      void ShowMatrix( );
};
#endif //COEFFICIENTS_H