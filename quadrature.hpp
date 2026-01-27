#include <iostream>
#include <cstdlib>
#include <cmath>

#ifndef QUADRATURE
#define QUADRATURE

class quadrature
{
public:
  // Case of Gauss quadrature on \RR with respec to weight exp(-x^2)
  int _nG; // Number of quadrature_points
  double * _X; // Integration points for Gauss Lobatto in [-1,1]
  double * _W; // Weights for Gauss Lobatto
  

public:
   quadrature(int nG=50); // Build the quadrature (50 points --> 2x49 +1)
  ~quadrature(); //Destructor
  //double integrate_1d(double,double, double (*)(double));

};
#endif
