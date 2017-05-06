#ifndef CCRESP_H
#define CCRESP_H

#include "ccpert.h"

using namespace std;

namespace psi { namespace ugacc {

class CCLinResp {
public:
  CCLinResp(shared_ptr<CCPert> X, shared_ptr<CCPert> Y);
  //CCLinResp(shared_ptr<CCPert> X);
  ~CCLinResp();
  double linresp(shared_ptr<CCPert> X, shared_ptr<CCPert> Y);
  double  check_quadratic(shared_ptr<CCPert> X, shared_ptr<CCPert> Y);
  void check_linear(shared_ptr<CCPert> X);

protected:
  shared_ptr<Hamiltonian> H_;
  shared_ptr<CCWfn> CC_;
  shared_ptr<HBAR> HBAR_;
  shared_ptr<CCLambda> CCLambda_;

  int no_;
  int nv_;
  shared_ptr<CCPert> Y_;
  shared_ptr<CCPert> X_;

  double **   X1_x_;
  double **   Y1_x_;
  double **** X2_x_;
  double **** Y2_x_;

 
  double **   X1_y_;
  double **   Y1_y_;
  double **** X2_y_;
  double **** Y2_y_;

  double ** pert_x_;
  double ** pert_y_;

  // Similarity transformed perturbation operator components
  double **    Aov_x_;
  double **    Aoo_x_;
  double **    Avv_x_;
  double **    Avo_x_;
  double ****Aovoo_x_;
  double ****Avvvo_x_;
  double ****Avvoo_x_;


  double **    Aov_y_;
  double **    Aoo_y_;
  double **    Avv_y_;
  double **    Avo_y_;
  double ****Aovoo_y_;
  double ****Avvvo_y_;
  double ****Avvoo_y_;

};

}} // psi::ugaccc

#endif
