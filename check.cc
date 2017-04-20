#include "ccpert.h"
#include "cclinresp.h"
#include "hbar.h"
#include "cclambda.h"
#include <psi4/libciomr/libciomr.h>
#include <psi4/libqt/qt.h>
#include "array.h"

using namespace std;

namespace psi { namespace ugacc {
void  CCLinResp::check_linear(shared_ptr<CCPert> X)
{


 double **l1_ = CCLambda_->l1_;
 double ****l2_ = CCLambda_->l2_; 

// g, m,n,o,p
  int no = no_;
  int nv = nv_;
  double **Hvv = HBAR_->Hvv_;
  double **Hoo = HBAR_->Hoo_;
  double **Hov = HBAR_->Hov_;
  double ****Hvvvo = HBAR_->Hvvvo_;
  double ****Hovoo = HBAR_->Hovoo_;
  double ****Hovvo = HBAR_->Hovvo_;
  double ****Hovov = HBAR_->Hovov_;
  double ****Hvvvv = HBAR_->Hvvvv_;
  double ****Hoooo = HBAR_->Hoooo_;
  double ****Hvovv = HBAR_->Hvovv_;
  double ****Hooov = HBAR_->Hooov_;
  double ****ints =  H_->ints_;

 double ****R2 = init_4d_array(no, no, nv, nv); // residual
 double ****R2_ = init_4d_array(no, no, nv, nv); // residual
 double **R1 = block_matrix(no, nv); // residual

 double **Aov_ = X->Aov_;
 double **Aoo_ = X->Aoo_;
 double **Avo_ = X->Avo_;
 double **Avv_ = X->Avv_;
 double ****Aovoo_ = X->Aovoo_;
 double ****Avvvo_ = X->Avvvo_;
 double ****Avvoo_ = X->Avvoo_;
 double **pert_ = X->pert_;

 double **X1_ = X->X1_;
 double ****X2_ = X->X2_;


     for(int i=0; i < no; i++)
        for(int a=0; a < nv; a++){

         R1[i][a] = 2*pert_[i][a+no] ;         // (g)

      for(int m=0; m < no; m++)
          R1[i][a] -= l1_[m][a] * Aoo_[i][m];  // (m)

      for(int e=0; e < nv; e++)
          R1[i][a] += l1_[i][e] * Avv_[e][a];  // (m)

      for(int m=0; m < no; m++)
        for(int e=0; e < nv; e++)
          for(int f=0; f < nv; f++)
            R1[i][a] += l2_[i][m][f][e] * Avvvo_[f][e][a][m];     // (n)

      for(int m=0; m < no; m++)
        for(int n=0; n < no; n++)
          for(int e=0; e < nv; e++)
          {
            R1[i][a] -= 0.5 * l2_[m][n][e][a] * Aovoo_[i][e][n][m];    //  (n)
            R1[i][a] -= 0.5 * l2_[m][n][a][e] * Aovoo_[i][e][m][n];    //  (n)
          }

       }

     for(int i=0; i < no; i++)
         for(int j=0; j < no; j++)
           for(int a=0; a < nv; a++)
             for(int b=0; b < nv; b++) {

          R2[i][j][a][b] += 2*l1_[j][b] * pert_[i][a+no] - l1_[i][b] * pert_[j][a+no];  // (o)

          for(int e=0; e < nv; e++)
            R2[i][j][a][b] += l2_[i][j][e][b]*Avv_[e][a];  // (p)

          for(int m=0; m < no; m++)
            R2[i][j][a][b] -= l2_[m][j][a][b]*Aoo_[i][m];  // (p)
        }

        for(int i=0; i < no; i++)
         for(int j=0; j < no; j++)
           for(int a=0; a < nv; a++)
             for(int b=0; b < nv; b++) {
                R2_[i][j][a][b] = R2[i][j][a][b] + R2[j][i][b][a] ;
               }

      double value1=0, value2=0;

for(int i=0; i < no; i++)
           for(int a=0; a < nv; a++){
               value1 += R1[i][a] * X1_[i][a] ;
               for(int j=0; j < no; j++)
                  for(int b=0; b < nv; b++) {
                      value2 += R2_[i][j][a][b] * 0.5 * X2_[i][j][a][b] ;
        }
    }
        outfile->Printf("\n LCX1: %20.14lf \n", value1);
        outfile->Printf("\n LCX2: %20.14lf \n", value2);
        outfile->Printf("\n LCX: %20.14lf \n", value1 + value2);




}



void  CCLinResp::check_quadratic(shared_ptr<CCPert> X, shared_ptr<CCPert> Y)
{
  int no = no_;
  int nv = nv_;
  double **Hvv = HBAR_->Hvv_;
  double **Hoo = HBAR_->Hoo_;
  double **Hov = HBAR_->Hov_;
  double ****Hvvvo = HBAR_->Hvvvo_;
  double ****Hovoo = HBAR_->Hovoo_;
  double ****Hovvo = HBAR_->Hovvo_;
  double ****Hovov = HBAR_->Hovov_;
  double ****Hvvvv = HBAR_->Hvvvv_;
  double ****Hoooo = HBAR_->Hoooo_;
  double ****Hvovv = HBAR_->Hvovv_;
  double ****Hooov = HBAR_->Hooov_;
  double ****ints =  H_->ints_;

  double **l1_ = CCLambda_->l1_;
  double ****l2_ = CCLambda_->l2_;

 double ****R2 = init_4d_array(no, no, nv, nv); // residual
 double ****R2_ = init_4d_array(no, no, nv, nv); // residual
 double **R1 = block_matrix(no, nv); // residual


 double **X1_ = X->X1_;
 double ****X2_ = X->X2_;

 double **X1_y = Y->X1_;
 double ****X2_y = Y->X2_;

   for(int i=0; i < no; i++)
     for(int a=0; a < nv; a++){

   for(int m=0; m < no; m++)
     for(int e=0; e < nv; e++)
       R1[i][a] += 2.0 * X1_[m][e] * H_->L_[i][m][a+no][e+no];   // (i)


       for(int m=0; m < no; m++)
        for(int e=0; e < nv; e++){
          R1[i][a] -= X1_[m][e] * Hov[m][a] * l1_[i][e];  // (q)       
          R1[i][a] -= X1_[m][e] * Hov[i][e] * l1_[m][a];  // (q)       
          for(int n=0; n < no; n++){
            R1[i][a] -= X1_[m][e] * (2*Hooov[m][i][n][a] - Hooov[i][m][n][a]) * l1_[n][e];  // (q)  
            R1[i][a] -= X1_[m][e] * (2*Hooov[i][m][n][e] - Hooov[m][i][n][e]) * l1_[n][a];  // (q) 
        }
          for(int f=0; f < nv; f++){
            R1[i][a] += X1_[m][e] * (2*Hvovv[f][m][a][e] - Hvovv[f][m][e][a]) * l1_[i][f];   // (q)
            R1[i][a] += X1_[m][e] * (2*Hvovv[f][i][e][a] - Hvovv[f][i][a][e]) * l1_[m][f];   // (q) 
        }
   }
    for(int m=0; m < no; m++)
      for(int n=0; n < no; n++)
        for(int e=0; e < nv; e++)
          for(int f=0; f < nv; f++){
            R1[i][a] += 2.0 * X2_[m][n][e][f] * H_->L_[i][m][a+no][e+no] * l1_[n][f];   // (r)
            R1[i][a] -= X2_[m][n][e][f] * H_->L_[i][m][a+no][f+no] * l1_[n][e];         // (r)
            R1[i][a] -= X2_[m][n][e][f] * H_->L_[m][i][e+no][f+no] * l1_[n][a];         // (r)
            R1[i][a] -= X2_[m][n][e][f] * H_->L_[n][m][f+no][a+no] * l1_[i][e];         // (r)
}

      for(int m=0; m < no; m++)
        for(int e=0; e < nv; e++){
          for(int n=0; n < no; n++)
            for(int f=0; f < nv; f++){
              R1[i][a] -= X1_[m][e] * Hovov[m][f][n][a] * l2_[n][i][e][f];     // (s)
              R1[i][a] -= X1_[m][e] * Hovov[i][f][n][e] * l2_[n][m][a][f];     // (s)
              R1[i][a] -= X1_[m][e] * Hovvo[m][f][a][n] * l2_[n][i][f][e];     // (s)  
              R1[i][a] -= X1_[m][e] * Hovvo[i][f][e][n] * l2_[n][m][f][a];     // (s)
                }
            for(int f=0; f < nv; f++)
              for(int g=0; g < nv; g++)
                R1[i][a] += X1_[m][e] * Hvvvv[f][g][a][e] * l2_[i][m][f][g];   // (s)
            for(int n=0; n < no; n++)
              for(int o=0; o < no; o++)
                R1[i][a] += X1_[m][e] * Hoooo[i][m][n][o] * l2_[n][o][a][e];   // (s)
         }

  //double **Zvv = block_matrix(nv,nv);
  //double ****t2 = CC_->t2_;

   ////    for(int f=0; f < nv; f++)
   ////       for(int b=0; b < nv; b++)
   ////          for(int j=0; j < no; j++)
   ////              for(int m=0; m < no; m++)
   ////                  for(int e=0; e < nv; e++)
   ////                    Zvv[f][b] += t2[j][m][f][e] * l2_[j][m][b][e];

   ////     double **Zoo = block_matrix(no,no);

   ////    for(int n=0; n < no; n++)
   ////       for(int j=0; j < no; j++)
   ////          for(int b=0; b < nv; b++)
   ////            for(int f=0; f < nv; f++) {
   ////                       Zoo[n][i] += t2[j][n][b][f] * l2_[j][i][b][f];
   ////     }

   ////    for(int n=0; n < no; n++)
   ////       for(int b=0; b < nv; b++)
   ////          for(int f=0; f < nv; f++){
   ////           R1[i][a] -= 0.5 * H_->L_[i][n][a+no][f+no] * Zvv[f][b] * X1_[n][b];
   ////           R1[i][a] -= 0.5 * H_->L_[n][i][b+no][f+no] * Zvv[f][a] * X1_[n][b];
   ////     }

   ////     for(int m=0; m < no; m++)
   ////        for(int n=0; n < no; n++)
   ////           for(int e=0; e < nv; e++){
   ////              R1[i][a] -= 0.5 * H_->L_[m][n][e+no][a+no] * Zoo[n][i] * X1_[m][e];
   ////              R1[i][a] -= 0.5 * H_->L_[i][n][a+no][e+no] * Zoo[n][m] * X1_[m][e];
//     }

//      for(int m=0; m < no; m++)
//        for(int n=0; n < no; n++)
//          for(int e=0; e < nv; e++)
//            for(int f=0; f < nv; f++)
//              for(int j=0; j < no; j++)
//                 for(int b=0; b < nv; b++){
//                    R1[i][a] -=  2.0 * H_->L_[i][n][a+no][f+no] * X1_[n][b] * t2[j][m][f][e] * l2_[j][m][b][e] ;
//                    R1[i][a] -=  2.0 * H_->L_[i][n][a+no][f+no] * X1_[j][f] * t2[m][n][e][b] * l2_[m][j][e][b] ;
//
//                    //R1[i][a] -=  H_->L_[i][n][a+no][f+no] * X1_[n][b] * t2[j][m][f][e] * l2_[j][m][b][e] ;
//                    //R1[i][a] -=  H_->L_[m][i][e+no][f+no] * X1_[m][e] * t2[j][n][f][b] * l2_[j][n][a][b] ;
//
//                    //R1[i][a] -=   H_->L_[j][n][f+no][a+no] * X1_[m][e] * t2[j][n][f][b] * l2_[j][i][f][b] ;
//                    //R1[i][a] -=   H_->L_[m][n][e+no][f+no] * X1_[j][f] * t2[m][n][e][b] * l2_[m][j][e][b] ;
//        }


      for(int m=0; m < no; m++)
        for(int n=0; n < no; n++)
          for(int e=0; e < nv; e++)
            for(int f=0; f < nv; f++){
              R1[i][a] -=  X2_[m][n][e][f] * Hov[m][a] * l2_[n][i][f][e];   // (t)
              R1[i][a] -=  X2_[m][n][e][f] * Hov[i][e] * l2_[n][m][f][a];   // (t)
              for(int g=0; g < nv; g++){
                R1[i][a] -=  X2_[m][n][e][f] * Hvovv[g][n][e][a] * l2_[i][m][f][g] ;   //(t)
                R1[i][a] -=  X2_[m][n][e][f] * Hvovv[g][n][a][e] * l2_[m][i][f][g] ;   //(t)
                R1[i][a] -=  X2_[m][n][e][f] * Hvovv[g][i][e][f] * l2_[n][m][a][g] ;      //(t) 
                R1[i][a] +=  X2_[m][n][e][f] * (2*Hvovv[g][m][a][e] - Hvovv[g][m][e][a]) * l2_[n][i][f][g] ;   //(t)
                R1[i][a] +=  X2_[m][n][e][f] * (2*Hvovv[g][i][e][a] - Hvovv[g][i][a][e]) * l2_[n][m][f][g] ;      //(t)
              }

        for(int o=0; o < no; o++){
                R1[i][a] +=  X2_[m][n][e][f] * Hooov[m][n][o][a] * l2_[o][i][e][f] ; // (t)    
                R1[i][a] +=  X2_[m][n][e][f] * Hooov[i][n][o][e] * l2_[o][m][a][f] ; // (t)    
                R1[i][a] +=  X2_[m][n][e][f] * Hooov[m][i][o][f] * l2_[o][n][e][a] ; // (t)    
                R1[i][a] -=  X2_[m][n][e][f] * (2*Hooov[m][i][o][a] - Hooov[i][m][o][a]) * l2_[n][o][f][e] ;   // (t)
                R1[i][a] -=  X2_[m][n][e][f] * (2*Hooov[i][m][o][e] - Hooov[m][i][o][e]) * l2_[n][o][f][a] ;      // (t)
              }
       }
}
//
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++) {

        for(int m=0; m < no; m++)
            for(int e=0; e < nv; e++) {
              R2[i][j][a][b] -=  X1_[m][e] * l1_[j][a] * H_->L_[m][i][e+no][b+no]; // (u)
              R2[i][j][a][b] -=  X1_[m][e] * l1_[m][b] * H_->L_[i][j][a+no][e+no]; // (u)
              R2[i][j][a][b] -=  X1_[m][e] * l1_[i][e] * H_->L_[j][m][b+no][a+no]; // (u)
              R2[i][j][a][b] += 2.0 * X1_[m][e] * l1_[j][b] * H_->L_[i][m][a+no][e+no]; // (u) 
         }

          for(int m=0; m < no; m++)
            for(int e=0; e < nv; e++) {
              R2[i][j][a][b] -=  Hov[m][a] * l2_[j][i][b][e] * X1_[m][e]; // (w)
              R2[i][j][a][b] -=  Hov[i][e] * l2_[j][m][b][a] * X1_[m][e]; // (w)
              for(int f=0; f < nv; f++){
                R2[i][j][a][b] -=  Hvovv[f][m][b][a] * l2_[j][i][f][e] * X1_[m][e];   // (w)
                R2[i][j][a][b] -=  Hvovv[f][j][e][a] * l2_[m][i][f][b] * X1_[m][e];   // (w)
                R2[i][j][a][b] -=  Hvovv[f][i][b][e] * l2_[j][m][f][a] * X1_[m][e];   // (w)
                R2[i][j][a][b] +=  (2 * Hvovv[f][m][a][e] - Hvovv[f][m][e][a]) * l2_[j][i][b][f] * X1_[m][e];  // (w)
                R2[i][j][a][b] +=  (2 * Hvovv[f][i][e][a] - Hvovv[f][i][a][e]) * l2_[j][m][b][f] * X1_[m][e];  // (w)
              }
              for(int n=0; n < no; n++){
                 R2[i][j][a][b] += Hooov[j][m][n][a] * l2_[n][i][b][e] * X1_[m][e];   // (w)
                 R2[i][j][a][b] += Hooov[m][j][n][a] * l2_[n][i][e][b] * X1_[m][e];   // (w)
                 R2[i][j][a][b] += Hooov[j][i][n][e] * l2_[n][m][b][a] * X1_[m][e];   // (w)
                 R2[i][j][a][b] -= (2 * Hooov[m][i][n][a] - Hooov[i][m][n][a]) * l2_[j][n][b][e] * X1_[m][e];  // (w)
                 R2[i][j][a][b] -= (2 * Hooov[i][m][n][e] - Hooov[m][i][n][e]) * l2_[j][n][b][a] * X1_[m][e];  // (w)
              }
        }
for(int m=0; m < no; m++)
            for(int n=0; n < no; n++)
              for(int e=0; e < nv; e++)
                for(int f=0; f < nv; f++){

                        /* x terms */
                  R2[i][j][a][b] += 0.25 * ints[n][m][a+no][b+no] * l2_[j][i][e][f] * X2_[m][n][e][f] ;
                  R2[i][j][a][b] += 0.25 * ints[m][n][a+no][b+no] * l2_[j][i][f][e] * X2_[m][n][e][f] ;

                  R2[i][j][a][b] += 0.25 * ints[j][n][a+no][e+no] * l2_[m][i][f][b] * X2_[m][n][e][f] ;
                  R2[i][j][a][b] += 0.25 * ints[j][m][a+no][f+no] * l2_[n][i][e][b] * X2_[m][n][e][f] ;

                  R2[i][j][a][b] += 0.25 * ints[n][j][a+no][e+no] * l2_[m][i][b][f] * X2_[m][n][e][f] ;
                  R2[i][j][a][b] += 0.25 * ints[m][j][a+no][f+no] * l2_[n][i][b][e] * X2_[m][n][e][f] ;

                  R2[i][j][a][b] += 0.25 * ints[j][n][e+no][a+no] * l2_[i][m][f][b] * X2_[m][n][e][f] ;
                  R2[i][j][a][b] += 0.25 * ints[j][m][f+no][a+no] * l2_[i][n][e][b] * X2_[m][n][e][f] ;

                  R2[i][j][a][b] += 0.25 * ints[n][j][e+no][a+no] * l2_[i][m][b][f] * X2_[m][n][e][f] ;
                  R2[i][j][a][b] += 0.25 * ints[m][j][f+no][a+no] * l2_[i][n][b][e] * X2_[m][n][e][f] ;

                  R2[i][j][a][b] += 0.25 * ints[j][i][e+no][f+no] * l2_[n][m][a][b] * X2_[m][n][e][f] ;
                  R2[i][j][a][b] += 0.25 * ints[j][i][f+no][e+no] * l2_[m][n][a][b] * X2_[m][n][e][f] ;

                  R2[i][j][a][b] -=  0.5 * H_->L_[i][m][a+no][f+no] * l2_[j][n][b][e] * X2_[m][n][e][f] ;
                  R2[i][j][a][b] -=  0.5 * H_->L_[i][n][a+no][e+no] * l2_[j][m][b][f] * X2_[m][n][e][f] ;
                  R2[i][j][a][b] -=  0.5 * H_->L_[m][i][e+no][f+no] * l2_[j][n][b][a] * X2_[m][n][e][f] ;
                  R2[i][j][a][b] -=  0.5 * H_->L_[m][n][e+no][a+no] * l2_[j][i][b][f] * X2_[m][n][e][f] ;
                  R2[i][j][a][b] -=  0.5 * H_->L_[n][m][f+no][a+no] * l2_[j][i][b][e] * X2_[m][n][e][f] ;
                  R2[i][j][a][b] -=  0.5 * H_->L_[n][i][f+no][e+no] * l2_[j][m][b][a] * X2_[m][n][e][f] ;

                  R2[i][j][a][b] -=   H_->L_[i][j][a+no][e+no] * l2_[n][m][f][b] * X2_[m][n][e][f] ;
                  R2[i][j][a][b] -=   H_->L_[i][m][a+no][b+no] * l2_[n][j][f][e] * X2_[m][n][e][f] ;
                  R2[i][j][a][b] -=   H_->L_[m][j][e+no][a+no] * l2_[n][i][f][b] * X2_[m][n][e][f] ;

                  R2[i][j][a][b] += 2.0 * H_->L_[i][m][a+no][e+no] * l2_[n][j][f][b] * X2_[m][n][e][f] ;

         }

     }

      for(int i=0; i < no; i++)
         for(int j=0; j < no; j++)
           for(int a=0; a < nv; a++)
             for(int b=0; b < nv; b++) {
                R2_[i][j][a][b] = R2[i][j][a][b] + R2[j][i][b][a] ;
               }

        double value1=0, value2=0;
      for(int i=0; i < no; i++)
           for(int a=0; a < nv; a++){
               value1 += R1[i][a] * X1_y[i][a] ;
               for(int j=0; j < no; j++)
                  for(int b=0; b < nv; b++) {
                      value2 += R2_[i][j][a][b] * 0.5 * X2_y[i][j][a][b] ;
        }
 }

      double **Zvv = block_matrix(nv,nv);

       for(int f=0; f < nv; f++)
          for(int b=0; b < nv; b++) 
             for(int m=0; m < no; m++) 
                      for(int n=0; n < no; n++)
                   for(int e=0; e < nv; e++) {
                     Zvv[f][b] += (X1_[m][e] * X1_y[n][b] + X1_[n][b] * X1_y[m][e]) * H_->L_[m][n][e+no][f+no];
      } 

      double **Zoo = block_matrix(no,no);

       for(int n=0; n < no; n++)
          for(int j=0; j < no; j++)
             for(int m=0; m < no; m++)
                   for(int e=0; e < nv; e++)
                      for(int f=0; f < nv; f++) {
                          Zoo[n][j] += (X1_[m][e] * X1_y[j][f] + X1_[j][f] * X1_y[m][e]) * H_->L_[m][n][e+no][f+no];
        }

 double ****t2 = CC_->t2_;

 for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++)
           for(int f=0; f < nv; f++){
            value1 -= l2_[i][j][a][b] * t2[i][j][a][f] * Zvv[f][b] ;
      }

for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++)
           for(int n=0; n < no; n++){
              value1 -= l2_[i][j][a][b] * t2[i][n][a][b] * Zoo[n][j] ;
        }




        outfile->Printf("\n LHX1Y1 ....: %20.14lf \n", value1);
        outfile->Printf("\n ..LHX1Y2 + LHX2Y1 + LHX2Y2: %20.14lf \n", value2);
        outfile->Printf("\n Total Quadratic term: %20.14lf \n", value1 + value2);

}
}}
