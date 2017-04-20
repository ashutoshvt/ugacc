#include "ccpert.h"
#include "cclinresp.h"

using namespace std;

namespace psi { namespace ugacc {


CCLinResp::CCLinResp(shared_ptr<CCPert> X, shared_ptr<CCPert> Y)
{

   H_ = X->H_;
   CC_ = X->CC_;
   HBAR_ = X->HBAR_;
   CCLambda_ = X->CCLambda_;
   
   no_ = CC_->no_;
   nv_ = CC_->nv_;
   
   X_ = X;
   Y_ = Y;
   
   X1_x_ = X_->X1_ ;
   Y1_x_ = X_->Y1_ ;
   X2_x_ = X_->X2_ ;
   Y2_x_ = X_->Y2_ ;
   
   
   X1_y_ = Y_->X1_ ;
   Y1_y_ = Y_->Y1_ ;
   X2_y_ = Y_->X2_ ;
   Y2_y_ = Y_->Y2_ ;

   pert_x_ = X_->pert_ ;
   pert_y_ = Y_->pert_ ;

// Similarity transformed perturbation operator components
   Aov_x_ =   X_->Aov_;
   Aoo_x_ =   X_->Aoo_;
   Avv_x_ =   X->Avv_;
   Avo_x_ =   X->Avo_;
   Aovoo_x_ = X->Aovoo_;
   Avvvo_x_ = X->Avvvo_;
   Avvoo_x_ = X->Avvoo_;

   Aov_y_ =   Y_->Aov_;
   Aoo_y_ =   Y_->Aoo_;
   Avv_y_ =   Y->Avv_;
   Avo_y_ =   Y->Avo_;
   Aovoo_y_ = Y->Aovoo_;
   Avvvo_y_ = Y->Avvvo_;
   Avvoo_y_ = Y->Avvoo_;


}

//CCLinResp::CCLinResp(shared_ptr<CCPert> X)
//{
//
//H_ = X->H_;
//CC_ = X->CC_;
//HBAR_ = X->HBAR_;
//CCLambda_ = X->CCLambda_;
//
//no_ = CC_->no_;
//nv_ = CC_->nv_;
//
//X_ = X;
//
//}




CCLinResp::~CCLinResp()
{

}

double  CCLinResp::linresp()
{

 int no = no_ ;
 int nv = nv_ ;

 double **t1 = CC_->t1_;
 double ****t2 = CC_->t2_;

 double **l1 =  CCLambda_->l1_;
 double ****l2 = CCLambda_->l2_;

 double **X1_x = X1_x_ ;
 double **Y1_x = Y1_x_ ;

 double ****X2_x = X2_x_ ;
 double ****Y2_x = Y2_x_ ;

 double **X1_y = X1_y_ ;
 double **Y1_y = Y1_y_ ;

 double ****X2_y = X2_y_ ;
 double ****Y2_y = Y2_y_ ;

 double **pert_x = pert_x_;
 double **pert_y = pert_y_;

 double **Aov_x = Aov_x_;
 double **Aoo_x = Aoo_x_;
 double **Avv_x = Avv_x_;
 double **Avo_x = Avo_x_;

 double ****Aovoo_x = Aovoo_x_;
 double ****Avvvo_x = Avvvo_x_;
 double ****Avvoo_x = Avvoo_x_;

 double **Aov_y = Aov_y_;
 double **Aoo_y = Aoo_y_;
 double **Avv_y = Avv_y_;
 double **Avo_y = Avo_y_;

 double ****Aovoo_y = Aovoo_y_;
 double ****Avvvo_y = Avvvo_y_;
 double ****Avvoo_y = Avvoo_y_;



double polar1=0, polar2=0, polar=0;
for (int i=0; i< no; i++)
  for (int a=0; a< nv; a++){
    polar1 += Avo_x[a][i] * Y1_y[i][a] ;
        for (int j=0; j< no; j++) 
          for (int b=0; b< nv; b++){
	    polar1 += (0.50) * (Avvoo_x[a][b][i][j] + Avvoo_x[b][a][j][i])* Y2_y[i][j][a][b] ;
    } 
}

for (int i=0; i< no; i++)
  for (int a=0; a< nv; a++){
        for (int j=0; j< no; j++) 
          for (int b=0; b< nv; b++){
		polar2 +=  l1[i][a] * pert_x[j][b+no] * (2.0 * X2_y[i][j][a][b] - X2_y[i][j][b][a]);  // 2. checked
     }
          for (int c=0; c< nv; c++)
		polar2 +=  l1[i][a] * Avv_x[a][c] * X1_y[i][c];   // 2. checked
          for (int k=0; k< no; k++)
                polar2 -=  l1[i][a] * Aoo_x[k][i] * X1_y[k][a];   // 2. checked


          for (int j=0; j< no; j++)
              for (int b=0; b< nv; b++){         
                 for (int c=0; c< nv; c++)
                       polar2 += l2[i][j][b][c] * X1_y[i][a] * Avvvo_x[b][c][a][j] ;   // 3. checked
	
          for (int k=0; k< no; k++){
	
                polar2 -= 0.5 * l2[i][j][a][b] * X1_y[k][a] * Aovoo_x[k][b][i][j] ;    // 4. checked
                polar2 -= 0.5 * l2[i][j][a][b] * X1_y[k][b] * Aovoo_x[k][a][j][i] ;    // 4. checked
               
		polar2 -= 0.5 * l2[i][j][a][b] * (X2_y[k][j][a][b]) * Aoo_x[k][i] ; // 4. checked
		polar2 -= 0.5 * l2[i][j][a][b] * (X2_y[k][i][b][a]) * Aoo_x[k][j] ; // 4. checked 
        }

           for (int c=0; c< nv; c++){
                polar2 += 0.5 * l2[i][j][a][b] * (X2_y[i][j][a][c]) * Avv_x[b][c] ;  // 4. checked
                polar2 += 0.5 * l2[i][j][a][b] * (X2_y[j][i][b][c]) * Avv_x[a][c] ;  // 4. checked
        }

  }
		polar2 += 2.0 * pert_x[i][a+no] * X1_y[i][a]; // 1. checked
}


  outfile->Printf("\n polarizability first term: %20.14lf\n", polar1);
  outfile->Printf("\n polarizability second term: %20.14lf\n", polar2);
  //outfile->Printf("\n polarizability: %20.14lf\n", polar1 + polar2);
  return -1.0 * (polar1 + polar2) ;


 }

}} // psi::ugaccc

