#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ugacc {

void N3_ijk(double ***N3, int i, int j, int k, double ****t2, double **t1, double **fock, double ****ints)
{
  int no = moinfo.no;
  int nv = moinfo.nv;

  for(int a=0; a < nv; a++)
    for(int b=0; b < nv; b++)
      for(int c=0; c < nv; c++) {
        double value = 0.0;

        value = ints[i][j][a+no][b+no] * t1[k][c]
              + ints[i][k][a+no][c+no] * t1[j][b]
              + ints[j][k][b+no][c+no] * t1[i][a]
              + t2[i][j][a][b] * fock[k][c+no]
              + t2[i][k][a][c] * fock[j][b+no]
              + t2[j][k][b][c] * fock[i][a+no];

        double denom = fock[i][i] + fock[j][j] + fock[k][k];
        denom -= fock[a+no][a+no] + fock[b+no][b+no] + fock[c+no][c+no];

        N3[a][b][c] = value/denom;
      } // abc

  return;
}

}} // namespace psi::ugacc
