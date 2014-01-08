#include <libciomr/libciomr.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ugacc {

void cleanup(void)
{
  int no = moinfo.no; 
  int nv = moinfo.nv;
  int nact = moinfo.nact;

  free_block(moinfo.fock);
  free_4d_array(moinfo.ints, nact, nact, nact);
  free_4d_array(moinfo.L, nact, nact, nact);
  free_block(moinfo.D1);
  free_4d_array(moinfo.D2, no, no, nv);

  free_block(moinfo.t1);
  free_block(moinfo.t1old);
  free_4d_array(moinfo.t2, no, no, nv);
  free_4d_array(moinfo.t2old, no, no, nv);
  free_4d_array(moinfo.tau, no, no, nv);
  free_4d_array(moinfo.ttau, no, no, nv);

  free_block(moinfo.Fme);
  free_block(moinfo.Fmi);
  free_block(moinfo.Fae);
  free_4d_array(moinfo.Wmnij, no, no, no);
  free_4d_array(moinfo.Wmbej, no, nv, nv);
  free_4d_array(moinfo.Wmbje, no, nv, no);

  free_block(moinfo.Hoo);
  free_block(moinfo.Hvv);
  free_block(moinfo.Hov);
  free_4d_array(moinfo.Hoooo, no, no, no);
  free_4d_array(moinfo.Hvvvv, nv, nv, nv);
  free_4d_array(moinfo.Hovov, no, nv, no);
  free_4d_array(moinfo.Hovvo, no, nv, nv);
  free_4d_array(moinfo.Hvovv, nv, no, nv);
  free_4d_array(moinfo.Hooov, no, no, no);
  free_4d_array(moinfo.Hovoo, no, nv, no);
  free_4d_array(moinfo.Hvvvo, nv, nv, nv);

  free_block(moinfo.l1);
  free_block(moinfo.l1old);
  free_4d_array(moinfo.l2, no, no, nv);
  free_4d_array(moinfo.l2old, no, no, nv);

  if(params.wfn == "CCSD_T") {
    free_block(moinfo.s1);
    free_4d_array(moinfo.s2, no, no, nv);
  }
}

}} // namespace psi::ugacc