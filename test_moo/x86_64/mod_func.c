#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;

extern void _kv_reg(void);
extern void _na_reg(void);
extern void _NMDA_reg(void);
extern void _ProbAMPANMDA2_ratio_reg(void);
extern void _ProbAMPANMDA_reg(void);

void modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," \"mechanisms//kv.mod\"");
    fprintf(stderr," \"mechanisms//na.mod\"");
    fprintf(stderr," \"mechanisms//NMDA.mod\"");
    fprintf(stderr," \"mechanisms//ProbAMPANMDA2_ratio.mod\"");
    fprintf(stderr," \"mechanisms//ProbAMPANMDA.mod\"");
    fprintf(stderr, "\n");
  }
  _kv_reg();
  _na_reg();
  _NMDA_reg();
  _ProbAMPANMDA2_ratio_reg();
  _ProbAMPANMDA_reg();
}
