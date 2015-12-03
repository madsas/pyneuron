#include <stdio.h>
#include "hocdec.h"
#define IMPORT extern __declspec(dllimport)
IMPORT int nrnmpi_myid, nrn_nobanner_;

extern void _capump_reg();
extern void _spike_rat_reg();

modl_reg(){
	//nrn_mswindll_stdio(stdin, stdout, stderr);
    if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
	fprintf(stderr, "Additional mechanisms from files\n");

fprintf(stderr," capump.mod");
fprintf(stderr," spike_rat.mod");
fprintf(stderr, "\n");
    }
_capump_reg();
_spike_rat_reg();
}
