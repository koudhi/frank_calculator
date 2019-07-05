#ifndef _ENERGY_EVAL_H_
#define _ENERGY_EVAL_H_
#include <iostream>
#include "../include/dir.h"

int elastic_energy_eval (int Nx, int Ny, int Nz, vector *n, double *energy, vector lat_size);
int energy_integrator (double *K, int Nx, int Ny, int Nz);
#endif /*energy_eval.h*/
