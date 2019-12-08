#ifndef FORWARD_FUNC_H
#define FORWARD_FUNC_H
#include "ShengShen_head.h"
#include "Parameter_func.h"
#include "Model_func.h"
#include "RickerWavelet_func.h"
#include "IO_func.h"

//double forward(int NX, int NZ, int Nt, int PML, float dx, float dz, float dt);
double forward_func(PARAMETER* param);

#endif
