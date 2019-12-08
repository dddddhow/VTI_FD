#ifndef MODEL_FUNC_H
#define MODEL_FUNC_H
#include "ShengShen_head.h"

//double model_func(int NX, int NZ, int PML, float **V,float **Vv,float **DEN,float epsilon,float delt,float yita);
double model_func(struct PARAMETER *param, float **V, float **Vv, float **DEN);

#endif
