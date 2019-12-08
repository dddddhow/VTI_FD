#ifndef IO_FUNC_H
#define IO_FUNC_H
#include"ShengShen_head.h"

//double IO_func(int NX, int NZ, int PML, int Nt, float dx, float dz, float dt, int nx_location, int nz_location, float **record);

double IO_func(struct PARAMETER *param, float **record);



#endif

