#ifndef RIGHT_FUNC_H
#define RIGHT_FUNC_H

#include "ShengShen_head.h"

//计算Vx,Vz
double right01_func(float **P_now,float **P_pre,float **P_aft,			     
				 float **Vx_now,float **Vx_pre,float **Vx_aft,
				 float **Vz_now,float **Vz_pre,float **Vz_aft,
				 float **DEN,float **absorbx,float **absorbz,struct PARAMETER *param);
				 
//计算K
double right02_func(float **Vx_now,float **Vx_pre,float **Vx_aft,
				 float **K_now,float **K_pre,float **K_aft,
				 float **DEN,struct PARAMETER *param);	
				 
//计算Psa				 
double right03_func(float **K_now,float **K_pre,float **K_aft,
				 float **Psa_now,float **Psa_pre,float **Psa_aft,
				 float **absorbx,float **absorbz,struct PARAMETER *param);
				 
//计算Cita
double right04_func(float **Psa_now,float **Psa_pre,float **Psa_aft,
				 float **Cita_now,float **Cita_pre,float **Cita_aft,
				 float **absorbx,float **absorbz,struct PARAMETER *param);	
				 
//计算P
double right05_func(float **P_now,float **P_pre,float **P_aft,
			     float **Px_now,float **Px_pre,float **Px_aft,
				 float **Pz_now,float **Pz_pre,float **Pz_aft,
				 float **Vz_now,float **Vz_pre,float **Vz_aft,
				 float **K_now,float **K_pre,float **K_aft,
				 float **Psa_now,float **Psa_pre,float **Psa_aft,
				 float **Cita_now,float **Cita_pre,float **Cita_aft,
				 float **V,float **Vv,float **DEN,
				 float **absorbx,float **absorbz,struct PARAMETER *param,
				 float wavelet);	 
#endif