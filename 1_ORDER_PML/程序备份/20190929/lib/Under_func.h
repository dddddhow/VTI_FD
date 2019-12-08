#ifndef UNDER_FUNC_H
#define UNDER_FUNC_H

#include "ShengShen_head.h"

double under_func(float **P_now,float **P_pre,float **P_aft,
			     float **Px_now,float **Px_pre,float **Px_aft,
				 float **Pz_now,float **Pz_pre,float **Pz_aft,
				 float **Vx_now,float **Vx_pre,float **Vx_aft,
				 float **Vz_now,float **Vz_pre,float **Vz_aft,
				 float **K_now,float **K_pre,float **K_aft,
				 float **Psa_now,float **Psa_pre,float **Psa_aft,
				 float **Cita_now,float **Cita_pre,float **Cita_aft,
				 float **V,float **Vv,float **DEN,
				 float **absorbx,float **absorbz,struct PARAMETER *param,
				 float wavelet);

double under01_func(float **P_now,float **P_pre,float **P_aft,
			     float **Px_now,float **Px_pre,float **Px_aft,
				 float **Pz_now,float **Pz_pre,float **Pz_aft,
				 float **Vx_now,float **Vx_pre,float **Vx_aft,
				 float **Vz_now,float **Vz_pre,float **Vz_aft,
				 float **K_now,float **K_pre,float **K_aft,
				 float **Psa_now,float **Psa_pre,float **Psa_aft,
				 float **Cita_now,float **Cita_pre,float **Cita_aft,
				 float **V,float **Vv,float **DEN,
				 float **absorbx,float **absorbz,struct PARAMETER *param,
				 float wavelet);
				 
double under02_func(float **P_now,float **P_pre,float **P_aft,
			     float **Px_now,float **Px_pre,float **Px_aft,
				 float **Pz_now,float **Pz_pre,float **Pz_aft,
				 float **Vx_now,float **Vx_pre,float **Vx_aft,
				 float **Vz_now,float **Vz_pre,float **Vz_aft,
				 float **K_now,float **K_pre,float **K_aft,
				 float **Psa_now,float **Psa_pre,float **Psa_aft,
				 float **Cita_now,float **Cita_pre,float **Cita_aft,
				 float **V,float **Vv,float **DEN,
				 float **absorbx,float **absorbz,struct PARAMETER *param,
				 float wavelet);				 
double under03_func(float **P_now,float **P_pre,float **P_aft,
			     float **Px_now,float **Px_pre,float **Px_aft,
				 float **Pz_now,float **Pz_pre,float **Pz_aft,
				 float **Vx_now,float **Vx_pre,float **Vx_aft,
				 float **Vz_now,float **Vz_pre,float **Vz_aft,
				 float **K_now,float **K_pre,float **K_aft,
				 float **Psa_now,float **Psa_pre,float **Psa_aft,
				 float **Cita_now,float **Cita_pre,float **Cita_aft,
				 float **V,float **Vv,float **DEN,
				 float **absorbx,float **absorbz,struct PARAMETER *param,
				 float wavelet);				 
double under04_func(float **P_now,float **P_pre,float **P_aft,
			     float **Px_now,float **Px_pre,float **Px_aft,
				 float **Pz_now,float **Pz_pre,float **Pz_aft,
				 float **Vx_now,float **Vx_pre,float **Vx_aft,
				 float **Vz_now,float **Vz_pre,float **Vz_aft,
				 float **K_now,float **K_pre,float **K_aft,
				 float **Psa_now,float **Psa_pre,float **Psa_aft,
				 float **Cita_now,float **Cita_pre,float **Cita_aft,
				 float **V,float **Vv,float **DEN,
				 float **absorbx,float **absorbz,struct PARAMETER *param,
				 float wavelet);
double under05_func(float **P_now,float **P_pre,float **P_aft,
			     float **Px_now,float **Px_pre,float **Px_aft,
				 float **Pz_now,float **Pz_pre,float **Pz_aft,
				 float **Vx_now,float **Vx_pre,float **Vx_aft,
				 float **Vz_now,float **Vz_pre,float **Vz_aft,
				 float **K_now,float **K_pre,float **K_aft,
				 float **Psa_now,float **Psa_pre,float **Psa_aft,
				 float **Cita_now,float **Cita_pre,float **Cita_aft,
				 float **V,float **Vv,float **DEN,
				 float **absorbx,float **absorbz,struct PARAMETER *param,
				 float wavelet);				 
#endif