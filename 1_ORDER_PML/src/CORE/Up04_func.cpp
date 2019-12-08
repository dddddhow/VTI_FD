#include "../../lib/Up_func.h"

double up04_func(float **P_now,float **P_pre,float **P_aft,
        float **Px_now,float **Px_pre,float **Px_aft,
        float **Pz_now,float **Pz_pre,float **Pz_aft,
        float **Vx_now,float **Vx_pre,float **Vx_aft,
        float **Vz_now,float **Vz_pre,float **Vz_aft,
        float **K_now,float **K_pre,float **K_aft,
        float **Psa_now,float **Psa_pre,float **Psa_aft,
        float **Cita_now,float **Cita_pre,float **Cita_aft,
        float **V,float **Vv,float **DEN,
        float **absorbx,float **absorbz,struct PARAMETER *param,
        float wavelet)
{
    int i,j;


    for(i=8;i<param->Nx-8; i++)
    {
        for(j=8; j<param->PML; j++)
        {
			/*
			Cita_aft[i][j]=Cita_now[i][j]+param->dt*Psa_aft[i][j];
			*/
            Cita_aft[i][j]=1.0/(1.0+0.5*param->dt*absorbz[i][j])*(
                    (1.0-0.5*param->dt*absorbz[i][j])
					*Cita_now[i][j]+param->dt*Psa_aft[i][j]);
					
        }
    }


return 0.0;
}
