#include "../../lib/Left_func.h"

double left02_func(float **P_now,float **P_pre,float **P_aft,
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
    for(i=8; i<param->PML; i++)
    {
        for(j=8; j<param->Nz-8; j++)
        {
			/*
			K_now[i][j]=-DEN[i][j]*1.0/param->dx*(
                            param->C1*(Vx_now[i][j]-Vx_now[i-1][j])+
                            param->C2*(Vx_now[i+1][j]-Vx_now[i-2][j])+
                            param->C3*(Vx_now[i+2][j]-Vx_now[i-3][j])+
                            param->C4*(Vx_now[i+3][j]-Vx_now[i-4][j])+
                            param->C5*(Vx_now[i+4][j]-Vx_now[i-5][j])+
                            param->C6*(Vx_now[i+5][j]-Vx_now[i-6][j])+
                            param->C7*(Vx_now[i+6][j]-Vx_now[i-7][j])+
                            param->C8*(Vx_now[i+7][j]-Vx_now[i-8][j]));
							*/
            K_now[i][j]=-DEN[i][j]*1.0/param->dx*(
                    param->C1*(Vx_now[i][j]-Vx_now[i-1][j])+
                    param->C2*(Vx_now[i+1][j]-Vx_now[i-2][j])+
                    param->C3*(Vx_now[i+2][j]-Vx_now[i-3][j])+
                    param->C4*(Vx_now[i+3][j]-Vx_now[i-4][j])+
                    param->C5*(Vx_now[i+4][j]-Vx_now[i-5][j])+
                    param->C6*(Vx_now[i+5][j]-Vx_now[i-6][j])+
                    param->C7*(Vx_now[i+6][j]-Vx_now[i-7][j])+
                    param->C8*(Vx_now[i+7][j]-Vx_now[i-8][j]));
					
        }
    }

return 0.0;

}
