#include "../lib/Under_func.h"

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
        float wavelet)
{
    int i,j;

    for(i=8;i<param->Nx-8; i++)
    {
        for(j=param->NZ+param->PML; j<param->Nz-8; j++)
        {
			/*
			Vx_now[i][j]=Vx_pre[i][j]-1.0/DEN[i][j]*(param->dt*1.0/param->dx)*(
                            param->C1*(P_now[i+1][j]-P_now[i][j])+
                            param->C2*(P_now[i+2][j]-P_now[i-1][j])+
                            param->C3*(P_now[i+3][j]-P_now[i-2][j])+
                            param->C4*(P_now[i+4][j]-P_now[i-3][j])+
                            param->C5*(P_now[i+5][j]-P_now[i-4][j])+
                            param->C6*(P_now[i+6][j]-P_now[i-5][j])+
                            param->C7*(P_now[i+7][j]-P_now[i-6][j])+
                            param->C8*(P_now[i+8][j]-P_now[i-7][j]));
							*/
            Vx_now[i][j]=1.0/(1.0+0.5*param->dt*absorbx[i][j])*(
                    (1.0-0.5*param->dt*absorbx[i][j])*Vx_pre[i][j]
                    -1.0/DEN[i][j]*(param->dt*1.0/param->dx)*(
                        param->C1*(P_now[i+1][j]-P_now[i][j])+
                        param->C2*(P_now[i+2][j]-P_now[i-1][j])+
                        param->C3*(P_now[i+3][j]-P_now[i-2][j])+
                        param->C4*(P_now[i+4][j]-P_now[i-3][j])+
                        param->C5*(P_now[i+5][j]-P_now[i-4][j])+
                        param->C6*(P_now[i+6][j]-P_now[i-5][j])+
                        param->C7*(P_now[i+7][j]-P_now[i-6][j])+
                        param->C8*(P_now[i+8][j]-P_now[i-7][j])));
						
        }
    }

    for(i=8;i<param->Nx-8; i++)
    {
        for(j=param->NZ+param->PML; j<param->Nz-8; j++)
        {
			/*
			Vz_now[i][j]=Vz_pre[i][j]-(1.0/DEN[i][j])*(param->dt*1.0/param->dz)*(
                            param->C1*(P_now[i][j+1]-P_now[i][j])+
                            param->C2*(P_now[i][j+2]-P_now[i][j-1])+
                            param->C3*(P_now[i][j+3]-P_now[i][j-2])+
                            param->C4*(P_now[i][j+4]-P_now[i][j-3])+
                            param->C5*(P_now[i][j+5]-P_now[i][j-4])+
                            param->C6*(P_now[i][j+6]-P_now[i][j-5])+
                            param->C7*(P_now[i][j+7]-P_now[i][j-6])+
                            param->C8*(P_now[i][j+8]-P_now[i][j-7]));
							*/
            Vz_now[i][j]=1.0/(1.0+0.5*param->dt*absorbz[i][j])*(
                    (1.0-0.5*param->dt*absorbz[i][j])*Vz_pre[i][j]
                    -(1.0/DEN[i][j])*(param->dt*1.0/param->dz)*(
                        param->C1*(P_now[i][j+1]-P_now[i][j])+
                        param->C2*(P_now[i][j+2]-P_now[i][j-1])+
                        param->C3*(P_now[i][j+3]-P_now[i][j-2])+
                        param->C4*(P_now[i][j+4]-P_now[i][j-3])+
                        param->C5*(P_now[i][j+5]-P_now[i][j-4])+
                        param->C6*(P_now[i][j+6]-P_now[i][j-5])+
                        param->C7*(P_now[i][j+7]-P_now[i][j-6])+
                        param->C8*(P_now[i][j+8]-P_now[i][j-7])));
						

        }
    }


	return 0.0;
}
