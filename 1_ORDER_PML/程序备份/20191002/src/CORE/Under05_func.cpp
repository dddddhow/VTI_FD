#include "../../lib/Under_func.h"

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
        float wavelet)
{
    int i,j;
	float temp=0;


    for(i=8;i<param->Nx-8; i++)
    {
        for(j=param->NZ+param->PML; j<param->Nz-8; j++)
        {
            Px_aft[i][j]=1.0/(1.0+0.5*param->dt*absorbx[i][j])*((1.0-0.5*param->dt*absorbx[i][j])*Px_now[i][j]
                    +param->dt*(1.0+2.0*param->yita)*V[i][j]*V[i][j]*K_now[i][j]);
        }
    }

    for(i=8;i<param->Nx-8; i++)
    {
        for(j=param->NZ+param->PML; j<param->Nz-8; j++)
        {
            Pz_aft[i][j]=1.0/(1.0+0.5*param->dt*absorbz[i][j])*((1.0-0.5*param->dt*absorbz[i][j])*Pz_now[i][j]
                    -Vv[i][j]*Vv[i][j]*DEN[i][j]*param->dt*1.0/param->dz*(
                        param->C1*(Vz_now[i][j]-Vz_now[i][j-1])+
                        param->C2*(Vz_now[i][j+1]-Vz_now[i][j-2])+
                        param->C3*(Vz_now[i][j+2]-Vz_now[i][j-3])+
                        param->C4*(Vz_now[i][j+3]-Vz_now[i][j-4])+
                        param->C5*(Vz_now[i][j+4]-Vz_now[i][j-5])+
                        param->C6*(Vz_now[i][j+5]-Vz_now[i][j-6])+
                        param->C7*(Vz_now[i][j+6]-Vz_now[i][j-7])+
                        param->C8*(Vz_now[i][j+7]-Vz_now[i][j-8]))-
                    2.0*param->yita*V[i][j]*V[i][j]*Vv[i][j]*Vv[i][j]*param->dt*1.0/param->dz*(
                        param->C1*(Cita_now[i][j]-Cita_now[i][j-1])+
                        param->C2*(Cita_now[i][j+1]-Cita_now[i][j-2])+
                        param->C3*(Cita_now[i][j+2]-Cita_now[i][j-3])+
                        param->C4*(Cita_now[i][j+3]-Cita_now[i][j-4])+
                        param->C5*(Cita_now[i][j+4]-Cita_now[i][j-5])+
                        param->C6*(Cita_now[i][j+5]-Cita_now[i][j-6])+
                        param->C7*(Cita_now[i][j+6]-Cita_now[i][j-7])+
                        param->C8*(Cita_now[i][j+7]-Cita_now[i][j-8])));
        }
    }
	
    for(i=8;i<param->Nx-8; i++)
    {
        for(j=param->NZ+param->PML; j<param->Nz-8; j++)
        {
			/*
			 temp=(1+2.0*param->yita)*V[i][j]*V[i][j]*K_now[i][j]-
                        Vv[i][j]*Vv[i][j]*DEN[i][j]*1.0/param->dz*(
                                param->C1*(Vz_now[i][j]-Vz_now[i][j-1])+
                                param->C2*(Vz_now[i][j+1]-Vz_now[i][j-2])+
                                param->C3*(Vz_now[i][j+2]-Vz_now[i][j-3])+
                                param->C4*(Vz_now[i][j+3]-Vz_now[i][j-4])+
                                param->C5*(Vz_now[i][j+4]-Vz_now[i][j-5])+
                                param->C6*(Vz_now[i][j+5]-Vz_now[i][j-6])+
                                param->C7*(Vz_now[i][j+6]-Vz_now[i][j-7])+
                                param->C8*(Vz_now[i][j+7]-Vz_now[i][j-8]))-
                        2.0*param->yita*V[i][j]*V[i][j]*Vv[i][j]*Vv[i][j]*1.0/param->dx*(
                                param->C1*(Cita_now[i][j]-Cita_now[i][j-1])+
                                param->C2*(Cita_now[i][j+1]-Cita_now[i][j-2])+
                                param->C3*(Cita_now[i][j+2]-Cita_now[i][j-3])+
                                param->C4*(Cita_now[i][j+3]-Cita_now[i][j-4])+
                                param->C5*(Cita_now[i][j+4]-Cita_now[i][j-5])+
                                param->C6*(Cita_now[i][j+5]-Cita_now[i][j-6])+
                                param->C7*(Cita_now[i][j+6]-Cita_now[i][j-7])+
                                param->C8*(Cita_now[i][j+7]-Cita_now[i][j-8]));


                    if(i==param->nx_location && j==param->nz_location)
                    {
                        temp=temp+wavelet;
                    }

                    P_aft[i][j]=P_now[i][j]+temp*param->dt;
			*/
            P_aft[i][j] = Px_aft[i][j] + Pz_aft[i][j];
            if(i==param->nx_location && j==param->nz_location)
            {
                P_aft[i][j]=P_aft[i][j]+wavelet*param->dt;
            }		
					
        }
    }	
	
	return 0.0;
}
