#include "../lib/Model_func.h"

double model_func(struct PARAMETER *param, float **V, float **Vv, float **DEN)
{
    /********************参数定义********************/
    int i, j;
    FILE  *fp;
    /********************模型赋值********************/
 
    for(i=0;i<param->Nx;i++)
    {
        for(j=0;j<param->Nz;j++)
        {
            Vv[i][j]=2000.0;
            V[i][j]=Vv[i][j]*sqrt(1.0+2.0*param->delt);
            DEN[i][j]=1.0;
        }
    }

    for(i=0;i<param->Nx;i++)
    {
        for(j=param->PML+param->NZ/2;j<param->Nz;j++)
        {
            Vv[i][j]=4000.0;
            V[i][j]=Vv[i][j]*sqrt(1.0+2.0*param->delt);
        }
    }
/*
    for(i=0;i<param->Nx;i++)
    {
        for(j=0;j<param->Nz;j++)
        {
        Vv[i][j]=2000.0;
        V[i][j]=Vv[i][j]*sqrt(1.0+2.0*param->delt);
        DEN[i][j]=1.0;
        }
    }

    for(i=0;i<param->Nx;i++)
    {
        for(j=param->PML+40;j<param->Nz;j++)
        {
            Vv[i][j]=3000.0;
            V[i][j]=Vv[i][j]*sqrt(1.0+2.0*param->delt);
        }
    }

    for(i=0;i<param->Nx;i++)
    {
        for(j=param->PML+80;j<param->Nz;j++)
        {
            Vv[i][j]=2500.0;
            V[i][j]=Vv[i][j]*sqrt(1.0+2.0*param->delt);
        }
    }

	for(i=0;i<param->Nx;i++)
    {
        for(j=param->PML+160;j<param->Nz;j++)
        {
            Vv[i][j]=3500.0;
            V[i][j]=Vv[i][j]*sqrt(1.0+2.0*param->delt);
        }
    }
*/
    /********************模型保存********************/
    if((fp = fopen ("../file/Vp.dat", "wb"))!=NULL)
    {

        for (i=param->PML;i<param->Nx-param->PML;i++)
        {
            for (j=param->PML;j<param->Nz-param->PML;j++)
            {
                fwrite (&V[i][j] , sizeof(float), 1, fp);

            }
        }
        fclose (fp);
    }

}
