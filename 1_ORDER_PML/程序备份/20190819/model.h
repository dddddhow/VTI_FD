#ifndef MODEL_H
#define MODEL_H

double model(int Nx,int Nz,float **V,float **Vv,float **DEN,float epsilon,float delt,float yita)
{
    /********************参数定义********************/
    int i,j;
    FILE  *fp;
    /********************模型赋值********************/

    for(i=0; i<Nx; i++)
    {

        for(j=0; j<PML+100;j++)
        {
            Vv[i][j]=3000.0;
        }

        for(j=PML+100; j<PML+150;j++)
        {
            Vv[i][j]=4000.0;
        }

        for(j=PML+150; j<PML+300;j++)
        {
            Vv[i][j]=5000.0;
        }

        for(j=PML+300; j<Nz;j++)
        {
            Vv[i][j]=5000.0;
        }

        for(j=0; j<Nz;j++)
        {
            V[i][j]=Vv[i][j]*sqrt(1.0+2.0*delt);
            DEN[i][j]=2400.0;
        }

    }
    /********************模型保存********************/
    if((fp = fopen ("../file/Vp.dat", "wb"))!=NULL)
    {

        for (i=PML;i<Nx-PML;i++)
        {
            for (j=PML;j<Nz-PML;j++)
            {
                fwrite (&V[i][j] , sizeof(float), 1, fp);

            }
        }
        fclose (fp);
    }


}
#endif
