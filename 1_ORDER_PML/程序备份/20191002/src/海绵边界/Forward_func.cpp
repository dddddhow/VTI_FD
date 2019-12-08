#include "../lib/Forward_func.h"


double forward_func(struct PARAMETER *param)
{

    //---------------------------------------------//
    //              变量申请及赋值                 //
    //---------------------------------------------//
    FILE  *fpp, *fpr;
    float **P_now,**P_pre,**P_aft,**Vx_now,**Vx_pre,**Vx_aft,
          **Vz_now,**Vz_pre,**Vz_aft;
    float **K_now,**K_pre,**K_aft,**Psa_now,**Psa_pre,**Psa_aft,**Cita_now,
          **Cita_pre,**Cita_aft;
    float **Px_now, **Px_pre, **Px_aft;
    float **Pz_now, **Pz_pre, **Pz_aft;
    float **V,**Vv,**DEN,**record;
	float **absorbx,**absorbz;
    float *signal;

    float wavelet;
    float temp=0.0;
    int i, ii, j, jj, k;

    //*开辟动态数组*//
    V  = new float* [param->Nx];
    Vv = new float* [param->Nx];
    DEN = new float* [param->Nx];

    P_now  = new float* [param->Nx];
    P_pre  = new float* [param->Nx];
    P_aft  = new float* [param->Nx];

    Px_now = new float* [param->Nx];
    Pz_now = new float* [param->Nx];
    Px_pre = new float* [param->Nx];
    Pz_pre = new float* [param->Nx];
    Px_aft = new float* [param->Nx];
    Pz_aft = new float* [param->Nx];

    Vx_now = new float* [param->Nx];
    Vx_aft = new float* [param->Nx];
    Vx_pre = new float* [param->Nx];

    Vz_now = new float* [param->Nx];
    Vz_aft = new float* [param->Nx];
    Vz_pre = new float* [param->Nx];

    K_now = new float* [param->Nx];
    K_pre = new float* [param->Nx];
    K_aft = new float* [param->Nx];

    Psa_now = new float* [param->Nx];
    Psa_pre = new float* [param->Nx];
    Psa_aft = new float* [param->Nx];

    Cita_now = new float* [param->Nx];
    Cita_pre = new float* [param->Nx];
    Cita_aft = new float* [param->Nx];
    record  = new float* [param->Nx];

	absorbx  = new float* [param->Nx];
	absorbz  = new float* [param->Nx];

    for(i=0;i<param->Nx;i++){V[i]=new float [param->Nz];}
    for(i=0;i<param->Nx;i++){Vv[i]=new float [param->Nz];}
    for(i=0;i<param->Nx;i++){DEN[i]=new float [param->Nz];}


    for(i=0;i<param->Nx;i++){P_now[i]=new float [param->Nz];}
    for(i=0;i<param->Nx;i++){P_pre[i]=new float [param->Nz];}
    for(i=0;i<param->Nx;i++){P_aft[i]=new float [param->Nz];}

    for(i=0;i<param->Nx;i++){Px_now[i]=new float [param->Nz];}
    for(i=0;i<param->Nx;i++){Px_pre[i]=new float [param->Nz];}
    for(i=0;i<param->Nx;i++){Px_aft[i]=new float [param->Nz];}

	for(i=0;i<param->Nx;i++){Pz_now[i]=new float [param->Nz];}
    for(i=0;i<param->Nx;i++){Pz_pre[i]=new float [param->Nz];}
    for(i=0;i<param->Nx;i++){Pz_aft[i]=new float [param->Nz];}

    for(i=0;i<param->Nx;i++){Vx_now[i]=new float [param->Nz];}
    for(i=0;i<param->Nx;i++){Vx_pre[i]=new float [param->Nz];}
    for(i=0;i<param->Nx;i++){Vx_aft[i]=new float [param->Nz];}

    for(i=0;i<param->Nx;i++){Vz_now[i]=new float [param->Nz];}
    for(i=0;i<param->Nx;i++){Vz_pre[i]=new float [param->Nz];}
    for(i=0;i<param->Nx;i++){Vz_aft[i]=new float [param->Nz];}

    for(i=0;i<param->Nx;i++){K_now[i]=new float [param->Nz];}
    for(i=0;i<param->Nx;i++){K_pre[i]=new float [param->Nz];}
    for(i=0;i<param->Nx;i++){K_aft[i]=new float [param->Nz];}

    for(i=0;i<param->Nx;i++){Psa_now[i]=new float [param->Nz];}
    for(i=0;i<param->Nx;i++){Psa_pre[i]=new float [param->Nz];}
    for(i=0;i<param->Nx;i++){Psa_aft[i]=new float [param->Nz];}

    for(i=0;i<param->Nx;i++){Cita_now[i]=new float [param->Nz];}
    for(i=0;i<param->Nx;i++){Cita_pre[i]=new float [param->Nz];}
    for(i=0;i<param->Nx;i++){Cita_aft[i]=new float [param->Nz];}

    for(i=0;i<param->Nx;i++){record[i]=new float [param->Nt];}

	for(i=0;i<param->Nx;i++){absorbx[i]=new float [param->Nz];}
	for(i=0;i<param->Nx;i++){absorbz[i]=new float [param->Nz];}

    //*初始化*//

       for(i=0; i<param->Nx; i++)
       {
       for(j=0; j<param->Nz;j++)
       {
       P_now[i][j]=0.0;
       P_pre[i][j]=0.0;
       P_aft[i][j]=0.0;

       Px_now[i][j]=0.0;
       Px_pre[i][j]=0.0;
       Px_aft[i][j]=0.0;

       Pz_now[i][j]=0.0;
       Pz_pre[i][j]=0.0;
       Pz_aft[i][j]=0.0;

       Vx_now[i][j]=0.0;
       Vx_pre[i][j]=0.0;
       Vx_aft[i][j]=0.0;

       Vz_now[i][j]=0.0;
       Vz_pre[i][j]=0.0;
       Vz_aft[i][j]=0.0;

       K_now[i][j]=0.0;
       K_pre[i][j]=0.0;
       K_aft[i][j]=0.0;

       Psa_now[i][j]=0.0;
       Psa_pre[i][j]=0.0;
       Psa_aft[i][j]=0.0;

       Cita_now[i][j]=0.0;
       Cita_pre[i][j]=0.0;
       Cita_aft[i][j]=0.0;
       }
       }

       for(i=0; i<param->Nx; i++)
       {
       for(j=0; j<param->Nt;j++)
       {
       record[i][j]=0.0;
       }
       }

    for(i=0; i<param->Nx; i++)
    {
        for(j=0; j<param->Nz;j++)
        {
            absorbx[i][j]=0.0;//吸收衰减函数初始值为1.0
			absorbz[i][j]=0.0;//吸收衰减函数初始值为1.0
        }
    }

    //------------------------------------//
    //            模型设计                //
    //------------------------------------//

    model_func(param, V, Vv, DEN);

    //------------------------------------//
    //            衰减模型                //
    //------------------------------------//
/*
    //左边界//
    for(i=0;i<param->PML;i++)
    {
        for(j=0;j<param->Nz;j++)
        {
            absorbx[i][j]=3.0*V[i][j]*1.0/(2.0*param->PML)*log(1.0/0.001)
                *((param->PML-i)*1.0/param->PML)*((param->PML-i)*1.0/param->PML);
        }
    }
    //右边界//
    for(i=param->NX+param->PML;i<param->Nx;i++)
    {
        for(j=0;j<param->Nz;j++)
        {
            absorbx[i][j]=3.0*V[i][j]*1.0/(2.0*param->PML)*log(1.0/0.001)
                *((param->Nx-i)*1.0/param->PML)*((param->Nx-i)*1.0/param->PML);
        }
    }
    //上边界//
    for(i=0;i<param->Nx;i++)
    {
        for(j=0;j<param->PML;j++)
        {
            absorbz[i][j]=3.0*V[i][j]*1.0/(2.0*param->PML)*log(1.0/0.001)
                *((param->PML-j)*1.0/param->PML)*((param->PML-j)*1.0/param->PML);
        }
    }
    //下边界//
    for(i=0;i<param->Nx;i++)
    {
        for(j=param->NZ+param->PML;j<param->Nz;j++)
        {
            absorbz[i][j]=3.0*V[i][j]*1.0/(2.0*param->PML)*log(1.0/0.001)
                *((param->Nz-i)*1.0/param->PML)*((param->Nz-i)*1.0/param->PML);
        }
    }

*/
    //------------------------------------//
    //            子波设计                //
    //------------------------------------//

    signal = new float [param->Lw];
    rickerwavelet_func(param, signal);

    //-----------------------------------------//
    //          进度条显示设置                 //
    //-----------------------------------------//
    int i_timebar = 0;
    char timebar[52];
    const char *timebar_lable = "|/-\\";
    timebar[0] = 0;
    //-----------------------------------------//
    //          炮点位置                       //
    //-----------------------------------------//

    int num_shot=-1;
    int ds=floor(param->NX / (param->Ns-1));

    //-----------------------------------------//
    //          炮循环                         //
    //-----------------------------------------//
	/*
    for(param->nx_location = param->PML;
            param->nx_location <= param->PML + param->NX;
            param->nx_location=param->nx_location+ds)
    {
		*/

        //初始化
        for(i=0; i<param->Nx; i++)
        {
            for(j=0; j<param->Nz;j++)
            {
                P_now[i][j]=0.0;
                P_pre[i][j]=0.0;
                P_aft[i][j]=0.0;
                Px_now[i][j]=0.0;
                Px_pre[i][j]=0.0;
                Px_aft[i][j]=0.0;
                Pz_now[i][j]=0.0;
                Pz_pre[i][j]=0.0;
                Pz_aft[i][j]=0.0;
                Vx_now[i][j]=0.0;
                Vx_pre[i][j]=0.0;
                Vx_aft[i][j]=0.0;
                Vz_now[i][j]=0.0;
                Vz_pre[i][j]=0.0;
                Vz_aft[i][j]=0.0;
                K_now[i][j]=0.0;
                K_pre[i][j]=0.0;
                K_aft[i][j]=0.0;
                Psa_now[i][j]=0.0;
                Psa_pre[i][j]=0.0;
                Psa_aft[i][j]=0.0;
                Cita_now[i][j]=0.0;
                Cita_pre[i][j]=0.0;
                Cita_aft[i][j]=0.0;
            }
        }
        for(i=0; i<param->Nx; i++)
        {
            for(j=0; j<param->Nt;j++)
            {
                record[i][j]=0.0;
            }
        }


        //炮术显示设置
        num_shot++;//炮数
        printf("第%d/%d炮\t",num_shot+1,param->Ns);
        printf("炮点位置：[%f米,%f米]\n",(param->nx_location-param->PML)
                *param->dx,(param->nz_location-param->PML)*param->dz);


        //-----------------------------------------//
        //            波场时间递推计算             //
        //----------------------------------------//

        for(k=0; k<param->Nt; k++) //波场时间递推计算开始
        {

            //*进度条显示*
            i_timebar=k*50/param->Nt;
            printf("[%-50s][%d%%][%c]\r", timebar, (i_timebar+1)*2, timebar_lable[k%4]);
            fflush(stdout);
            timebar[i_timebar] = '#';
            timebar[i_timebar+1] = 0;

            //*震源积分*//
            wavelet=0;
            if(k<param->Lw)  //震源的长度param->Lw
            {
                for(ii=0;ii<k;ii++)
                {
                    wavelet=wavelet+signal[ii]*param->dt;
                }
            }

            //*计算P*//

			//左边//
			left_func(P_now,P_pre,P_aft,
			     Px_now,Px_pre,Px_aft,
				 Pz_now,Pz_pre,Pz_aft,
				 Vx_now,Vx_pre,Vx_aft,
				 Vz_now,Vz_pre,Vz_aft,
				 K_now,K_pre,K_aft,
				 Psa_now,Psa_pre,Psa_aft,
				 Cita_now,Cita_pre,Cita_aft,
				 V,Vv,DEN,
				 absorbx,absorbz,param,
				 wavelet);
			//右边//
			right_func(P_now,P_pre,P_aft,
			     Px_now,Px_pre,Px_aft,
				 Pz_now,Pz_pre,Pz_aft,
				 Vx_now,Vx_pre,Vx_aft,
				 Vz_now,Vz_pre,Vz_aft,
				 K_now,K_pre,K_aft,
				 Psa_now,Psa_pre,Psa_aft,
				 Cita_now,Cita_pre,Cita_aft,
				 V,Vv,DEN,
				 absorbx,absorbz,param,
				 wavelet);
			//上边//
			up_func(P_now,P_pre,P_aft,
			     Px_now,Px_pre,Px_aft,
				 Pz_now,Pz_pre,Pz_aft,
				 Vx_now,Vx_pre,Vx_aft,
				 Vz_now,Vz_pre,Vz_aft,
				 K_now,K_pre,K_aft,
				 Psa_now,Psa_pre,Psa_aft,
				 Cita_now,Cita_pre,Cita_aft,
				 V,Vv,DEN,
				 absorbx,absorbz,param,
				 wavelet);
			//下边//
			under_func(P_now,P_pre,P_aft,
			     Px_now,Px_pre,Px_aft,
				 Pz_now,Pz_pre,Pz_aft,
				 Vx_now,Vx_pre,Vx_aft,
				 Vz_now,Vz_pre,Vz_aft,
				 K_now,K_pre,K_aft,
				 Psa_now,Psa_pre,Psa_aft,
				 Cita_now,Cita_pre,Cita_aft,
				 V,Vv,DEN,
				 absorbx,absorbz,param,
				 wavelet);

			//中间
			center_func(P_now,P_pre,P_aft,
			     Px_now,Px_pre,Px_aft,
				 Pz_now,Pz_pre,Pz_aft,
				 Vx_now,Vx_pre,Vx_aft,
				 Vz_now,Vz_pre,Vz_aft,
				 K_now,K_pre,K_aft,
				 Psa_now,Psa_pre,Psa_aft,
				 Cita_now,Cita_pre,Cita_aft,
				 V,Vv,DEN,
				 absorbx,absorbz,param,
				 wavelet);


            //*时间片波场值替换*//
            for(i=0; i<param->Nx; i++)
            {
                for(j=0; j<param->Nz; j++)
                {
                    Vz_pre[i][j]=Vz_now[i][j];
                    Vx_pre[i][j]=Vx_now[i][j];
                    P_pre[i][j]=P_now[i][j];
                    P_now[i][j]=P_aft[i][j];
                    Px_pre[i][j]=Px_now[i][j];
                    Px_now[i][j]=Px_aft[i][j];
                    Pz_pre[i][j]=Pz_now[i][j];
                    Pz_now[i][j]=Pz_aft[i][j];
                    Psa_now[i][j]=Psa_aft[i][j];
                    Cita_now[i][j]=Cita_aft[i][j];
                }
            }

            //-------------------------------------------//
            //               时间切片                    //
            //-------------------------------------------//
        if(k%20==0)
        {
            if((fpr = fopen ("../file/P", "a+"))!=NULL)
            {
                for(ii=param->PML; ii<param->Nx - param->PML; ii++)
                    for(jj=param->PML; jj<param->Nz - param->PML; jj++)
                    {
                        fwrite (&P_now[ii][jj] , sizeof(float), 1, fpr);

                    }
            }
            fclose (fpr);
        }


            //-------------------------------------------//
            //               地震记录record              //
            //-------------------------------------------//


            for(i=param->PML; i<param->Nx-param->PML; i++)
            {
                record[i][k]=P_now[i][param->nz_location];
            }

        }
        printf("\n");

        //*波场时间递推计算结束*//

        //-------------------------------------------//
        //       存储地表地震记录record              //
        //-------------------------------------------//
        IO_func(param, record);
		
    //}



    printf("\n计算结束啦\n");
    //---------------------------------------------//
    //                 释放内存                    //
    //---------------------------------------------//
    for(i=0;i<param->Nx;i++)
        delete []	V[i];
    delete []V;
    for(i=0;i<param->Nx;i++)
        delete []	Vv[i];
    delete []Vv;
    for(i=0;i<param->Nx;i++)
        delete []	DEN[i];
    delete []DEN;
    for(i=0;i<param->Nx;i++)
        delete []	record[i];
    delete []record;

    for(i=0;i<param->Nx;i++)
        delete []	P_now[i];
    delete []P_now;
    for(i=0;i<param->Nx;i++)
        delete []	P_pre[i];
    delete []P_pre;
    for(i=0;i<param->Nx;i++)
        delete []	P_aft[i];
    delete []P_aft;

    for(i=0;i<param->Nx;i++)
        delete []	Vx_now[i];
    delete []Vx_now;
    for(i=0;i<param->Nx;i++)
        delete []	Vx_pre[i];
    delete []Vx_pre;
    for(i=0;i<param->Nx;i++)
        delete []	Vx_aft[i];
    delete []Vx_aft;

    for(i=0;i<param->Nx;i++)
        delete []	Vz_now[i];
    delete []Vz_now;
    for(i=0;i<param->Nx;i++)
        delete []	Vz_aft[i];
    delete []Vz_aft;
    for(i=0;i<param->Nx;i++)
        delete []	Vz_pre[i];
    delete []Vz_pre;

    for(i=0;i<param->Nx;i++)
        delete []	K_now[i];
    delete []K_now;
    for(i=0;i<param->Nx;i++)
        delete []	K_aft[i];
    delete []K_aft;
    for(i=0;i<param->Nx;i++)
        delete []	K_pre[i];
    delete []K_pre;

    for(i=0;i<param->Nx;i++)
        delete []	Cita_now[i];
    delete []Cita_now;
    for(i=0;i<param->Nx;i++)
        delete []	Cita_aft[i];
    delete []Cita_aft;
    for(i=0;i<param->Nx;i++)
        delete []	Cita_pre[i];
    delete []Cita_pre;

    for(i=0;i<param->Nx;i++)
        delete []	Psa_now[i];
    delete []Psa_now;
    for(i=0;i<param->Nx;i++)
        delete []	Psa_aft[i];
    delete []Psa_aft;
    for(i=0;i<param->Nx;i++)
        delete []	Psa_pre[i];
    delete []Psa_pre;

    return 0;
}


