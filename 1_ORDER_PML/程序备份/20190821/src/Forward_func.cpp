#include "../lib/Forward_func.h"

double forward(PARAMETER* param)
{
    //---------------------------------------------//
    //              变量申请及赋值                 //
    //---------------------------------------------//
    float **P_now,**P_pre,**P_aft,**Vx_now,**Vx_pre,**Vx_aft,**Vz_now,**Vz_pre,**Vz_aft;
    float **K_now,**K_pre,**K_aft,**Psa_now,**Psa_pre,**Psa_aft,**Cita_now,**Cita_pre,**Cita_aft;
    float **V,**Vv,**DEN,**record;
	float *signal;

    float absorb[Nx][Nz],wavelet;
    float temp=0.0;
    int i, ii, j, k;
    //*开辟动态数组*//

    V  = new float* [Nx];
    Vv = new float* [Nx];
    DEN = new float* [Nx];

    P_now  = new float* [Nx];
    P_pre  = new float* [Nx];
    P_aft  = new float* [Nx];

    Vx_now = new float* [Nx];
    Vx_aft = new float* [Nx];
    Vx_pre = new float* [Nx];

    Vz_now = new float* [Nx];
    Vz_aft = new float* [Nx];
    Vz_pre = new float* [Nx];

    K_now = new float* [Nx];
    K_pre = new float* [Nx];
    K_aft = new float* [Nx];

    Psa_now = new float* [Nx];
    Psa_pre = new float* [Nx];
    Psa_aft = new float* [Nx];

    Cita_now = new float* [Nx];
    Cita_pre = new float* [Nx];
    Cita_aft = new float* [Nx];

    record  = new float* [Nx];


    for(i=0;i<Nx;i++){V[i]=new float [Nz];}
    for(i=0;i<Nx;i++){Vv[i]=new float [Nz];}
    for(i=0;i<Nx;i++){DEN[i]=new float [Nz];}

    for(i=0;i<Nx;i++){P_now[i]=new float [Nz];}
    for(i=0;i<Nx;i++){P_pre[i]=new float [Nz];}
    for(i=0;i<Nx;i++){P_aft[i]=new float [Nz];}

    for(i=0;i<Nx;i++){Vx_now[i]=new float [Nz];}
    for(i=0;i<Nx;i++){Vx_pre[i]=new float [Nz];}
    for(i=0;i<Nx;i++){Vx_aft[i]=new float [Nz];}

    for(i=0;i<Nx;i++){Vz_now[i]=new float [Nz];}
    for(i=0;i<Nx;i++){Vz_pre[i]=new float [Nz];}
    for(i=0;i<Nx;i++){Vz_aft[i]=new float [Nz];}

    for(i=0;i<Nx;i++){K_now[i]=new float [Nz];}
    for(i=0;i<Nx;i++){K_pre[i]=new float [Nz];}
    for(i=0;i<Nx;i++){K_aft[i]=new float [Nz];}

    for(i=0;i<Nx;i++){Psa_now[i]=new float [Nz];}
    for(i=0;i<Nx;i++){Psa_pre[i]=new float [Nz];}
    for(i=0;i<Nx;i++){Psa_aft[i]=new float [Nz];}

    for(i=0;i<Nx;i++){Cita_now[i]=new float [Nz];}
    for(i=0;i<Nx;i++){Cita_pre[i]=new float [Nz];}
    for(i=0;i<Nx;i++){Cita_aft[i]=new float [Nz];}

    for(i=0;i<Nx;i++){record[i]=new float [Nt];}

    //*初始化*//
    for(i=0; i<Nx; i++)
    {
        for(j=0; j<Nz;j++)
        {
            P_now[i][j]=0.0;
            P_pre[i][j]=0.0;
            P_aft[i][j]=0.0;

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

    for(i=0; i<Nx; i++)
    {
        for(j=0; j<Nt;j++)
        {
            record[i][j]=0.0;
        }
    }

    for(i=0; i<Nx; i++)
    {
        for(j=0; j<Nz;j++)
        {
            absorb[i][j]=1.0;//吸收衰减函数初始值为1.0
        }
    }

    //------------------------------------//
    //            模型设计                //
    //------------------------------------//
    model(NX, NZ, PML,V, Vv, DEN, epsilon, delt, yita);

    //------------------------------------//
    //            衰减模型                //
    //------------------------------------//
    //左边界//
    for(i=PML;i<Nx-PML;i++)
    {
        for(j=0;j<PML;j++)
        {
            absorb[i][j]=-0.38*(PML-j-1)*(PML-j-1)/(1.0*PML*PML)+1;
        }
    }
    //右边界//
    for(i=PML;i<Nx-PML;i++)
    {
        for(j=Nz-PML;j<Nz;j++)
        {
            absorb[i][j]=-0.38*(j-Nz+PML)*(j-Nz+PML)/(1.0*PML*PML)+1;
        }
    }
    //上边界//
    for(i=0;i<PML;i++)
    {
        for(j=0;j<Nz;j++)
        {
            absorb[i][j]=-0.38*(PML-i-1)*(PML-i-1)/(1.0*PML*PML)+1;
        }
    }
    //下边界//
    for(i=Nx-PML;i<Nx;i++)
    {
        for(j=0;j<Nx;j++)
        {
            absorb[i][j]=-0.38*(i-Nx+PML)*(i-Nx+PML)/(1.0*PML*PML)+1;
        }
    }


    //------------------------------------//
    //            子波设计                //
    //------------------------------------//
	signal = new float [Lw];
	ricker_wavelet(Lw, freq, dt, signal);

    //-----------------------------------------//
    //          进度条显示设置                 //
    //-----------------------------------------//
    int i_timebar = 0;
    char timebar[102];
    const char *timebar_lable = "|/-\\";
    timebar[0] = 0;

    //-----------------------------------------//
    //          炮点位置                       //
    //-----------------------------------------//
    int num_shot=-1;
    int ds=NX/Ns;

    for(nx_location=PML;nx_location<PML+NX;nx_location=nx_location+ds)
    {
        num_shot++;//炮数
        printf("第%d/%d炮\t",num_shot+1,Ns);
        printf("炮点位置：[%f米,%f米]\n",(nx_location-PML)*dx,(nz_location-PML)*dz);

        //-----------------------------------------//
        //            波场时间递推计算             //
        //----------------------------------------//

        for(k=0; k<Nt; k++) //波场时间递推计算开始
        {
            //*进度条显示*//
            i_timebar=k*100/Nt;
            printf("[%-100s][%d%%][%c]\r", timebar, i_timebar+1, timebar_lable[k%4]);
            fflush(stdout);
            timebar[i_timebar] = '#';
            timebar[i_timebar+1] = 0;
            usleep(100000);



            //*震源积分*//
            wavelet=0;
            if(k<Lw)  //震源的长度Lw
            {
                for(ii=0;ii<k;ii++)
                {
                    wavelet=wavelet+signal[ii]*dt;
                }
            }

            //*计算P*//
            for(i=8; i<Nx-8; i++)
            {
                for(j=8; j<Nz-8; j++)
                {
                    Vx_now[i][j]=Vx_pre[i][j]-1.0/DEN[i][j]*(dt*1.0/dx)*(
                            C1*(P_now[i+1][j]-P_now[i][j])+
                            C2*(P_now[i+2][j]-P_now[i-1][j])+
                            C3*(P_now[i+3][j]-P_now[i-2][j])+
                            C4*(P_now[i+4][j]-P_now[i-3][j])+
                            C5*(P_now[i+5][j]-P_now[i-4][j])+
                            C6*(P_now[i+6][j]-P_now[i-5][j])+
                            C7*(P_now[i+7][j]-P_now[i-6][j])+
                            C8*(P_now[i+8][j]-P_now[i-7][j]));
                }
            }

            for(i=8; i<Nx-8; i++)
            {
                for(j=8; j<Nz-8; j++)
                {
                    Vz_now[i][j]=Vz_pre[i][j]-(1.0/DEN[i][j])*(dt*1.0/dz)*(
                            C1*(P_now[i][j+1]-P_now[i][j])+
                            C2*(P_now[i][j+2]-P_now[i][j-1])+
                            C3*(P_now[i][j+3]-P_now[i][j-2])+
                            C4*(P_now[i][j+4]-P_now[i][j-3])+
                            C5*(P_now[i][j+5]-P_now[i][j-4])+
                            C6*(P_now[i][j+6]-P_now[i][j-5])+
                            C7*(P_now[i][j+7]-P_now[i][j-6])+
                            C8*(P_now[i][j+8]-P_now[i][j-7]));
                }
            }

            for(i=8; i<Nx-8; i++)
            {
                for(j=8; j<Nz-8; j++)
                {

                    K_now[i][j]=-DEN[i][j]*1.0/dx*(
                            C1*(Vx_now[i][j]-Vx_now[i-1][j])+
                            C2*(Vx_now[i+1][j]-Vx_now[i-2][j])+
                            C3*(Vx_now[i+2][j]-Vx_now[i-3][j])+
                            C4*(Vx_now[i+3][j]-Vx_now[i-4][j])+
                            C5*(Vx_now[i+4][j]-Vx_now[i-5][j])+
                            C6*(Vx_now[i+5][j]-Vx_now[i-6][j])+
                            C7*(Vx_now[i+6][j]-Vx_now[i-7][j])+
                            C8*(Vx_now[i+7][j]-Vx_now[i-8][j]));
                }
            }

            for(i=8; i<Nx-8; i++)
            {
                for(j=8; j<Nz-8; j++)
                {
                    Psa_aft[i][j]=Psa_now[i][j]+dt*1.0/dz*(
                            C1*(K_now[i][j+1]-K_now[i][j])+
                            C2*(K_now[i][j+2]-K_now[i][j-1])+
                            C3*(K_now[i][j+3]-K_now[i][j-2])+
                            C4*(K_now[i][j+4]-K_now[i][j-3])+
                            C5*(K_now[i][j+5]-K_now[i][j-4])+
                            C6*(K_now[i][j+6]-K_now[i][j-5])+
                            C7*(K_now[i][j+7]-K_now[i][j-6])+
                            C8*(K_now[i][j+8]-K_now[i][j-7]));
                }
            }

            for(i=8; i<Nx-8; i++)
            {
                for(j=8; j<Nz-8; j++)
                {
                    Cita_aft[i][j]=Cita_now[i][j]+dt*Psa_aft[i][j];
                }
            }

            for(i=8; i<Nx-8; i++)
            {
                for(j=8; j<Nz-8; j++)
                {

                    temp=(1+2.0*yita)*V[i][j]*V[i][j]*K_now[i][j]-
                        Vv[i][j]*Vv[i][j]*DEN[i][j]*1.0/dz*(
                                C1*(Vz_now[i][j]-Vz_now[i][j-1])+
                                C2*(Vz_now[i][j+1]-Vz_now[i][j-2])+
                                C3*(Vz_now[i][j+2]-Vz_now[i][j-3])+
                                C4*(Vz_now[i][j+3]-Vz_now[i][j-4])+
                                C5*(Vz_now[i][j+4]-Vz_now[i][j-5])+
                                C6*(Vz_now[i][j+5]-Vz_now[i][j-6])+
                                C7*(Vz_now[i][j+6]-Vz_now[i][j-7])+
                                C8*(Vz_now[i][j+7]-Vz_now[i][j-8]))-
                        2.0*yita*V[i][j]*V[i][j]*Vv[i][j]*Vv[i][j]*1.0/dx*(
                                C1*(Cita_now[i][j]-Cita_now[i][j-1])+
                                C2*(Cita_now[i][j+1]-Cita_now[i][j-2])+
                                C3*(Cita_now[i][j+2]-Cita_now[i][j-3])+
                                C4*(Cita_now[i][j+3]-Cita_now[i][j-4])+
                                C5*(Cita_now[i][j+4]-Cita_now[i][j-5])+
                                C6*(Cita_now[i][j+5]-Cita_now[i][j-6])+
                                C7*(Cita_now[i][j+6]-Cita_now[i][j-7])+
                                C8*(Cita_now[i][j+7]-Cita_now[i][j-8]));


                    if(i==nx_location && j==nz_location)
                    {
                        temp=temp+wavelet;
                    }

                    P_aft[i][j]=P_now[i][j]+temp*dt;
                }
            }


            //*时间片波场值替换*//
            for(i=0; i<Nx; i++)
            {
                for(j=0; j<Nz; j++)
                {
                    Vz_pre[i][j]=Vz_now[i][j];
                    Vx_pre[i][j]=Vx_now[i][j];
                    P_pre[i][j]=P_now[i][j];
                    P_now[i][j]=P_aft[i][j];
                    Psa_now[i][j]=Psa_aft[i][j];
                    Cita_now[i][j]=Cita_aft[i][j];
                }
            }


            //----------------------------------------------//
            //            应用吸收边界条件                  //
            //----------------------------------------------//

            //*左边界*//
            for(i=PML;i<Nx-PML;i++)
                for(j=0;j<PML;j++)
                {
                    Vx_pre[i][j]=Vx_pre[i][j]*absorb[i][j];
                    Vx_now[i][j]=Vx_now[i][j]*absorb[i][j];
                    Vx_aft[i][j]=Vx_aft[i][j]*absorb[i][j];
                    Vz_pre[i][j]=Vz_pre[i][j]*absorb[i][j];
                    Vz_now[i][j]=Vz_now[i][j]*absorb[i][j];
                    Vz_aft[i][j]=Vz_aft[i][j]*absorb[i][j];
                    P_pre[i][j]=P_pre[i][j]*absorb[i][j];
                    P_now[i][j]=P_now[i][j]*absorb[i][j];
                    P_aft[i][j]=P_aft[i][j]*absorb[i][j];
                }
            //*右边界*//
            for(i=PML;i<Nx-PML;i++)
                for(j=Nz-PML;j<Nz;j++)
                {
                    Vx_pre[i][j]=Vx_pre[i][j]*absorb[i][j];
                    Vx_now[i][j]=Vx_now[i][j]*absorb[i][j];
                    Vx_aft[i][j]=Vx_aft[i][j]*absorb[i][j];
                    Vz_pre[i][j]=Vz_pre[i][j]*absorb[i][j];
                    Vz_now[i][j]=Vz_now[i][j]*absorb[i][j];
                    Vz_aft[i][j]=Vz_aft[i][j]*absorb[i][j];
                    P_pre[i][j]=P_pre[i][j]*absorb[i][j];
                    P_now[i][j]=P_now[i][j]*absorb[i][j];
                    P_aft[i][j]=P_aft[i][j]*absorb[i][j];
                }
            //*上边界*//
            for(i=0;i<PML;i++)
                for(j=0;j<Nz;j++)
                {
                    Vx_pre[i][j]=Vx_pre[i][j]*absorb[i][j];
                    Vx_now[i][j]=Vx_now[i][j]*absorb[i][j];
                    Vx_aft[i][j]=Vx_aft[i][j]*absorb[i][j];
                    Vz_pre[i][j]=Vz_pre[i][j]*absorb[i][j];
                    Vz_now[i][j]=Vz_now[i][j]*absorb[i][j];
                    Vz_aft[i][j]=Vz_aft[i][j]*absorb[i][j];
                    P_pre[i][j]=P_pre[i][j]*absorb[i][j];
                    P_now[i][j]=P_now[i][j]*absorb[i][j];
                    P_aft[i][j]=P_aft[i][j]*absorb[i][j];
                }
            //*下边界*//
            for(i=Nx-PML;i<Nx;i++)
                for(j=0;j<Nx;j++)
                {
                    Vx_pre[i][j]=Vx_pre[i][j]*absorb[i][j];
                    Vx_now[i][j]=Vx_now[i][j]*absorb[i][j];
                    Vx_aft[i][j]=Vx_aft[i][j]*absorb[i][j];
                    Vz_pre[i][j]=Vz_pre[i][j]*absorb[i][j];
                    Vz_now[i][j]=Vz_now[i][j]*absorb[i][j];
                    Vz_aft[i][j]=Vz_aft[i][j]*absorb[i][j];
                    P_pre[i][j]=P_pre[i][j]*absorb[i][j];
                    P_now[i][j]=P_now[i][j]*absorb[i][j];
                    P_aft[i][j]=P_aft[i][j]*absorb[i][j];
                }


            //-------------------------------------------//
            //               地震记录record              //
            //-------------------------------------------//
            for(i=PML; i<Nx-PML; i++)
            {
                record[i][k]=P_now[i][nz_location];
            }

        }
        printf("\n");


        //*波场时间递推计算结束*//


        //-------------------------------------------//
        //       存储地表地震记录record              //
        //-------------------------------------------//
        IO_func(NX, NZ, PML, Nt, dx, dz, dt, nx_location-PML, nz_location-PML, record);
    }



    printf("\n计算结束啦\n");
    //---------------------------------------------//
    //                 释放内存                    //
    //---------------------------------------------//
    for(i=0;i<Nx;i++)
        delete []	V[i];
    delete []V;
    for(i=0;i<Nx;i++)
        delete []	Vv[i];
    delete []Vv;
    for(i=0;i<Nx;i++)
        delete []	DEN[i];
    delete []DEN;
    for(i=0;i<Nx;i++)
        delete []	record[i];
    delete []record;

    for(i=0;i<Nx;i++)
        delete []	P_now[i];
    delete []P_now;
    for(i=0;i<Nx;i++)
        delete []	P_pre[i];
    delete []P_pre;
    for(i=0;i<Nx;i++)
        delete []	P_aft[i];
    delete []P_aft;

    for(i=0;i<Nx;i++)
        delete []	Vx_now[i];
    delete []Vx_now;
    for(i=0;i<Nx;i++)
        delete []	Vx_pre[i];
    delete []Vx_pre;
    for(i=0;i<Nx;i++)
        delete []	Vx_aft[i];
    delete []Vx_aft;

    for(i=0;i<Nx;i++)
        delete []	Vz_now[i];
    delete []Vz_now;
    for(i=0;i<Nx;i++)
        delete []	Vz_aft[i];
    delete []Vz_aft;
    for(i=0;i<Nx;i++)
        delete []	Vz_pre[i];
    delete []Vz_pre;

    for(i=0;i<Nx;i++)
        delete []	K_now[i];
    delete []K_now;
    for(i=0;i<Nx;i++)
        delete []	K_aft[i];
    delete []K_aft;
    for(i=0;i<Nx;i++)
        delete []	K_pre[i];
    delete []K_pre;

    for(i=0;i<Nx;i++)
        delete []	Cita_now[i];
    delete []Cita_now;
    for(i=0;i<Nx;i++)
        delete []	Cita_aft[i];
    delete []Cita_aft;
    for(i=0;i<Nx;i++)
        delete []	Cita_pre[i];
    delete []Cita_pre;

    for(i=0;i<Nx;i++)
        delete []	Psa_now[i];
    delete []Psa_now;
    for(i=0;i<Nx;i++)
        delete []	Psa_aft[i];
    delete []Psa_aft;
    for(i=0;i<Nx;i++)
        delete []	Psa_pre[i];
    delete []Psa_pre;

    return 0;
}


