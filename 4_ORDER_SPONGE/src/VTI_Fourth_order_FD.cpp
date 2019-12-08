#include<iostream>
#include<cmath>

using namespace std;

#define PI 3.1415926
#define nx 400
#define nz 400
#define nt 1400
#define dx 10
#define dz 10
#define dt 0.0005
#define npmlx 200
#define npmlz 200

//子函数
void **alloc2 (int n1, int n2, int size)
{
    int i2;
    void **p;

    if ((p=(void**)malloc(n2*sizeof(void*)))==NULL)                    //申请数组p[n2]，数组内都是指针，由p指针指向p(0)首地址
        return NULL;
    if ((p[0]=(void*)malloc(n2*n1*size))==NULL)
    {                                                                  //p(0)为指针指向一个n2*n1*size大小的段的首地址
        free(p);
        return NULL;
    }
    for (i2=0; i2<n2; i2++)
        p[i2] = (char*)p[0]+size*n1*i2;                                //p(1),p(2)...指向1列，2列...的首地址
    return p;
}


void exchange23(float **a1,float **a2,float **a3)

{
    int i,j;


    for (i=0;i<nz+npmlz*2;i++)
    {
        for (j=0;j<nx+npmlx*2;j++)
        {

            a1[i][j]=a2[i][j];
            a2[i][j]=a3[i][j];


        }
    }

}




//主函数
int main()
{
    FILE *fp;
    int nstartx,nstartz,nendx,nendz,npmlsx,npmlex,npmlsz,npmlez;
    int i,j,l,k,nx_location,nz_location;

    float d1,d2,d3,d4,fre,ft[nt],t,w0,w1,w2,w3,w4,w5,t1,par,nzall,nxall;
    float **location,**f_pre,**f_now,**f_next,**p_pre,**p_now,**p_next;
    //float **time_left;
    float **dxx,**dzz,**a,**b,**vp,**v;
    float **axx,**azz,**ax,**zx,**bx,**bz;
    float record[nt][nx+npmlx*2];

    float temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,temp9,temp10;
    float a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,b1,b2,b3,b4,b5,b6,c,d,g1,g2,g3,e,f;
    a1=0.32673;
    a2=0.097676;
    a3=0.018563;
    a4=0.0072016;
    a5=-0.0017692;
    a6=0.0016338;
    a7=-0.000002584;
    a8=0.000020026;
    a9=-0.000184;
    a10=0.000090819;
    b1=0.447;
    b2=0.12637;
    b3=0.02139;
    b4=-0.0048025;
    b5=0.006032;
    b6=-0.0023837;
    c=0.69263;
    d=0.45944;
    g1=1.1273;
    g2=-0.087043;
    g3=-0.04665;
    e=0.49089;
    f=0.38731;
    //各向异性参数

    float epsilon,delt,yita,lamda;
    epsilon=0.2;
    delt=0.2;
    lamda=1+2.0*epsilon;
    yita=epsilon-delt;



    //差分边界
    nzall=nz+npmlz*2;    //总的计算大小 [nzall,nxall]
    nxall=nx+npmlx*2;
    nstartx=5;           //差分计算开始的地方 nstartx nstartz
    nstartz=5;
    nendx=nxall-nstartx;//差分计算结束的地方 nendx nendz
    nendz=nzall-nstartz;
    npmlsx=npmlx;        //PML边界
    npmlsz=npmlz;
    npmlex=nx+npmlsx;
    npmlez=nz+npmlsz;
    par=0.004;

    //二阶导数的10阶有限差分系数
    w0=-5.8544444444;
    w1=3.3333333333;
    w2=-0.4761904762;
    w3=0.0793650794;
    w4=-0.0099206349;
    w5=0.0006349206;

    //生成二维动态数组
    location=(float **)alloc2(nzall,nxall,sizeof(float));

    p_pre=(float **)alloc2(nzall,nxall,sizeof(float));
    p_now=(float **)alloc2(nzall,nxall,sizeof(float));
    p_next=(float **)alloc2(nzall,nxall,sizeof(float));

    f_pre=(float **)alloc2(nzall,nxall,sizeof(float));
    f_now=(float **)alloc2(nzall,nxall,sizeof(float));
    f_next=(float **)alloc2(nzall,nxall,sizeof(float));

    dxx=(float **)alloc2(nzall,nxall,sizeof(float));
    dzz=(float **)alloc2(nzall,nxall,sizeof(float));

    vp=(float **)alloc2(nzall,nxall,sizeof(float));
    //v=(float **)alloc2(1601,601,s//izeof(float));


    //震源设置
    for (i=0;i<nzall;i++)
    {
        for(j=0;j<nxall;j++)
        {
            location[i][j]=0;
            f_pre[i][j]=0;
            f_now[i][j]=0;
            f_next[i][j]=0;

            p_pre[i][j]=0;
            p_now[i][j]=0;
            p_next[i][j]=0;

            dxx[i][j]=0;
            dzz[i][j]=0;

            vp[i][j]=3000;

        }
    }




/*
       for (i=nz/2+npmlz;i<nzall;i++)
       {

       for(j=0;j<nxall;j++)
       {
       vp[i][j]=5000;
       }
       }

*/
/*
    FILE *fpvp;
    fpvp=fopen("../file/sigsbee_1601x601.dat","rb");
    for(i=0;i<1601;i++)
    {
        fread(v[i],sizeof(float),601,fpvp);
    }
    fclose(fpvp);
    for (i=0;i<nxall;i++)
    {

        for(j=0;j<nzall;j++)
        {

            vp[i][j]=v[i][j];
        }
    }

*/

       if((fp = fopen ("../file/vp.dat", "wb"))!=NULL)
    {

        for (j=npmlsx;j<npmlex;j++)
        {
            for (i=npmlsz;i<npmlez;i++)
            {
                fwrite (&vp[i][j] , sizeof(float), 1, fp);

            }
        }
        fclose (fp);
    }
    cout<<"Finsh vp save"<<endl;






    for (i=0;i<nt;i++)
    {
        for(j=0;j<nx;j++)
        {
            record[i][j]=0;
        }
    }



    //震源加载设置
    nz_location=nzall/2;
    nx_location=nxall/2;
    location[nz_location][nx_location]=1;

	/*
    for(i=nz_location;i<=nz_location+4;i++){
        for(j=nx_location;j<=nx_location+4;j++){
            location[i][j]=1.0/(sqrt((i-nz_location)*(i-nz_location)*dz*dz+(j-nx_location)*(j-nx_location)*dx*dx)+1);
            cout<<location[i][j]<<endl;
        }
    }
	*/




    //雷克子波，f为频率,ft(t)为雷克子波

    fre=28;
    for(i=0;i<nt;i++)
    {
        t=dt*i;
        t1=1/fre;//双边雷克子波
        //t1=0;//单边雷克子波
        ft[i]=(1-2*PI*PI*fre*fre*(t-t1)*(t-t1))*exp(-PI*PI*fre*fre*(t-t1)*(t-t1));
    }





    //PML
    for (i=0;i<nzall;i++)
    {
        for (j=0;j<nxall;j++)
        {
            if (i<npmlsz)
            {

                dxx[i][j]=par*(npmlsz-i)*(npmlsz-i);

            }

            if (i>=npmlez)
            {
                dxx[i][j]=par*(i-npmlez)*(i-npmlez);

            }
        }


    }

    for (i=0;i<nzall;i++)
    {
        for (j=0;j<nxall;j++)
        {
            if (j<npmlsx)
            {

                dzz[i][j]=par*(npmlsx-j)*(npmlsx-j);

            }

            if (j>=npmlez)
            {
                dzz[i][j]=par*(j-npmlex)*(j-npmlex);

            }
        }


    }




    //差分的计算



    for(l=0;l<nt;l++)
    {
        if(l%100==0)
        {
            printf("l=%d\n",l);
        }


        for (i=nstartz;i<nendz;i++)
        {
            for (j=nstartx;j<nendx;j++)
            {
                //p对xx方向的差分
                /*
                   temp1=b1*(c*(p_now[i+1][j]-2*p_now[i][j]+p_now[i-1][j])+d*1.0/4*(p_now[i+2][j]-2*p_now[i][j]+p_now[i-2][j]));
                   temp2=b2*(c*(p_now[i+1][j]-2*p_now[i][j]+p_now[i-1][j])+d*1.0/4*(p_now[i+2][j]-2*p_now[i][j]+p_now[i-2][j]));
                   temp2=temp2+b2*(c*(p_now[i+1][j]-2*p_now[i][j]+p_now[i-1][j])+d*1.0/4*(p_now[i+2][j]-2*p_now[i][j]+p_now[i-2][j]));
                   temp2=temp2+b2*(c*(p_now[i+1][j+1]-2*p_now[i][j+1]+p_now[i-1][j+1])+d*1.0/4*(p_now[i+2][j+1]-2*p_now[i][j+1]+p_now[i-2][j+1]));
                   temp2=temp2+b2*(c*(p_now[i+1][j-1]-2*p_now[i][j-1]+p_now[i-1][j-1])+d*1.0/4*(p_now[i+2][j-1]-2*p_now[i][j-1]+p_now[i-2][j-1]));
                   temp3=b3*(c*(p_now[i+1][j-1]-2*p_now[i][j-1]+p_now[i-1][j-1])+d*1.0/4*(p_now[i+2][j-1]-2*p_now[i][j-1]+p_now[i-2][j-1]));
                   temp3=temp3+b3*(c*(p_now[i+1][j+1]-2*p_now[i][j+1]+p_now[i-1][j+1])+d*1.0/4*(p_now[i+2][j+1]-2*p_now[i][j+1]+p_now[i-2][j+1]));
                   temp3=temp3+b3*(c*(p_now[i+1][j-1]-2*p_now[i][j-1]+p_now[i-1][j-1])+d*1.0/4*(p_now[i+2][j-1]-2*p_now[i][j-1]+p_now[i-2][j-1]));
                   temp3=temp3+b3*(c*(p_now[i+1][j+1]-2*p_now[i][j+1]+p_now[i-1][j+1])+d*1.0/4*(p_now[i+2][j+1]-2*p_now[i][j+1]+p_now[i-2][j+1]));
                   temp4=b4*(c*(p_now[i+1][j+2]-2*p_now[i][j+2]+p_now[i-1][j+2])+d*1.0/4*(p_now[i+2][j+2]-2*p_now[i][j+2]+p_now[i-2][j+2]));
                   temp4=temp4+b4*(c*(p_now[i+1][j]-2*p_now[i][j]+p_now[i-1][j])+d*1.0/4*(p_now[i+2][j]-2*p_now[i][j]+p_now[i-2][j]));
                   temp4=temp4+b4*(c*(p_now[i+1][j-2]-2*p_now[i][j-2]+p_now[i-1][j-2])+d*1.0/4*(p_now[i+2][j-2]-2*p_now[i][j-2]+p_now[i-2][j-2]));
                   temp4=temp4+b4*(c*(p_now[i+1][j]-2*p_now[i][j]+p_now[i-1][j])+d*1.0/4*(p_now[i+2][j]-2*p_now[i][j]+p_now[i-2][j]));
                   temp5=b5*(c*(p_now[i+1][j+2]-2*p_now[i][j+2]+p_now[i-1][j+2])+d*1.0/4*(p_now[i+2][j+2]-2*p_now[i][j+2]+p_now[i-2][j+2]));
                   temp5=temp5+b5*(c*(p_now[i+1][j+2]-2*p_now[i][j+2]+p_now[i-1][j+2])+d*1.0/4*(p_now[i+2][j+2]-2*p_now[i][j+2]+p_now[i-2][j+2]));
                   temp5=temp5+b5*(c*(p_now[i+1][j+1]-2*p_now[i][j+1]+p_now[i-1][j+1])+d*1.0/4*(p_now[i+2][j+1]-2*p_now[i][j+1]+p_now[i-2][j+1]));
                   temp5=temp5+b5*(c*(p_now[i+1][j+1]-2*p_now[i][j+1]+p_now[i-1][j+1])+d*1.0/4*(p_now[i+2][j+1]-2*p_now[i][j+1]+p_now[i-2][j+1]));
                   temp5=temp5+b5*(c*(p_now[i+1][j-2]-2*p_now[i][j-2]+p_now[i-1][j-2])+d*1.0/4*(p_now[i+2][j-2]-2*p_now[i][j-2]+p_now[i-2][j-2]));
                   temp5=temp5+b5*(c*(p_now[i+1][j-2]-2*p_now[i][j-2]+p_now[i-1][j-2])+d*1.0/4*(p_now[i+2][j-2]-2*p_now[i][j-2]+p_now[i-2][j-2]));
                   temp5=temp5+b5*(c*(p_now[i+1][j-1]-2*p_now[i][j-1]+p_now[i-1][j-1])+d*1.0/4*(p_now[i+2][j-1]-2*p_now[i][j-1]+p_now[i-2][j-1]));
                   temp5=temp5+b5*(c*(p_now[i+1][j-1]-2*p_now[i][j-1]+p_now[i-1][j-1])+d*1.0/4*(p_now[i+2][j-1]-2*p_now[i][j-1]+p_now[i-2][j-1]));
                   temp6=b6*(c*(p_now[i+1][j-2]-2*p_now[i][j-2]+p_now[i-1][j-2])+d*1.0/4*(p_now[i+2][j-2]-2*p_now[i][j-2]+p_now[i-2][j-2]));
                   temp6=temp6+b6*(c*(p_now[i+1][j-2]-2*p_now[i][j-2]+p_now[i-1][j-2])+d*1.0/4*(p_now[i+2][j-2]-2*p_now[i][j-2]+p_now[i-2][j-2]));
                   temp6=temp6+b6*(c*(p_now[i+1][j+2]-2*p_now[i][j+2]+p_now[i-1][j+2])+d*1.0/4*(p_now[i+2][j+2]-2*p_now[i][j+2]+p_now[i-2][j+2]));
                   temp6=temp6+b6*(c*(p_now[i+1][j+2]-2*p_now[i][j+2]+p_now[i-1][j+2])+d*1.0/4*(p_now[i+2][j+2]-2*p_now[i][j+2]+p_now[i-2][j+2]));
                   d1=1.0/pow(dx,2)*(temp1+temp2+temp3+temp4+temp5+temp6);
                   */
                d1=1.0/(dx*dx)*(w0*p_now[i][j]+w1*(p_now[i+1][j]+p_now[i-1][j])+w2*(p_now[i+2][j]+p_now[i-2][j])+w3*(p_now[i+3][j]+p_now[i-3][j])+w4*(p_now[i+4][j]+p_now[i-4][j])+w5*(p_now[i+5][j]+p_now[i-5][j]));

                //p对yy方向的差分
                /*
                   temp1=b1*(c*(p_now[i][j+1]-2*p_now[i][j]+p_now[i][j-1])+d*1.0/4*(p_now[i][j+2]-2*p_now[i][j]+p_now[i][j-2]));
                   temp2=b2*(c*(p_now[i][j+1]-2*p_now[i][j]+p_now[i][j-1])+d*1.0/4*(p_now[i][j+2]-2*p_now[i][j]+p_now[i][j-2]));
                   temp2=temp2+b2*(c*(p_now[i][j+1]-2*p_now[i][j]+p_now[i][j-1])+d*1.0/4*(p_now[i][j+2]-2*p_now[i][j]+p_now[i][j-2]));
                   temp2=temp2+b2*(c*(p_now[i+1][j+1]-2*p_now[i+1][j]+p_now[i+1][j-1])+d*1.0/4*(p_now[i+1][j+2]-2*p_now[i+1][j]+p_now[i+1][j-2]));
                   temp2=temp2+b2*(c*(p_now[i-1][j+1]-2*p_now[i-1][j]+p_now[i-1][j-1])+d*1.0/4*(p_now[i-1][j+2]-2*p_now[i-1][j]+p_now[i-1][j-2]));
                   temp3=b3*(c*(p_now[i-1][j+1]-2*p_now[i-1][j]+p_now[i-1][j-1])+d*1.0/4*(p_now[i-1][j+2]-2*p_now[i-1][j]+p_now[i-1][j-2]));
                   temp3=temp3+b3*(c*(p_now[i+1][j+1]-2*p_now[i+1][j]+p_now[i+1][j-1])+d*1.0/4*(p_now[i+1][j+2]-2*p_now[i+1][j]+p_now[i+1][j-2]));
                   temp3=temp3+b3*(c*(p_now[i-1][j+1]-2*p_now[i-1][j]+p_now[i-1][j-1])+d*1.0/4*(p_now[i-1][j+2]-2*p_now[i-1][j]+p_now[i-1][j-2]));
                   temp3=temp3+b3*(c*(p_now[i+1][j+1]-2*p_now[i+1][j]+p_now[i+1][j-1])+d*1.0/4*(p_now[i+1][j+2]-2*p_now[i+1][j]+p_now[i+1][j-2]));
                   temp4=b4*(c*(p_now[i+2][j+1]-2*p_now[i+2][j]+p_now[i+2][j-1])+d*1.0/4*(p_now[i+2][j+2]-2*p_now[i+2][j]+p_now[i+2][j-2]));
                   temp4=temp4+b4*(c*(p_now[i][j+1]-2*p_now[i][j]+p_now[i][j-1])+d*1.0/4*(p_now[i][j+2]-2*p_now[i][j]+p_now[i][j-2]));
                   temp4=temp4+b4*(c*(p_now[i-2][j+1]-2*p_now[i-2][j]+p_now[i-2][j-1])+d*1.0/4*(p_now[i-2][j+2]-2*p_now[i-2][j]+p_now[i-2][j-2]));
                   temp4=temp4+b4*(c*(p_now[i][j+1]-2*p_now[i][j]+p_now[i][j-1])+d*1.0/4*(p_now[i][j+2]-2*p_now[i][j]+p_now[i][j-2]));
                   temp5=b5*(c*(p_now[i+2][j+1]-2*p_now[i+2][j]+p_now[i+2][j-1])+d*1.0/4*(p_now[i+2][j+2]-2*p_now[i+2][j]+p_now[i+2][j-2]));
                   temp5=temp5+b5*(c*(p_now[i+2][j+1]-2*p_now[i+2][j]+p_now[i+2][j-1])+d*1.0/4*(p_now[i+2][j+2]-2*p_now[i+2][j]+p_now[i+2][j-2]));
                   temp5=temp5+b5*(c*(p_now[i+1][j+1]-2*p_now[i+1][j]+p_now[i+1][j-1])+d*1.0/4*(p_now[i+1][j+2]-2*p_now[i+1][j]+p_now[i+1][j-2]));
                   temp5=temp5+b5*(c*(p_now[i+1][j+1]-2*p_now[i+1][j]+p_now[i+1][j-1])+d*1.0/4*(p_now[i+1][j+2]-2*p_now[i+1][j]+p_now[i+1][j-2]));
                   temp5=temp5+b5*(c*(p_now[i-2][j+1]-2*p_now[i-2][j]+p_now[i-2][j-1])+d*1.0/4*(p_now[i-2][j+2]-2*p_now[i-2][j]+p_now[i-2][j-2]));
                   temp5=temp5+b5*(c*(p_now[i-2][j+1]-2*p_now[i-2][j]+p_now[i-2][j-1])+d*1.0/4*(p_now[i-2][j+2]-2*p_now[i-2][j]+p_now[i-2][j-2]));
                   temp5=temp5+b5*(c*(p_now[i-1][j+1]-2*p_now[i-1][j]+p_now[i-1][j-1])+d*1.0/4*(p_now[i-1][j+2]-2*p_now[i-1][j]+p_now[i-1][j-2]));
                   temp5=temp5+b5*(c*(p_now[i-1][j+1]-2*p_now[i-1][j]+p_now[i-1][j-1])+d*1.0/4*(p_now[i-1][j+2]-2*p_now[i-1][j]+p_now[i-1][j-2]));
                   temp6=b6*(c*(p_now[i-2][j+1]-2*p_now[i-2][j]+p_now[i-2][j-1])+d*1.0/4*(p_now[i-2][j+2]-2*p_now[i-2][j]+p_now[i-2][j-2]));
                   temp6=temp6+b6*(c*(p_now[i-2][j+1]-2*p_now[i-2][j]+p_now[i-2][j-1])+d*1.0/4*(p_now[i-2][j+2]-2*p_now[i-2][j]+p_now[i-2][j-2]));
                   temp6=temp6+b6*(c*(p_now[i+2][j+1]-2*p_now[i+2][j]+p_now[i+2][j-1])+d*1.0/4*(p_now[i+2][j+2]-2*p_now[i+2][j]+p_now[i+2][j-2]));
                   temp6=temp6+b6*(c*(p_now[i+2][j+1]-2*p_now[i+2][j]+p_now[i+2][j-1])+d*1.0/4*(p_now[i+2][j+2]-2*p_now[i+2][j]+p_now[i+2][j-2]));
                   d2=1.0/pow(dz,2)*(temp1+temp2+temp3+temp4+temp5+temp6);
                   */
                d2=1.0/(dz*dz)*(w0*p_now[i][j]+w1*(p_now[i][j+1]+p_now[i][j-1])+w2*(p_now[i][j+2]+p_now[i][j-2])+w3*(p_now[i][j+3]+p_now[i][j-3])+w4*(p_now[i][j+4]+p_now[i][j-4])+w5*(p_now[i][j+5]+p_now[i][j-5]));

                //F对xy方向的差分
                /*
                   temp1=e*1.0/(pow(dx,2)*pow(dz,2))*(f_now[i+1][j]-2.0*f_now[i][j]+f_now[i-1][j]-2.0*f_now[i+1][j]+4.0*f_now[i][j]-2.0*f_now[i-1][j]+f_now[i+1][j]-2.0*f_now[i][j]+f_now[i-1][j]);
                   temp1=temp1+f*1.0/(16.0*pow(dx,2)*pow(dz,2))*(f_now[i+2][j]-2.0*f_now[i][j]+f_now[i-2][j]-2.0*f_now[i+2][j]+4.0*f_now[i][j]-2.0*f_now[i-2][j]+f_now[i+2][j]-2.0*f_now[i][j]+f_now[i-2][j]);
                   temp2=e*1.0/(pow(dx,2)*pow(dz,2))*(f_now[i+1][j-1]-2.0*f_now[i][j-1]+f_now[i-1][j-1]-2.0*f_now[i+1][j-1]+4.0*f_now[i][j-1]-2.0*f_now[i-1][j-1]+f_now[i+1][j-1]-2.0*f_now[i][j-1]+f_now[i-1][j-1]);
                   temp2=temp2+f*1.0/(16.0*pow(dx,2)*pow(dz,2))*(f_now[i+2][j-1]-2.0*f_now[i][j-1]+f_now[i-2][j-1]-2.0*f_now[i+2][j-1]+4.0*f_now[i][j-1]-2.0*f_now[i-2][j-1]+f_now[i+2][j-1]-2.0*f_now[i][j-1]+f_now[i-2][j-1]);
                   temp3=e*1.0/(pow(dx,2)*pow(dz,2))*(f_now[i+1][j+1]-2.0*f_now[i][j+1]+f_now[i-1][j+1]-2.0*f_now[i+1][j+1]+4.0*f_now[i][j+1]-2.0*f_now[i-1][j+1]+f_now[i+1][j+1]-2.0*f_now[i][j+1]+f_now[i-1][j+1]);
                   temp3=temp3+f*1.0/(16.0*pow(dx,2)*pow(dz,2))*(f_now[i+2][j+1]-2.0*f_now[i][j+1]+f_now[i-2][j+1]-2.0*f_now[i+2][j+1]+4.0*f_now[i][j+1]-2.0*f_now[i-2][j+1]+f_now[i+2][j+1]-2.0*f_now[i][j+1]+f_now[i-2][j+1]);
                   temp4=e*1.0/(pow(dx,2)*pow(dz,2))*(f_now[i+1][j-2]-2.0*f_now[i][j-2]+f_now[i-1][j-2]-2.0*f_now[i+1][j-2]+4.0*f_now[i][j-2]-2.0*f_now[i-1][j-2]+f_now[i+1][j-2]-2.0*f_now[i][j-2]+f_now[i-1][j-2]);
                   temp4=temp4+f*1.0/(16.0*pow(dx,2)*pow(dz,2))*(f_now[i+2][j-2]-2.0*f_now[i][j-2]+f_now[i-2][j-2]-2.0*f_now[i+2][j-2]+4.0*f_now[i][j-2]-2.0*f_now[i-2][j-2]+f_now[i+2][j-2]-2.0*f_now[i][j-2]+f_now[i-2][j-2]);
                   temp5=e*1.0/(pow(dx,2)*pow(dz,2))*(f_now[i+1][j+2]-2.0*f_now[i][j+2]+f_now[i-1][j+2]-2.0*f_now[i+1][j+2]+4.0*f_now[i][j+2]-2.0*f_now[i-1][j+2]+f_now[i+1][j+2]-2.0*f_now[i][j+2]+f_now[i-1][j+2]);
                   temp5=temp5+f*1.0/(16.0*pow(dx,2)*pow(dz,2))*(f_now[i+2][j+2]-2.0*f_now[i][j+2]+f_now[i-2][j+2]-2.0*f_now[i+2][j+2]+4.0*f_now[i][j+2]-2.0*f_now[i-2][j+2]+f_now[i+2][j+2]-2.0*f_now[i][j+2]+f_now[i-2][j+2]);
                   d3=g1*temp1+g2*temp2+g2*temp3+g3*temp4+g3*temp5;
                   */
                d3=1.0/(dx*dx*dz*dz)*((f_now[i+1][j+1]-2*f_now[i][j+1]+f_now[i-1][j+1])-2*(f_now[i+1][j]-2*f_now[i][j]+f_now[i-1][j])+(f_now[i+1][j-1]-2*f_now[i][j-1]+f_now[i-1][j-1]));

                //对p的时间差分运算
                //if(l==0){
				d4=lamda*pow(vp[i][j],2)*d1+pow(vp[i][j],2)*d2-2.0*yita*pow(vp[i][j],4)*d3+location[i][j]*ft[l];
				//}
               //else {d4=lamda*pow(vp[i][j],2)*d1+pow(vp[i][j],2)*d2-2.0*yita*pow(vp[i][j],4)*d3+ft[l];}

                p_next[i][j]=2.0*p_now[i][j]-p_pre[i][j]+dt*dt*d4;
                f_next[i][j]=2.0*f_now[i][j]-f_pre[i][j]+dt*dt*p_now[i][j];
                /******************计算F***************/
                /*
                   temp1=f_now[i][j];
                   temp2=f_now[i-1][j]+f_now[i+1][j]+f_now[i][j+1]+f_now[i][j-1]+f_now[i][j]+f_now[i][j];
                   temp3=f_now[i+1][j]+f_now[i+1][j]+f_now[i-1][j]+f_now[i-1][j]+f_now[i+1][j+1]+f_now[i+1][j-1]+f_now[i-1][j+1]+f_now[i-1][j-1]+f_now[i][j+1]+f_now[i][j-1]+f_now[i][j+1]+f_now[i][j-1];
                   temp4=f_now[i+1][j+1]+f_now[i+1][j-1]+f_now[i+1][j+1]+f_now[i-1][j+1]+f_now[i-1][j-1]+f_now[i-1][j+1]+f_now[i+1][j-1]+f_now[i-1][j-1];
                   temp5=f_now[i+2][j]+f_now[i-2][j]+f_now[i][j]+f_now[i][j]+f_now[i][j+2]+f_now[i][j-2];
                   temp6=f_now[i+2][j]+f_now[i+2][j]+f_now[i-2][j]+f_now[i-2][j]+f_now[i+2][j+1]+f_now[i+2][j-1]+f_now[i-2][j+1]+f_now[i-2][j-1];
                   temp6=temp6+f_now[i+1][j]+f_now[i+1][j]+f_now[i-1][j]+f_now[i-1][j]+f_now[i+1][j+2]+f_now[i+1][j-2]+f_now[i-1][j+2]+f_now[i-1][j-2];
                   temp6=temp6+f_now[i][j+1]+f_now[i][j-1]+f_now[i][j+1]+f_now[i][j-1]+f_now[i][j+2]+f_now[i][j-2]+f_now[i][j+2]+f_now[i][j-2];
                   temp7=f_now[i+1][j+1]+f_now[i+1][j-1]+f_now[i+1][j+1]+f_now[i+1][j-1]+f_now[i-1][j+1]+f_now[i-1][j-1]+f_now[i-1][j+1]+f_now[i-1][j-1]+f_now[i+1][j+2]+f_now[i+1][j-2]+f_now[i+1][j+2]+f_now[i+1][j-2];
                   temp7=temp7+f_now[i-1][j+2]+f_now[i-1][j-2]+f_now[i-1][j+2]+f_now[i-1][j-2]+f_now[i+2][j+1]+f_now[i+2][j-1]+f_now[i+2][j+1]+f_now[i+2][j-1]+f_now[i-2][j+1]+f_now[i-2][j-1]+f_now[i-2][j+1]+f_now[i-2][j-1];
                   temp8=f_now[i+2][j]+f_now[i+2][j]+f_now[i+2][j+2]+f_now[i+2][j-2]+f_now[i-2][j]+f_now[i-2][j]+f_now[i-2][j+2]+f_now[i-1][j-2]+f_now[i][j-2]+f_now[i][j+2]+f_now[i][j+2]+f_now[i][j-2];
                   temp9=f_now[i+2][j+2]+f_now[i+2][j-2]+f_now[i+2][j+2]+f_now[i+2][j-2]+f_now[i-2][j+2]+f_now[i-2][j-2]+f_now[i-2][j+2]+f_now[i-2][j-2]+f_now[i+1][j+2]+f_now[i+1][j-2]+f_now[i+1][j+2]+f_now[i+1][j-2];
                   temp9=temp9+f_now[i-1][j+2]+f_now[i-1][j-2]+f_now[i-1][j+2]+f_now[i-1][j-2]+f_now[i+2][j+1]+f_now[i+2][j-1]+f_now[i+2][j+1]+f_now[i+2][j-1]+f_now[i-2][j+1]+f_now[i-2][j-1]+f_now[i-2][j+1]+f_now[i-2][j-1];
                   temp10=f_now[i+2][j+2]+f_now[i+2][j-2]+f_now[i+2][j+2]+f_now[i+2][j-2]+f_now[i-2][j+2]+f_now[i-2][j-2]+f_now[i-2][j+2]+f_now[i-2][j-2];
                   f_now[i][j]=a1*temp1+a2*temp2+a3*temp3+a4*temp4+a5*temp5+a6*temp6+a7*temp7+a8*temp8+a9*temp9+a10*temp10;
                   */

            }
        }
        exchange23(f_pre,f_now,f_next);
        exchange23(p_pre,p_now,p_next);



        //时间为700ms时的模型数值
        if(l==700)
        {
            if((fp = fopen ("../file/aha.dat", "wb"))!=NULL)
            {

                for (j=npmlsx;j<npmlex;j++)
                {
                    for (i=npmlsz;i<npmlez;i++)
                    {
                        fwrite (&f_now[i][j] , sizeof(float), 1, fp);

                    }
                }
                fclose (fp);
            }
            cout<<"Finsh times save"<<endl;
        }




        if(l%20==0)
        {
            FILE *fpm;
            if((fpm = fopen ("../file/allseistime.dat", "a+"))!=NULL)
            {

                for (j=npmlsx;j<npmlex;j++)
                {
                    for (i=npmlsz;i<npmlez;i++)
                    {
                        fwrite (&f_now[i][j] , sizeof(float), 1, fpm);

                    }
                }
                fclose (fpm);
            }
        }


        //cout<<"Finsh times save"<<endl;



        /*地震记录*/

        for (j=npmlsx;j<npmlex;j++)
        {
            record[l][j]=f_now[nz_location][j];
        }




    }


    if((fp = fopen ("../file/seis.dat", "wb"))!=NULL)
    {


        for (j=npmlsx;j<npmlex;j++)
        {
            for (i=0;i<nt;i++)
            {
                fwrite (&record[i][j] , sizeof(float), 1, fp);

            }
        }
        fclose (fp);
        cout<<"Finsh seisdata save"<<endl;
    }







}

























