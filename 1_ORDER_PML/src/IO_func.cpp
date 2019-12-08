#include"../lib/IO_func.h"

double IO_func(struct PARAMETER *param, float **record)
{
    FILE *fpr;
    int i;

    //fpr=fopen("../file/record.su","a+");

    fpr=fopen("../file/record","a+");
    SUHEAD hdr;
    memset(&hdr, 0, 240);//SUHEAD全部置零

    float * ptr;
    for(i=param->PML;i<param->Nx-param->PML;i++)
    {
        hdr.dt     = int(1e-6*param->dt);
        hdr.ns     = param->Nt;
        hdr.fldr   = param->nx_location-param->PML;
        hdr.tracf  = i-param->PML;
        hdr.sx     = (param->nx_location-param->PML)*param->dx;
        hdr.gx     = (i-param->PML)*param->dx;
        hdr.offset = -(hdr.sx - hdr.gx);
        hdr.cdp    = (hdr.sx + hdr.gx) /2.0f / (param->dx*0.5f);
        //hdr.cdp    = 1;
        ptr = &(record[i][0]);
        fwrite( &hdr,           240,  1, fpr);
        fwrite( ptr , sizeof(float), param->Nt, fpr);
    }
    fclose(fpr);
}

