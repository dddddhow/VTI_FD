#include"../lib/IO_func.h"


double IO_func(int NX, int NZ, int PML, int Nt, float dx, float dz, float dt, int nx_location, int nz_location, float **record)
{
    FILE *fpr;
    int i;
    int Nx, Nz;
    Nx=NX+2*PML;
    Nz=NZ+2*PML;

    fpr=fopen("../file/record.su","a+");

    SUHEAD hdr;
    memset(&hdr, 0, 240);//SUHEAD全部置零

    float * ptr;
    for(i=PML;i<Nx-PML;i++)
    {
        hdr.dt     = int(1e6*0.0005);
        hdr.ns     = Nt;
        hdr.fldr   = nx_location;
        hdr.tracf  = i-PML;
        hdr.sx     = nx_location*dx;
        hdr.gx     = (i-PML)*dx;
        hdr.offset = (hdr.sx-hdr.gx)/2.0f/(1.0f/2.0f*dx);

        ptr = &(record[i][0]);
        fwrite( &hdr,           240,  1, fpr);
        fwrite( ptr , sizeof(float), Nt, fpr);
    }
    fclose(fpr);
    return 0;
}

