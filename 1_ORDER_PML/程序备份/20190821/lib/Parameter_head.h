#pragma once
typedef struct {

    //**模型参数
    int NX;//x方向样点数
    int NZ;//z方向样点数
    int PML;//PML边界厚度
    int Nt;//模拟时间样点数
    float dx;//x方向取样
    float dz;//z方向取样
    float dt;//时间取样

    //**各向异性参数**Thomson表征形式
    float epsilon;
    float delt;
    float yita;

    //**有限差分系数
    float C1; 
    float C2; 
    float C3; 
    float C4;
    float C5; 
    float C6; 
    float C7; 
    float C8;//差分系数

    //计算空间大小
    int Nx;//x方向采样点数(x-z方向各留了PML行用于边界处理)
    int Nz;

    //子波参数设置
    int Lw;   //子波长度
    float freq; //主频

    //**观测系统
    int nx_location;
    int nz_location;
    int Ns;//炮数

}PARAMETER;
