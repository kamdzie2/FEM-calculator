#ifndef STRUKTURY_H
#define STRUKTURY_H

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <math.h>
#include <vector>
using namespace std;

const int N=4;

struct SOE
{
    double **mglobH;
    double *wglobP;
    double **mglobC;
    double *t0;
};

struct Global_data
{
    int SimulationTime;
    int SimulationStepTime;
    int Conductivity;
    int Alfa;
    int Tot;
    int InitialTemp;
    int Density;
    int SpecificHeat;
    int NodesNumber;
    int ElementsNumber;
};

struct Node
{
    double x;
    double y;
    double BC=0;
};

struct Element
{
    int ID[N];
    double H[N][N];
    double Hbc[N][N];
    double P[N];
    double HFULL[N][N];
    double C[N][N];
};

struct Grid
{
    Element *e;
    Node *n;

};

struct Nwezlow
{
    int n;
    int *k;
    double *wx;
    double *wa;
};

struct ElementUniwersalny
{
    double **dKsi;
    double **dEta;
    double **dx;
    double **dy;
    double **Nkszt;
};

struct TabXY
{

    double tx[N];
    double ty[N];
};

struct Jakobian
{
    double gl;   // jakobian[0][0]
    double gp;   // jakobian[0][1]
    double dl;   // jakobian[1][0]
    double dp;   // jakobian[1][1]
};

struct GaussIntegration
{
    static const int Nn=4;
    double w;
    double pc[Nn][Nn];
    double p[Nn];
};



#endif // STRUKTURY_H
