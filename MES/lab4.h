#ifndef LAB4_H
#define LAB4_H
#include "struktury.h"
#include "lab3.h"

class lab4: public lab3
{


    public:
    TabXY x;
    TabXY y;
    Element *macierze;
    ElementUniwersalny valXY;
    Jakobian **jk;
    double **jakob_wyn;
     double **jakob_odw;
    double **tabWEKx;
    double **tabWEKy;
    double **tabWEKsuma;

   void konstruktor4();
    void wypelnijtab0(double **tablice,int wier, int kol);
   lab4();

   void jakobian();

   void mnozenie_pc_det();

   void wektor();

   void macierz_wagi();

   void fullzad();

   ~lab4();
};

#endif // LAB4_H
