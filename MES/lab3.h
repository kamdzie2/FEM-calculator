#ifndef LAB3_H
#define LAB3_H
#include "struktury.h"
#include "lab2.h"


class lab3: public lab2
{

public:
    int w=(Nw+1)*(Nw+1);
    int k=4;

    ElementUniwersalny val;

    double fksi1( double eta);

    double fksi2(double eta);

    double fksi3(double eta);

    double fksi4(double eta);

    double fketa1( double ksi);

    double fketa2(double ksi);

    double fketa3(double ksi);

    double fketa4(double ksi);

    void matrix();

    void wylicz();


    void wypisz();

    ~lab3();

};

#endif // LAB3_H
