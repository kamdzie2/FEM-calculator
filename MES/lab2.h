#ifndef LAB2_H
#define LAB2_H
#include "struktury.h"
#include "lab1.h"
class lab2 : public lab1
{

public:
     Nwezlow values;
    int Nw=1;
    double f1(double x);

    double f2(double x, double y);

    void read_N_number();

    void case_1d(int Nwez);

    void case_2d(int Nwez);

    ~lab2();
};

#endif // LAB2_H
