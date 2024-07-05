#ifndef LAB7_H
#define LAB7_H
#include "lab5i6.h"
#include "struktury.h"

class lab7:public lab5i6
{
public:
    double tmin, tmax;
    double **A;
    double *b;
    lab7();
    ~lab7();
    SOE macierzag;
    void konstruktor7();
    void macierzC();
    void gaussianElimination(double** A, double* B, double* X, int N);
    void uklad_rownan();
    void agregacja();
    void sumowanieH();
    double* AiB(double t0[]);
};

#endif // LAB7_H
