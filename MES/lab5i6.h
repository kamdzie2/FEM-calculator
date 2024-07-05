#ifndef LAB5I6_H
#define LAB5I6_H
#include "struktury.h"
#include "lab4.h"


class lab5i6: public lab4
{
    GaussIntegration pc;
    vector<GaussIntegration> SumaWektorowdoHBC;
    vector<GaussIntegration> WektoryP;
    double **detj;
    vector <double>wektordetj;

 public:

    lab5i6();


    double N(double ksi, double eta, int i);

    void konstruktor5();
    void wyzerujhbc();
    void wyzerujP();

    void tabelapc();

    void wektorP();

    void wypiszwektor();

    void licz_detj();

};

#endif // LAB5I6_H
