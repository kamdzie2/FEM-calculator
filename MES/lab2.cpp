#include "lab2.h"

double lab2::f1(double x)
{
    return 5*pow(x,2)+3*x+6;
}

double lab2::f2(double x, double y)
{
    return 5*pow(x,2)*pow(y,2)+3*x*y+6;
}

void lab2:: read_N_number()
{
    if(Nw>0 && Nw<4)
    {
        if(Nw==1)
        {
            int lp=2;
            values.k=new int[lp];
            values.k[0]=0;
            values.k[1]=1;

            values.wx=new double[lp];
            values.wx[0]=-(1.0/sqrt(3.0));
            values.wx[1]=(1.0/sqrt(3.0));

            values.wa=new double[lp];
            values.wa[0]=1.0;
            values.wa[1]=1.0;
        }
        else if(Nw==2)
        {
            int lp=3;
            values.k=new int[lp];
            values.k[0]=0.0;
            values.k[1]=1;
            values.k[2]=2;

            values.wx=new double[lp];
            values.wx[0]=-(sqrt(3.0/5.0));
            values.wx[1]=0.0;
            values.wx[2]=sqrt(3.0/5.0);

            values.wa=new double[lp];
            values.wa[0]=5.0/9.0;
            values.wa[1]=8.0/9.0;
            values.wa[2]=5.0/9.0;
        }
        else if(Nw==3)
        {
            int lp=4;
            values.k=new int[lp];
            values.k[0]=0;
            values.k[1]=1;
            values.k[2]=2;
            values.k[3]=3;

            values.wx=new double[lp];
            values.wx[0]=sqrt((3.0 / 7.0) + (2.0 / 7.0) * sqrt(6.0 / 5.0)) * -1.0;
            values.wx[1]=sqrt((3.0 / 7.0) - (2.0 / 7.0) *sqrt(6.0 / 5.0)) * -1.0;
            values.wx[2]=sqrt((3.0 / 7.0) - (2.0 / 7.0) * sqrt(6.0 / 5.0));
            values.wx[3]=sqrt((3.0 / 7.0) + (2.0 / 7.0) * sqrt(6.0 / 5.0));

            values.wa=new double[lp];
            values.wa[0]=(18.0 - sqrt(30.0)) / 36;
            values.wa[1]=(18.0 + sqrt(30.0)) / 36;
            values.wa[2]=(18.0 + sqrt(30.0)) / 36;
            values.wa[3]=(18.0 - sqrt(30.0)) / 36;
        }
    }

}
void lab2:: case_1d(int Nwez)
{
    double result_f1=0;
        for(int i=0; i<Nwez; i++)
        {
            result_f1+=values.wa[i]*f1(values.wx[i]);
        }
   cout<<"Wynik f1 dla "<<Nwez<<" wezlow= "<<result_f1<<endl;
}

void lab2:: case_2d(int Nwez)
{
    double result_f2=0;
    for (int i = 0; i < Nwez; i++)
    {
        for (int j = 0; j < Nwez; j++)
        {
            result_f2 += values.wa[i] * values.wa[j] * f2(values.wx[i], values.wx[j]);
        }
    }
    cout<<"Wynik f2 dla "<<Nwez<<" wezlow= "<<result_f2<<endl;
}

lab2:: ~lab2()
{
    delete[] values.k;
    delete[] values.wa;
    delete[] values.wx;
}
