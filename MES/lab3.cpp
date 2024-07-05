#include "lab3.h"

double lab3::fksi1( double eta)
{
    return ((-1.0/4.0)*(1.0-eta));
}
double  lab3::fksi2(double eta)
{
    return ((1.0/4.0)*(1.0-eta));
}
double  lab3::fksi3(double eta)
{
    return ((1.0/4.0)*(1.0+eta));
}
double  lab3::fksi4(double eta)
{
    return ((-1.0/4.0)*(1.0+eta));
}
double  lab3::fketa1( double ksi)
{
    return -(1.0/4.0)*(1.0-ksi);
}
double  lab3::fketa2(double ksi)
{
    return -(1.0/4.0)*(1.0+ksi);
}
double  lab3::fketa3(double ksi)
{
    return (1.0/4.0)*(1.0+ksi);
}
double  lab3::fketa4(double ksi)
{
    return (1.0/4)*(1.0-ksi);
}


void  lab3::matrix()
{
        val.dKsi= new double*[w];
        val.dEta= new double*[w];
        val.Nkszt= new double *[w];
        for(int i=0;i<w;i++)
        {
            val.dKsi[i]= new double [k];
            val.dEta[i]= new double [k];
            val.Nkszt[i]= new double [k];
        }
}

void  lab3::wylicz()
{
    int pom=0;
    int licznik=0;

    for(int i=0; i<k;i++)
    {
        for(int j=0; j<w;j++)
        {
            if(i==0)
            {
                val.dKsi[j][i]=fksi1(values.wx[licznik]);
                pom++;
                if(pom==(Nw+1))
                {
                    licznik++;
                    pom=0;
                }
            }
            else if(i==1)
            {
                val.dKsi[j][i]=fksi2(values.wx[licznik]);
                pom++;
                if(pom==(Nw+1))
                {
                    licznik++;
                    pom=0;
                }
            }
            else if(i==2)
            {
                val.dKsi[j][i]=fksi3(values.wx[licznik]);
                pom++;
                if(pom==(Nw+1))
                {
                    licznik++;
                    pom=0;
                }
            }
            else if(i==3)
            {
                val.dKsi[j][i]=fksi4(values.wx[licznik]);
                pom++;
                if(pom==(Nw+1))
                {
                    licznik++;
                    pom=0;
                }
            }
        }
        licznik=0;
    }

    for(int i=0; i<k;i++)
    {
        for(int j=0; j<w;j++)
        {
            if(i==0)
            {
                val.dEta[j][i]=fketa1(values.wx[licznik]);
                licznik++;
                if(licznik==(Nw+1))
                {
                    licznik=0;
                }
            }
            else if(i==1)
            {
                val.dEta[j][i]=fketa2(values.wx[licznik]);
                licznik++;
                if(licznik==(Nw+1))
                {
                    licznik=0;
                }
            }
            else if(i==2)
            {
                val.dEta[j][i]=fketa3(values.wx[licznik]);
                licznik++;
                if(licznik==(Nw+1))
                {
                    licznik=0;
                }
            }
            else if(i==3)
            {
                val.dEta[j][i]=fketa4(values.wx[licznik]);
                licznik++;
                if(licznik==(Nw+1))
                {
                    licznik=0;
                }
            }
        }
        licznik=0;
    }
}

void  lab3::wypisz()
{
     cout<<"wyniki dla ksi:"<<endl;
    for(int i=0;i<w;i++)
    {
        for(int j=0;j<k;j++)
        {
            cout<<val.dKsi[i][j]<<' ';
        }
        cout<<endl;

    }

    cout<<"\n wyniki dla eta:"<<endl;
   for(int i=0;i<w;i++)
   {
       for(int j=0;j<k;j++)
       {
           cout<<val.dEta[i][j]<<' ';
       }
       cout<<endl;

   }
}
 lab3::~lab3()
{
     for (int i = 0; i < w; i++)
         {
             delete[] val.Nkszt[i];
         }
         delete[] val.Nkszt;


         for (int i = 0; i < w; i++)
         {
             delete[] val.dEta[i];
         }
         delete[] val.dEta;


         for (int i = 0; i < w; i++)
         {
             delete[] val.dKsi[i];
         }
         delete[] val.dKsi;
}
