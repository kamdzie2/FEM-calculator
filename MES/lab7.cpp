#include "lab7.h"


lab7::lab7()
{

}

void lab7::konstruktor7()
{
    macierzag.wglobP=new double[gb.NodesNumber];
     macierzag.t0=new double[gb.NodesNumber];
     b=new double[gb.NodesNumber];

     macierzag.mglobH= new double*[gb.NodesNumber];
     for(int i=0;i<gb.NodesNumber;i++)
     {
         macierzag.mglobH[i]= new double [gb.NodesNumber];
     }

     macierzag.mglobC= new double*[gb.NodesNumber];
     for(int i=0;i<gb.NodesNumber;i++)
     {
         macierzag.mglobC[i]= new double [gb.NodesNumber];
     }

     A= new double*[gb.NodesNumber];
     for(int i=0;i<gb.NodesNumber;i++)
     {
        A[i]= new double [gb.NodesNumber];
     }

     for(int i=0;i<gb.NodesNumber;i++)
     {
         for(int j=0;j<gb.NodesNumber;j++)
         {
             macierzag.mglobH[i][j] =0;
             macierzag.mglobC[i][j] =0;

         }
         macierzag.wglobP[i]=0;
         macierzag.t0[i]=gb.InitialTemp;
    }
}

void lab7::sumowanieH()
{
    for(int t=0;t<gb.ElementsNumber;t++)
    {
        for(int i=0;i<4;i++)
        {
            for(int j=0;j<4;j++)
            {
                macierze[t].HFULL[i][j]=macierze[t].H[i][j]+macierze[t].Hbc[i][j];
            }
        }
    }
   // cout<<endl;
   // cout<<"macierz full"<<endl;

    for(int t=0;t<gb.ElementsNumber;t++)
    {
        //cout<<"ELEMENT "<<t<<endl;
        for(int i=0;i<4;i++)
        {
            for(int j=0;j<4;j++)
            {
              //  cout<<fixed<<setprecision(3)<< macierze[t].HFULL[i][j]<<' ';
            }
          //  cout<<endl;
        }
    }
}

void lab7::macierzC()
{
   // cout<<endl<<"El uni"<<endl;
    for(int i=0;i<w;i++)
    {
        for(int j=0;j<4;j++)
        {
            val.Nkszt[i][j]=N(values.wx[i%(Nw+1)],values.wx[(int)floor((double)i/(Nw+1))],j);
            //cout<<val.Nkszt[i][j]<<' ';
        }
       // cout<<endl;
    }

    for(int ttt=0;ttt<gb.ElementsNumber;ttt++)
    {
        for(int t=0;t<w;t++)
        {
             for(int i=0;i<4;i++)
             {
                 for(int j=0;j<4;j++)
                 {
                     macierze[ttt].C[i][j]+=gb.SpecificHeat*gb.Density*jakob_wyn[ttt][t]*val.Nkszt[t][i]*val.Nkszt[t][j]
                             *values.wa[t%(Nw+1)]*values.wa[(int)floor((double)t/(Nw+1))];
                 }
             }
         }
    }

   // cout<<endl<<"MACIERZ C"<<endl;

    for(int t=0;t<gb.ElementsNumber;t++)
    {
      //  cout<<"ELEMENT "<<t<<endl;
            for(int j=0;j<4;j++)
            {
               // cout<<'[';
                for(int jj=0;jj<4;jj++)
                {
               //      cout<<macierze[t].C[j][jj]<<' ';
                }
               // cout<<']'<<endl;
            }
    }
}

void lab7::agregacja()
{
    for(int t=0;t<gb.ElementsNumber;t++)
    {
        for(int i=0;i<4;i++)
        {
            for(int j=0;j<4;j++)
            {
               macierzag.mglobH[g.e[t].ID[i]-1][g.e[t].ID[j]-1]+=macierze[t].HFULL[i][j];
               macierzag.mglobC[g.e[t].ID[i]-1][g.e[t].ID[j]-1]+=macierze[t].C[i][j];
            }
            macierzag.wglobP[g.e[t].ID[i]-1]+=macierze[t].P[i];
        }
    }

    //cout<<"AGREGACJA H"<<endl;
    for(int i=0;i<gb.NodesNumber;i++)
    {
        for(int j=0;j<gb.NodesNumber;j++)
        {
            //cout<<fixed<<setprecision(3)<<macierzag.mglobH[i][j]<<' ';
        }
        //cout<<endl;
   }

   // cout<<"AGREGACJA C"<<endl;
    for(int i=0;i<gb.NodesNumber;i++)
    {
        for(int j=0;j<gb.NodesNumber;j++)
        {
           // cout<<fixed<<setprecision(3)<<macierzag.mglobC[i][j]<<' ';
        }
       // cout<<endl;
   }

    //cout<<"AGREGACJA P"<<endl;
    for(int i=0;i<gb.NodesNumber;i++)
    {
            //cout<<fixed<<setprecision(3)<<macierzag.wglobP[i]<<' ';
   }
}

void lab7::gaussianElimination(double** A, double* B, double* X, int N)
{
    for (int i = 0; i < N; ++i)
    {
        double maxEl = fabs(A[i][i]);
        int maxRow = i;
        for (int k = i + 1; k < N; ++k)
        {
            if (fabs(A[k][i]) > maxEl)
            {
                maxEl = fabs(A[k][i]);
                maxRow = k;
            }
        }
        for (int k = i; k < N; ++k)
        {
            swap(A[maxRow][k], A[i][k]);
        }
        swap(B[maxRow], B[i]);

        for (int k = i + 1; k < N; ++k)
        {
            double c = -A[k][i] / A[i][i];
            for (int j = i; j < N; ++j)
            {
                if (i == j)
                {
                    A[k][j] = 0;
                } else
                {
                    A[k][j] += c * A[i][j];
                }
            }
            B[k] += c * B[i];
        }
    }

    for (int i = N - 1; i >= 0; --i)
    {
        X[i] = B[i] / A[i][i];
        for (int k = i - 1; k >= 0; --k)
        {
            B[k] -= A[k][i] * X[i];
        }
    }
}

double* lab7::AiB(double t0[])
{
    int N=gb.NodesNumber;
    double** A = new double*[N];
        for (int i = 0; i < N; ++i)
        {
            A[i] = new double[N];
            for (int j = 0; j < N; ++j)
            {
                A[i][j] = macierzag.mglobH[i][j] + (macierzag.mglobC[i][j] / gb.SimulationStepTime);
            }
        }

        double* B = new double[N];
        for(int i=0;i<N;i++)
        {
             B[i] =macierzag.wglobP[i];
        }
        for (int i = 0; i < N; ++i) {
             for (int j = 0; j < N; ++j)
             {
                B[i] += (macierzag.mglobC[i][j]/gb.SimulationStepTime) * t0[j];

             }
        }

        for(int i=0;i<N;i++)
        {
             //cout<<macierzag.t0[i]<<' ';
        }
        double* t1 = new double[N];
        gaussianElimination(A, B, t1, N);
        tmin=t1[0];
        tmax=t1[0];
        for(int i=0;i<N;i++)
        {
            if(t1[i]<tmin)
            {
                tmin=t1[i];
            }
            else if(t1[i]>tmax)
            {
                tmax=t1[i];
            }
        }
       // cout << "Rozwiazanie t1:" << std::endl;
        for (int i = 0; i < N; ++i) {
           // cout << "t1[" << i << "] = " << t1[i] << std::endl;
        }
        return t1;
}


void lab7::uklad_rownan()
{
    int N=gb.NodesNumber;
    cout<<endl;
    for(int i=0;i<60;i++)
    {        
        string pom="foo";
        string rozszerzenie=".vtk";
        pom+=to_string(i);
        pom+=rozszerzenie;
        ofstream plik(pom);
        plik<<"# vtk DataFile Version 2.0"<<endl;
        plik<<"Unstructured Grid Example"<<endl;
        plik<<"ASCII"<<endl;
        plik<<"DATASET UNSTRUCTURED_GRID"<<endl<<endl;
        plik<<"POINTS "<<gb.NodesNumber<<" float"<<endl;
        for(int i=0;i<gb.NodesNumber;i++)
        {
            plik<<g.n[i].x<<' '<<g.n[i].y<<' '<<0<<endl;
        }
        plik<<endl;
        plik<<"CELLS "<<gb.ElementsNumber<<' '<<gb.ElementsNumber*5<<endl;
        for(int i=0; i<gb.ElementsNumber;i++)
        {
            plik<<4;
            for(int j=0; j<4;j++)
            {
                plik<<' '<<g.e[i].ID[j]-1;
            }
            plik<<endl;
        }
        plik<<endl;
        plik<<"CELL_TYPES "<<gb.ElementsNumber<<endl;
        for(int i=0;i<gb.ElementsNumber;i++)
        {
            plik<<9<<endl;
        }
        plik<<endl;
        plik<<"POINT_DATA "<<gb.NodesNumber<<endl;
        plik<<"SCALARS Temp float 1 "<<endl;
        plik<<"LOOKUP_TABLE default"<<endl;
        macierzag.t0=AiB(macierzag.t0);
        cout<<tmin<<' '<<tmax<<endl;
        for(int i=0;i<gb.NodesNumber;i++)
        {
            plik<<macierzag.t0[i]<<endl;
        }
        plik.close();
    }    
}


lab7::~lab7()
{
    delete[] macierzag.wglobP;
        macierzag.wglobP = nullptr;

        delete[] macierzag.t0;
        macierzag.t0 = nullptr;

        delete[] b;
        b = nullptr;

        for(int i=0; i<gb.NodesNumber; i++)
        {
            delete[] macierzag.mglobH[i];
        }
        delete[] macierzag.mglobH;
        macierzag.mglobH = nullptr;

        for(int i=0; i<gb.NodesNumber; i++)
        {
            delete[] macierzag.mglobC[i];
        }
        delete[] macierzag.mglobC;
        macierzag.mglobC = nullptr;

        for(int i=0; i<gb.NodesNumber; i++)
        {
            delete[] A[i];
        }
        delete[] A;
        A = nullptr;
}
