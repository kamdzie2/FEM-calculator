#include "lab4.h"


lab4::lab4()
{

}

void lab4::konstruktor4()
{
    macierze=new Element[gb.ElementsNumber];

    //tablica wartości n/dx i n/dy
   valXY.dx= new double*[w];
   valXY.dy= new double*[w];
   for(int i=0;i<w;i++)
   {
       valXY.dx[i]= new double [k];
       valXY.dy[i]= new double [k];
   }

   for(int i=0;i<w;i++)
   {
       for(int j=0;j<k;j++)
       {
           valXY.dx[i][j]= 0;
           valXY.dy[i][j] =0;
       }
  }

   //jakobian końcowy
   jakob_wyn= new double*[gb.ElementsNumber];
   for(int i=0;i<gb.ElementsNumber;i++)
   {
       jakob_wyn[i]= new double [w];
   }
   wypelnijtab0(jakob_wyn,gb.ElementsNumber,w);

   //jakobian
   jk= new Jakobian*[gb.ElementsNumber];
   for(int i=0;i<gb.ElementsNumber;i++)
   {
       jk[i]= new Jakobian [w];
   }

   for(int i=0;i<gb.ElementsNumber;i++)
   {
       for(int j=0;j<w;j++)
       {
           jk[i][j].dl= 0;
           jk[i][j].gl= 0;
           jk[i][j].dp= 0;
           jk[i][j].gp= 0;
       }
  }

   //jakobian odwrócony
   jakob_odw= new double*[gb.ElementsNumber];
   for(int i=0;i<gb.ElementsNumber;i++)
   {
       jakob_odw[i]= new double [w];
   }
   wypelnijtab0(jakob_odw,gb.ElementsNumber,w);

}

void lab4::wypelnijtab0(double **tablice,int wier, int kol)
{
    for(int i=0;i<wier;i++)
    {
        for(int j=0;j<kol;j++)
        {
            tablice[i][j]= 0;
        }
   }
}

void lab4::jakobian()
{
    for(int i=0;i<gb.ElementsNumber;i++)
    {
       // cout<<"Element "<<i<<endl;
        for(int j=0;j<w;j++)
        {
            for(int jj=0; jj<4; jj++)
            {
                jk[i][j].gl+=val.dKsi[j][jj]*g.n[g.e[i].ID[jj]-1].x;
                jk[i][j].gp+=val.dKsi[j][jj]*g.n[g.e[i].ID[jj]-1].y;
                jk[i][j].dl+=val.dEta[j][jj]*g.n[g.e[i].ID[jj]-1].x;
                jk[i][j].dp+=val.dEta[j][jj]*g.n[g.e[i].ID[jj]-1].y;
            }
        }
    }

    for(int i=0;i<gb.ElementsNumber;i++)
    {
        //cout<<"Element "<<i<<endl;
        for(int j=0;j<w;j++)
        {
          //  cout<<"pc"<<j+1<<endl;
          // cout<<"["<<jk[i][j].gl<<"]"<<"["<<jk[i][j].gp<<"]"<<endl;
          // cout<<"["<<jk[i][j].dl<<"]"<<"["<<jk[i][j].dp<<"]"<<endl;
           jakob_wyn[i][j]=jk[i][j].gl*jk[i][j].dp- jk[i][j].dl*jk[i][j].gp;
          // cout<<"Jakobian= "<<jakob_wyn[i][j]<<endl;
        }
    }
}

void lab4::mnozenie_pc_det()
{
    for(int i=0;i<gb.ElementsNumber;i++)
    {
        for(int j=0;j<w;j++)
            jakob_odw[i][j]=1/jakob_wyn[i][j];
    }

   for(int i=0;i<gb.ElementsNumber;i++)
   {
       for(int j=0;j<w;j++)
       {
           double jkgl=jk[i][j].gl;
           double jkgp=jk[i][j].gp;
           double jkdl=jk[i][j].dl;
           double jkdp=jk[i][j].dp;

           jk[i][j].gl=jkdp*jakob_odw[i][j];
           jk[i][j].gp=-(jkdl*jakob_odw[i][j]);
           jk[i][j].dl=-(jkgp*jakob_odw[i][j]);
           jk[i][j].dp=jkgl*jakob_odw[i][j];
       }
   }
   for(int i=0;i<gb.ElementsNumber;i++)
   {
      // cout<<"Element "<<i<<endl;
       for(int j=0;j<w;j++)
       {
          // cout<<"pc"<<i+1<<endl;
         // cout<<"["<<jk[i][j].gl<<"]"<<"["<<jk[i][j].gp<<"]"<<endl;
         // cout<<"["<<jk[i][j].dl<<"]"<<"["<<jk[i][j].dp<<"]"<<endl;
       }
   }

   for(int t=0;t<gb.ElementsNumber;t++)
   {
           for(int j=0;j<N;j++)
           {
               for(int jj=0;jj<N;jj++)
               {
                     macierze[t].H[j][jj]=0;
               }
           }
   }

   //cout<<"TABLICA XY"<<endl;
   for(int t=0;t<gb.ElementsNumber;t++)
   {
       for(int i=0;i<w;i++)
       {
           for(int j=0;j<k;j++)
           {
               valXY.dx[i][j]=jk[t][i].gl*val.dKsi[i][j]+jk[t][i].dl*val.dEta[i][j];
               valXY.dy[i][j]=jk[t][i].gp*val.dKsi[i][j]+jk[t][i].dp*val.dEta[i][j];
              // cout<<valXY.dx[i][j]<<' ';
           }
          // cout<<endl;
       }

       for(int i=0;i<w;i++)
       {
          // cout<<"wiersz "<<i<<endl;
           for(int j=0;j<4;j++)
           {
               for(int jj=0;jj<4;jj++)
               {
                   macierze[t].H[j][jj]+=(((valXY.dx[i][j]*valXY.dx[i][jj])+(valXY.dy[i][j]*valXY.dy[i][jj]))*gb.Conductivity*jakob_wyn[t][i]*values.wa[i%(Nw+1)]*values.wa[(int)floor((double)i/(Nw+1))]);
               }
               //cout<<"x= "<<valXY.dx[i][j]<<' ';
               //cout<<endl;
           }

       }
   }


   for(int t=0;t<gb.ElementsNumber;t++)
   {
      // cout<<"ELEMENT "<<t<<endl;
           for(int j=0;j<4;j++)
           {
               //cout<<'[';
               for(int jj=0;jj<4;jj++)
               {
                   // cout<<macierze[t].H[j][jj]<<' ';
               }
              // cout<<']'<<endl;
           }
   }   
}

void lab4::fullzad()
{   
    readfile();
    write();
    read_N_number();
    matrix();
    wylicz();
    konstruktor4();
    jakobian();
    mnozenie_pc_det();
}

lab4::~lab4()
{
    for (int i = 0; i < gb.ElementsNumber; i++)
    {
        delete[] jakob_odw[i];
    }
    delete[] jakob_odw;

    for (int i = 0; i < gb.ElementsNumber; i++)
    {
        delete[] jk[i];
    }
    delete[] jk;

    for (int i = 0; i < gb.ElementsNumber; i++)
    {
        delete[] jakob_wyn[i];
    }
    delete[] jakob_wyn;

    for (int i = 0; i < w; i++)
    {
        delete[] valXY.dy[i];
    }
    delete[] valXY.dy;

    for (int i = 0; i < w; i++)
    {
        delete[] valXY.dx[i];
    }
    delete[] valXY.dx;

    delete[] macierze;
}
