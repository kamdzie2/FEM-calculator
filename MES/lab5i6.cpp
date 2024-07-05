#include "lab5i6.h"


lab5i6::lab5i6()
{

}

void lab5i6::konstruktor5()
{
    detj= new double*[gb.ElementsNumber];
   for(int i=0;i<gb.ElementsNumber;i++)
   {
       detj[i]= new double [4];
   }

   for(int i=0;i<gb.ElementsNumber;i++)
   {
       for(int j=0;j<4;j++)
       {
          detj[i][j]= 0;
       }
  }

   for(int t=0;t<gb.ElementsNumber;t++)
   {

           for(int j=0;j<4;j++)
           {
               for(int jj=0;jj<4;jj++)
               {
                     macierze[t].Hbc[j][jj]=0;
                     macierze[t].C[j][jj]=0;
               }
           }
   }
}

double lab5i6::N(double ksi, double eta, int i)
{
    if (i==0)
    {
        return (1.0/4.0)*(1.0-ksi)*(1.0-eta);
    }
    else if (i==1)
    {
        return (1.0/4.0)*(1.0+ksi)*(1.0-eta);
    }
    else if (i==2)
    {
        return (1.0/4.0)*(1.0+ksi)*(1.0+eta);

    }else if (i==3)
    {
        return (1.0/4.0)*(1.0-ksi)*(1.0+eta);
    }
    else
    {
        cout<<"Blad"<<endl;
        return 0;
    }
}


void lab5i6::wyzerujhbc()
{
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
        {
            pc.pc[i][j]=0;
        }
    }
}
void lab5i6::wyzerujP()
{
    for(int i=0;i<4;i++)
    {
            pc.p[i]=0;
    }
}



void lab5i6::licz_detj()
{
    for(int t=0;t<gb.ElementsNumber;t++)
    {
        for(int i=0;i<4;i++)
        {

            int j=i;

            if(i==3)
            {
                j=-1;
            }
            if(g.n[g.e[t].ID[i]-1].BC!=0 && g.n[g.e[t].ID[j+1]-1].BC!=0 )
            {
                double pom=0;

                detj[t][i]=sqrt(pow((g.n[g.e[t].ID[i]-1].x-g.n[g.e[t].ID[j+1]-1].x),2)+pow((g.n[g.e[t].ID[i]-1].y-g.n[g.e[t].ID[j+1]-1].y),2));

                pom=detj[t][i]/2;
                detj[t][i]=pom;
            }
            else
            {
                detj[t][i]=0;
            }

           // cout<< detj[t][i]<<' ';

        }
        //cout<<endl;
    }
}



void lab5i6::tabelapc()
{
  // cout<<"MACIERZE HBC"<<endl;
   wyzerujhbc();
   for(int tt=0; tt<gb.ElementsNumber;tt++)
   for(int pw=0;pw<4;pw++)
   {
        for(int t=0;t<Nw+1;t++)
        {
            for(int i=0;i<4;i++)
            {
                for(int j=0;j<4;j++)
                {

                    if(pw==0)
                    {
                      pc.pc[i][j]+=(values.wa[t]*gb.Alfa*(N(values.wx[t],-1,i)*N(values.wx[t],-1,j)));
                      if(t==Nw)
                      {
                          pc.pc[i][j]*=detj[tt][pw];
                       macierze[tt].Hbc[i][j]+=pc.pc[i][j];
                      }
                    }
                    else if(pw==1)
                    {
                         pc.pc[i][j]+=(values.wa[t]*gb.Alfa*(N(1,values.wx[t],i)*N(1,values.wx[t],j)));
                         if(t==Nw)
                         {
                             pc.pc[i][j]*=detj[tt][pw];
                               macierze[tt].Hbc[i][j]+=pc.pc[i][j];
                         }
                    }
                    else if(pw==2)
                    {
                        pc.pc[i][j]+=(values.wa[t]*gb.Alfa*(N(values.wx[t],1,i)*N(values.wx[t],1,j)));
                            if(t==Nw)
                            {
                                pc.pc[i][j]*=detj[tt][pw];
                                  macierze[tt].Hbc[i][j]+=pc.pc[i][j];
                            }
                    }
                    else if(pw==3)
                    {
                        pc.pc[i][j]+=(values.wa[t]*gb.Alfa*(N(-1,values.wx[t],i)*N(-1,values.wx[t],j)));
                            if(t==Nw)
                            {
                                pc.pc[i][j]*=detj[tt][pw];
                                  macierze[tt].Hbc[i][j]+=pc.pc[i][j];
                            }
                    }
                }
            }
        }
        wyzerujhbc();
    }

   for(int t=0;t<gb.ElementsNumber;t++)
   {
       //cout<<"ELEMENT "<<t<<endl;
           for(int j=0;j<4;j++)
           {
               //cout<<'[';
               for(int jj=0;jj<4;jj++)
               {
                    //cout<<macierze[t].Hbc[j][jj]<<' ';
               }
              // cout<<']'<<endl;
           }
   }
}

void lab5i6::wektorP()
{
    //cout<<"Wektory P"<<endl;
    double temp=gb.Tot;
    wyzerujP();
    for(int tt=0;tt<gb.ElementsNumber;tt++)
    {
        for(int pw=0;pw<4;pw++)
        {
             for(int t=0;t<Nw+1;t++)
             {
                 for(int i=0;i<4;i++)
                 {
                     if(pw==0)
                     {
                       pc.p[i]+=(values.wa[t]*gb.Alfa*temp*(N(values.wx[t],-1,i)));
                       if(t==Nw)
                       {
                           pc.p[i]*=detj[tt][pw];
                        macierze[tt].P[i]+=pc.p[i];
                       }
                     }
                     else if(pw==1)
                     {
                          pc.p[i]+=(values.wa[t]*gb.Alfa*temp*(N(1,values.wx[t],i)));
                          if(t==Nw)
                          {
                              pc.p[i]*=detj[tt][pw];
                           macierze[tt].P[i]+=pc.p[i];
                          }
                     }
                     else if(pw==2)
                     {
                         pc.p[i]+=(values.wa[t]*gb.Alfa*temp*(N(values.wx[t],1,i)));
                         if(t==Nw)
                         {
                             pc.p[i]*=detj[tt][pw];
                          macierze[tt].P[i]+=pc.p[i];
                         }
                     }
                     else if(pw==3)
                     {
                         pc.p[i]+=(values.wa[t]*gb.Alfa*temp*(N(-1,values.wx[t],i)));
                         if(t==Nw)
                         {
                             pc.p[i]*=detj[tt][pw];
                          macierze[tt].P[i]+=pc.p[i];
                         }
                     }
                 }
             }
             wyzerujP();
         }
    }

    for(int t=0;t<gb.ElementsNumber;t++)
    {
        //cout<<"ELEMENT "<<t<<endl;
            for(int j=0;j<4;j++)
            {
               // cout<<'[';
                //cout<<macierze[t].P[j]<<' ';
                //cout<<']'<<endl;
            }
    }
}

void lab5i6::wypiszwektor()
{

    for(const auto& tablica:SumaWektorowdoHBC)
    {
        for(int j=0;j<GaussIntegration::Nn;j++)
        {
            for(int jj=0; jj<GaussIntegration::Nn; jj++)
            {
               // cout<<tablica.pc[j][jj]<<' ';
            }
           // cout<<endl;
        }
      //  cout<<endl;

    }
}
