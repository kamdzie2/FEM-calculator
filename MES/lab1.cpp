#include "lab1.h"


void lab1::readfile()
{
    ifstream file(filename2);
    string line;
    char smieciCHAR;
    int smieciINT;

    if (!file)
    {
    cout<< "Nie mozna otworzyc pliku" << endl;
    }
    else
    {
        while(file >> line)
        {
            if (line=="SimulationTime" )
            {
                file >> gb.SimulationTime;
            }
            else if (line=="SimulationStepTime")
            {
                file >> gb.SimulationStepTime;
            }
            else if (line=="Conductivity")
            {
                file >> gb.Conductivity;
            }
            else if (line=="Alfa")
            {
                file >> gb.Alfa;
            }
            else if (line=="Tot")
            {
                file >> gb.Tot;
            }
            else if (line=="InitialTemp")
            {
                file >> gb.InitialTemp;
            }
            else if (line=="Density")
            {
                file >> gb.Density;
            }
            else if (line=="SpecificHeat")
            {
                file >> gb.SpecificHeat;
            }
            else if (line=="Nodes")
            {
                file>>line;
                if(line=="number")
                file >> gb.NodesNumber;

            }
            else if (line=="Elements")
            {
                file>>line;
                if(line=="number")
                file >> gb.ElementsNumber;
            }
            else if(line=="*Node")
            {
               g.n = new Node[gb.NodesNumber];
               for(int i=0;i<gb.NodesNumber;i++)
               {
                   file>>smieciINT>>smieciCHAR>>g.n[i].x>>smieciCHAR>>g.n[i].y;
               }

            }
            else if(line=="*Element")
            {
            }
            else if(line=="type=DC2D4")
            {
                g.e = new Element[gb.ElementsNumber];
                for(int i=0;i<gb.ElementsNumber;i++)
                {
                    file>>smieciINT;
                    for(int j=0; j<N;j++)
                    {
                        file>>smieciCHAR>>g.e[i].ID[j];
                    }
                }
             }else if(line=="*BC")
            {
                int pom;
                file>>pom;
                if(pom==1)
                {
                    g.n[0].BC=pom;
                }
                for(int i=1;i<gb.NodesNumber;i++)
                {
                    int pom;
                    char smieciCHAR;
                    file>>smieciCHAR;
                    file>>pom;
                        g.n[pom-1].BC=pom;
                }

            }
        }
    }
    file.close();
}
void lab1:: write()
{
    cout << "SimulationTime: " << gb.SimulationTime << endl;
    cout << "SimulationStepTime: " << gb.SimulationStepTime << endl;
    cout << "Conductivity: " << gb.Conductivity << endl;
    cout << "Alfa: " << gb.Alfa << endl;
    cout << "Tot: " << gb.Tot << endl;
    cout << "InitialTemp: " << gb.InitialTemp << endl;
    cout << "Density: " << gb.Density << endl;
    cout << "SpecificHeat: " << gb.SpecificHeat << endl;
    cout << "Nodes number: " << gb.NodesNumber << endl;
    cout << "Elements number: " << gb.ElementsNumber << endl;

    cout << fixed << setprecision(11);
    cout<<endl;
    cout<<"struktura Grid:"<<endl;
    for(int i=0; i<gb.NodesNumber;i++)
    {
        cout<<"Node "<<i<<": ";
        cout<<"x"<<i<<"="<<g.n[i].x<<" y"<<i<<"="<<g.n[i].y;
        cout<<endl;
    }

    cout<<endl;
    for(int i=0; i<gb.ElementsNumber;i++)
    {
        cout<<"ID"<<i<<": {";
        for(int j=0; j<N;j++)
        {
            cout<<g.e[i].ID[j]<<' ';
        }
        cout<<"}";
        cout<<endl;
    }

    cout<<endl;
    cout<<"BC:"<<endl;
    for(int i=0; i<gb.NodesNumber;i++)
    {
        cout<<g.n[i].BC<<", ";
    }


}

lab1::~lab1()
{
    delete[] g.n;
    delete[] g.e;

}
