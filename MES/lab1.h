#ifndef LAB1_H
#define LAB1_H
#include "struktury.h"
using namespace std;



class lab1
{

public:
    Global_data gb;
    Grid g;
    string filename1="Test1_4_4.txt";
    string filename2="Test2_4_4_MixGrid.txt";
    string filename3="Test3_31_31_kwadrat.txt";
    string filename4="Test4_31_31_trapez.txt";



    void readfile();
    void write();

    ~lab1();
};

#endif // LAB1_H
