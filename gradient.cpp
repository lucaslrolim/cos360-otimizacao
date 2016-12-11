#include <iostream>
#include <math.h>

using namespace std;

double gradiente_f1(double x, double y,int variavel){
    double gradiente[2];
    gradiente[0] = -1*(x*(1-y)-2*x*y*(1-y))/(x*(1-x)*y*(1-y));
    gradiente[1] = -1*(y*(1-x)-2*y*x*(1-x))/(x*(1-x)*y*(1-y));
    return gradiente[variavel];
}

double gradiente_f2(double x, double y,int variavel){
    double gradiente[2];
    gradiente[0] = 1/2*(gradiente_f1(x, y, 0))/2*sqrt(log(x*(1-x)*y*(1-y)));
    gradiente[1] = 1/2*(gradiente_f1(x, y, 1))/2*sqrt(log(x*(1-x)*y*(1-y)));
    return gradiente[variavel];
}

int main()
{
 cout << gradiente_f1(13,20,1) << endl;

 //  Aplicação do método do gradiente para F1
 0;
}