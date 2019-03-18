
#include<math.h>
#include <conio.h>
#include <iostream>
#include <fstream>
#include <cstdio>
#define eps 1e-5
#define n 5

double mass_x[n];
double mass_y[n];
double t[n];
double step = 0.1;

using namespace std;

double f(double x) 
{
  return 3*pow((x+1)/(x-1)*(x-1),(1./3));
}

double error(double x, double y) 
{
  return fabs(f(x)-y);
}

double d6(double x)
{
  return 640*(91*pow(x,6)+1638*pow(x,5)+4095*pow(x,4)+7020*pow(x,3)+5265*x*x+2430*x+405)/(243*pow(pow((x+1),17),(1./3))*pow(pow((x-1),20),(1./3)));
}

double Lagrange(double x)
{
  double res = 0;            
    for (int i = 0; i < n; i++)
    {
        double F = 1;
        for (int j = 0; j < n; j++)
        {
            if (i != j) F *= (x - mass_x[j]) / (mass_x[i] - mass_x[j]);
        }
        res += mass_y[i] * F;
    }
    return res;
}

double NewtonForward(double x)
{
  double res = mass_y[0], P, znam;            
    for (int i = 1; i < n; i++)
    {
      P = 0;
        for (int j = 0; j <= i; j++) 
        {
            znam = 1;                    
            for (int k = 0; k <= i; k++) 
                        if (k != j)
                            znam *= (mass_x[j] - mass_x[k]);                    
            P += mass_y[j] / znam; 
        }
        for (int k = 0; k < i; k++)
          P *= (x - mass_x[k]); 
        res += P;
    }
    return res;
}

double NewtonBackward(double x)
{
  double res = mass_y[n - 1], P, znam;
    for (int i = n-2; i >= 0; i--)
    {
      P = 0;
        for (int j = n-1; j >= i; j--)
        {
          znam = 1;
            for (int k = n-1; k >= i; k--)
                if (k != j)
                  znam *= (mass_x[j] - mass_x[k]);
            P += mass_y[j] / znam;
        }
        for (int k = n-1; k > i; k--)
          P *= (x - mass_x[k]);
        res += P;
    }
    return res;
} 
       
double Spline(double x, double mas_x[]) {
  struct spline_part
  {
    double a, b, c, d;
  };
  spline_part *splines;
  splines = new spline_part[n];
  for (int i = 0; i < n; ++i)
    splines[i].a = f(mas_x[i]);
  splines[0].c = 0;
  double *alpha = new double[n - 1];
  double *beta = new double[n - 1];
  double A, B, C, F, h_i, h_i1, z;
  alpha[0] = beta[0] = 0;
  for (int i = 1; i < n - 1; ++i)
  {
    h_i = mas_x[i] - mas_x[i - 1];
    h_i1 = mas_x[i + 1] - mas_x[i];
    A = h_i;
    C = 2 * (h_i + h_i1);
    B = h_i1;
    F = 6 * ((f(mas_x[i + 1]) - f(mas_x[i])) / h_i1 - (f(mas_x[i]) - f(mas_x[i - 1])) / h_i);
    z = (A * alpha[i - 1] + C);
    alpha[i] = -B / z;
    beta[i] = (F - A * beta[i - 1]) / z;
  }
  splines[n - 1].c = (F - A * beta[n - 2]) / (C + A * alpha[n - 2]);
  for (int i = n - 2; i > 0; --i)
    splines[i].c = alpha[i] * splines[i + 1].c + beta[i];

  for (int i = n - 1; i > 0; --i)
  {
    double h_i = mas_x[i] - mas_x[i - 1];
    splines[i].d = (splines[i].c - splines[i - 1].c) / h_i;
    splines[i].b = h_i * (2 * splines[i].c + splines[i - 1].c) / 6 + (f(mas_x[i]) - f(mas_x[i - 1])) / h_i;
  }
  spline_part *s;
  double r;
  if (x <= mas_x[0]) {
    s = splines + 1;
    r = mas_x[1];
  }
  else if (x >= mas_x[n - 1]) {
    s = splines + n - 1;
    r = mas_x[n - 1];
  }
  else
  {
    int i = 0, j = n - 1;
    while (i + 1 < j)
    {
      int k = i + (j - i) / 2;
      if (x <= mas_x[k])
        j = k;
      else
        i = k;
    }
    s = splines + j;
    r = mas_x[j];
  }
  double dx = (x - r);
  return s->a + s->b * dx + s->c * dx * dx / 2 + s->d * dx *dx * dx / 6;
}
double Magor(double x)
{
    double	w=1;
    for (int k=0; k<5; k++) w*=x-mass_x[k];																																	
	double magor =(fabs((d6(x)*w))/720);																																																																					//	magor  = 3.9608e-005;
	return magor;	
}

int main()
{
  double x=5.1, temp=5.05, z, w;
  double o=5.05;
  for (int i = 0; i < n; i++)
    {                
        mass_y[i] = f(x);
    	mass_x[i] = x;
        x += step;
    }
    x=5.1;
    for (int i = 0; i < n; i++)
    {
      cout<<"x="<<x<<"\tf(x)="<<f(mass_x[i])<<endl;
      cout<<"\tLG = "<<Lagrange(x)<<"\tError = "<<error(x,Lagrange(x))<<endl;
      cout<<"\tNF = "<<NewtonForward(x)<<"\tError = "<<error(x,NewtonForward(x))<<endl;
      cout<<"\tNB = "<<NewtonBackward(x)<<"\tError = "<<error(x,NewtonBackward(x))<<endl;
      cout<<"\tSP = "<<Spline(x, mass_x)<<"\tError = "<<error(x,Spline(x, mass_x))<<endl<<endl;
      x += step;
  }
  for(int i=0; i<n; i++)
  {
    t[i] = o;
    
    o = o+0.01;
  }
  o=0;
  for (z=5.05; z<=6.5; z+=0.01)
  {
     w=1;
    for (int k=0; k<5; k++) w*=temp-mass_x[k];
  }
  do
  {
    cout<<"Error LG ("<<temp<<") = "<<error(temp,Lagrange(temp))<<"\t";
    cout<<"Error NF ("<<temp<<") = "<<error(temp,NewtonForward(temp))<<"\t";
    cout<<"Error NB("<<temp<<") = "<<error(temp,NewtonBackward(temp))<<"\t";
    cout<<"Error SP ("<<temp<<") = "<<error(temp,Spline(temp, mass_x))<<endl;
    cout<<"Major = "<<Magor(temp)<<endl<<endl;
    temp+=0.01;
  }while(temp<=6.5);
  cout<<endl; 
  return 0;
}
