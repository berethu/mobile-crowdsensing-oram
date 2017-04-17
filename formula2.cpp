//g++ -std=c++1y formula2.cpp -o formula2

//This is newton raphson's methof for solving numerical equation
#include <iostream>
#include <math.h>
#include <algorithm>
#include <boost/math/constants/constants.hpp>
using namespace std;

//what preduces nan is usually when the number of initial guess w0 is smaller or larger and therefore cannot be calculated
double pi = boost::math::double_constants::pi;
	// double s = 5; //standard deviation
	// int v = 100;

double f(double w, double s, int v)
{
	double a = (pow((pow((1-1/v),(v/(s*sqrt(2*pi))))),w)) ;
		double b = ((w* - log(1-1/v*(v/(s*sqrt(2*pi))) )) +1);
	return (a*b) -1; //log(w) - 3;   //input fx here
}

double df(double w, double s, int v)
{
	double a = pow((   pow (((v-1)/v),(v/(sqrt(2 * pi)*s)))  ), w) ;
	double b = ( log( pow(((v-1)/v), (v/(sqrt(2 * pi) *s))) )- log(1-(1/(sqrt(2 * pi) *s))) * (w* log(pow(((v-1)/v),(v/(sqrt(2 * pi)*s))) )+1));
	return a *b ; //1/w;   //input derivative here
}

int main()
{
	//Sample test variable is 2 0.0001 100

int maxIteration = 1000;
double fwDivFprimeW, w0, w1, tolerance = 0.0001, s=31.8332;
int v = 2;
cout<<"Input w0, allowed error, maximum iterations, s, v"<<endl;
cin>>w0;//>>tolerance>>maxIteration>>s>>v;

for(int itr = 1; itr<=maxIteration; itr++)
{
	fwDivFprimeW= f(w0, s, v)/df(w0, s, v);
	w1 = w0-fwDivFprimeW;
	cout<<"At iteration "<<itr<<" value of w = "<<w1 <<endl;;
	if (fabs(fwDivFprimeW) < tolerance)
	{
		cout<<"\nAfter "<<itr<<" iterations, w = "<<w1<<"\n"<<endl;
		// cout<<s<<endl;
		return 0; //success
	}
w0 = w1;	
}
cout<<"No convergence or insufficient number of iterations"<<endl;
return 1; //failure


} 