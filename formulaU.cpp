//Equation formula for the analytical model sPASS
//g++ -std=c++1y formulaU.cpp -o formulaU 

#include <iostream>
#include <math.h>
#include <algorithm>
int main()
{
	//n = total no of partitionsa
	//u = number of super pertitoins while v = number of sub partitions
	
	int u, v, n;
	double uStar, alpha, M, valComp1, valComp2, Cb, tORAM, tHE , b, Cs, eMinus;
	std::cout<<"Input Progression: u v b Cb Cs M tORAM tHE"<<std::endl;
	std::cin>>u>>v>>b>>Cb>>Cs>>M>>tORAM>>tHE;

	 n = u*v;
	alpha = ((Cb * tORAM) / tHE);
	
	eMinus = 1/ (exp(1+ sqrt(1-(1/alpha)))) ;
	
	valComp1 = n * eMinus; 
	
	valComp2 = M / (b * Cs * log(v));

	std::cout<<valComp1<<" "<<valComp2<<std::endl;
uStar = std::min(valComp1, valComp2);
std::cout<<uStar<<std::endl;
}
