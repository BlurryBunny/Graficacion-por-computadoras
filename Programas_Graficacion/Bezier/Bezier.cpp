#include <iostream>
#include <armadillo>
#include <cmath>

using namespace std;
int main(int argc, char *argv[]){

	//ecuacion parametrica de curvas de hermite
		// Q(t) = T * M * G .... 0 <= t <=1
	arma::fmat MB = {{-1,3,-3,1},
								{3,-6,3,0},
								{-3,3,0,0},
								{1,0,0,0}};

	arma::frowvec P1 = {0,0,0};
	arma::frowvec P2 = {3,3,0};
	arma::frowvec P3 = {6,3,0};
	arma::frowvec P4 = {10,0,0};
	
	arma::fmat GB (4,3);
	GB.row(0)=P1;
	GB.row(1)=P2;
	GB.row(2)=P3;
	GB.row(3)=P4;
	
	float dt=0.1;
	
	for(float t=0.0; t<=1.0+dt; t+=dt){
		arma::frowvec T = {powf(t,3), powf(t,2), t, 1};
		arma::frowvec QT =T * MB *GB;
		cout<< " t= "<< t<< "... punto: "<<QT;
	}

}
