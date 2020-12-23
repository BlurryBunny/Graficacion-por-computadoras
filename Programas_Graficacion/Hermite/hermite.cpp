#include <iostream>
#include <armadillo>
#include <cmath>

int main(int argc, char *argv[]){

	//ecuacion parametrica de curvas de hermite
		// Q(t) = T * M * G .... 0 <= t <=1
	arma::fmat MH {{2,-2,1,1},
								{-3,3,-2,-1},
								{0,0,1,0},
								{1,0,0,0}};

	arma::frowvec P1 = {10,7,5};
	arma::frowvec P4 = {-1,-3,20};
	arma::frowvec R1 = {3,3,-2};
	arma::frowvec R4 = {-3,3,2};
	
	arma::fmat GH {4,3}
	GH.row(0)=P1;
	GH.row(1)=P4;
	GH.row(2)=R1;
	GH.row(3)=R4;
	
	float dt=0.1;
	
	for(float t=0.0; t<=1.0+dt; t+=dt){
		arma::frowvec T = {powf(t,3), powf(t,2), t, 1};
		arma::frowvec QT =T * MH *GH;
		cout<< " t= "<< t<< "... punto: "<<QT;
	}

}
