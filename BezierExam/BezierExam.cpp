#include <iostream>
#include <armadillo>
#include <cmath>
#define PI 3.1416
using namespace std;
float calculate_distance(arma::fcolvec v1, arma::fcolvec v2);

int main(int argc, char *argv[]){

	//suponiendo que el esquema tiene un angulo de 45 grados
	arma::fcolvec v1 ={3,4,1,1};
	arma::fcolvec v2 ={7,10,1,1};

	
	float distance=calculate_distance(v1,v2);
	
	//dependiendo del angulo cambia
	float altura=cos(45 * PI /180)*distance;
	
	arma::frowvec pc1 = {v1[0]-distance,v1[1]-altura,v1[2]};
	arma::frowvec pc2 ={v1[0],v1[1]+(distance/4),v1[2]};
	arma::frowvec pc3 ={v2[0]+(distance/4),v2[1],v2[2]};
	arma::frowvec pc4 ={v1[0]+distance,v1[1]+altura,v1[2]};
	
	
		//ecuacion parametrica de curvas de hermite
		// Q(t) = T * M * G .... 0 <= t <=1
	arma::fmat MB = {{-1,3,-3,1},
								{3,-6,3,0},
								{-3,3,0,0},
								{1,0,0,0}};
								
	arma::fmat GB (4,3);
	GB.row(0)=pc1;
	GB.row(1)=pc2;
	GB.row(2)=pc3;
	GB.row(3)=pc4;
	
	float dt=0.1;
	
	for(float t=0.0; t<=1.0+dt; t+=dt){
		arma::frowvec T = {powf(t,3), powf(t,2), t, 1};
		arma::frowvec QT =T * MB *GB;
		cout<< " t= "<< t<< "... punto: "<<QT;
	}

}

float calculate_distance(arma::fcolvec v1, arma::fcolvec v2){
	return( abs(sqrt(powf(v2[0]-v1[0],2) + powf(v2[2]-v1[2],2) + powf(v2[2]-v1[2],2))) );
}

