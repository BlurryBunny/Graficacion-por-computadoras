#include <iostream>
#include <armadillo>
#include <cmath>

#define PI 3.1416

using namespace std;
arma::dmat rot_v1v2(arma::fcolvec v1,arma::fcolvec v2, double angle){
	
	//Paso 1
	//Comenzamos con la rotacion del objeto.
	//Rotacion en X
	arma::dmat T = {
									{1,0,0,-v1[0]},
									{0,1,0,-v1[1]},
									{0,0,1,-v1[2]},
									{0,0,0,1}
									};

	//Paso 2
	double D1= sqrt(powf(v2[0]-v1[0],2) + powf(v2[2]-v1[2],2));
	
	arma::dmat Ry = {
									{(v2[0]-v1[0])/D1,0,(v2[2]-v1[2])/D1,0},
									{0,1,0,0},
									{-((v2[2]-v1[2])/D1),0,(v2[0]-v1[0])/D1,0},
									{0,0,0,1}
									};
									
	//Paso 3	
	//Rz(-alpha)
	//D2 = sqtr(x2p^2)
	double D2 = sqrt(powf(v2[0]-v1[0],2) + powf(v2[1]-v1[1],2) + powf(v2[2]-v1[2],2));
	
	arma::dmat Rz = {
									{D1/D2,0,(v2[1]-v1[1])/D2,0},
									{-((v2[1]-v1[1])/D2),0,D1/D2,0},
									{0,0,1,0},
									{0,0,0,1}
									};

	//Paso 4
	//Rx(angle);
	
	arma::dmat Rx = {
									{1,0,0,0},
									{0,cos(angle*PI/ 180.0),-sin(angle*PI/ 180.0),0},
									{0,sin(angle*PI/ 180.0),cos(angle*PI/ 180.0),0},
									{0,0,0,1}
									};
	
	//Ahora tenemos que regresar todo a su lugar.
	//Paso 5
	//Rz_i(alpha) 
	
	//Paso 6
	//Ry_i(-theta)
	
	//Paso 7
	//T_i(x1,y1,z1)
	 
	//Matriz de rotacion compuesta.								
	arma:: dmat Comp = Rx* Rz * Ry * T;
	//arma::dmat Comp = T_i *Ry_i * Rz_i
	
	return(Comp);
		
}

int main(int argc, char *argv[]){

	arma::fcolvec v1 ={3,4,0,1};
	arma::fcolvec v2 ={7,10,3,1};
	arma::fcolvec p = {0,10,0,1};
	
	arma:: dmat Comp = rot_v1v2(v1,v2,90);
	
	arma::dcolvec v1p = Comp *v1;
	arma:: dcolvec v2p = Comp	*v2;		
	
	arma::dcolvec pp =Comp *p;
	
	cout<<v1p <<endl <<v2p <<endl << pp <<endl;
}


