#include <iostream>
#include <armadillo>
#include <cmath>
#include<time.h>
#include <vector>

#define PI 3.1416

using namespace std;
arma::dmat rot_v1v2(arma::dcolvec v1,arma::dcolvec v2, double angle){
	
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

	//Matriz de rotacion compuesta.								
	arma:: dmat Comp = Rx* Rz * Ry * T;
	//arma::dmat Comp = T_i *Ry_i * Rz_i
	
	return(Comp);
		
}

int main(int argc, char *argv[]){
	
	srand(time(NULL));
	
	double theta = rand()%360+1;
	std::vector <arma::dcolvec> v;
	
	arma::fcolvec aux;
	
	//puntos que forman la diagonal
 	 v.push_back(aux={0,0,0,1});
	 v.push_back(aux={1,1,1,1});
	 v.push_back(aux={0,0,1,1});
	 v.push_back(aux={1,0,0,1}); 
	 v.push_back(aux={1,0,1,1});
	 v.push_back(aux={0,1,0,1});
	 v.push_back(aux={1,1,0,1});
	 v.push_back(aux={0,1,1,1});
	
	
	arma:: dmat Comp = rot_v1v2(v[0],v[1],theta);
	
	std::vector <arma::dcolvec> vp;
	for(int i=0;i<v.size();i++){
		vp.push_back(v[i]*Comp);
	}
	
	cout<<"Puntos del cubo unitario con rotacion"<<endl;
	cout<<"Angulo elegido: "<<theta<<endl;
	for(int i=0;i<v.size();i++){
		cout<<vp[i]<<endl;
	}
	
}


