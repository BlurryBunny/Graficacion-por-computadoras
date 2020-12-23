#include <iostream>
#include <armadillo>
#include <cmath>

#define PI 3.1416

using namespace std;
arma::dmat rot_v1v2v3(arma::dcolvec v1,arma::dcolvec v2,arma::dcolvec v3){
	
	//Step 1
	//Translate P1 to origin
	
	//matrix to translate P1 and every point according to P1
	arma::dmat T = {
									{1,0,0,-v1[0]},
									{0,1,0,-v1[1]},
									{0,0,1,-v1[2]},
									{0,0,0,1}
									};

	//Apply T to each point 
	//arma::dcolvec p1_p = T * v1;
	//arma::dcolvec p2_p = T * v2;
	//arma::dcolvec p3_p = T * v3;

	//Step 2
	//Angle rotation is -(90 - alpha) = alpha - 90 which is equal to ((z2-z1)/D1)
	
	//D1 = sqrt(pow(z2_p,2)+pow(x2_p,2)
	double D1= sqrt(powf(v2[0]-v1[0],2) + powf(v2[2]-v1[2],2));
	
	//We need to rotate along the y axis 
	arma::dmat Ry = {
									{(v2[0]-v1[0])/D1,0,(v2[2]-v1[2])/D1,0},
									{0,1,0,0},
									{-((v2[2]-v1[2])/D1),0,(v2[0]-v1[0])/D1,0},
									{0,0,0,1}
									};
									
	//p2_pp = Ry (alpha - 90)
	//p2_pp = {0, y2-y1, D1, 1}^T.
	//arma::dcolvec p2_tosquare_1={0,v2[1]-v1[1],D1,1};
	//arma::dcolvec p2_pp = sqrt(p2_tosqueare_1,T);
	//arma::dcolvec p2_pp = p2_tosquare_1 * Ry;
	
	
	//Step 3
	//Rotation in x axis
	//D2 = |p1_pp*p2_pp|
	double D2 = sqrt(powf(v2[0]-v1[0],2) + powf(v2[1]-v1[1],2) + powf(v2[2]-v1[2],2));
	
		arma::dmat Rx = {
									{1,0,0,0},
									{0,(v2[2]-v1[2])/D2,-(v2[1]-v1[1])/D2,0},
									{0,(v2[1]-v1[1])/D2,(v2[2]-v1[2])/D2,0},
									{0,0,0,1}
									};

	//p2_ppp = {0, 0, |p1*p2|, 1}^T.
	//arma::dcolvec p2_tosquare_2={0,0,D2,1};
	//arma::dcolvec p2_ppp = sqrt(p2_tosqueare_2,T);
	//arma::dcolvec p2_ppp = p2_tosquare_2 * Rx;
	
	//Step 4
	//Rotate about the z axis
	//cos(lamda)= y3_ppp/D3			sin= x3_ppp/D3
	
	//D3 = sqrt(powf(x3_ppp,2) + powf(y3_ppp,2))
	double D3 = sqrt(powf(v3[0]-v1[0],2) + powf(v3[1]-v1[1],2));
	
		arma::dmat Rz = {
									{(v3[1]-v1[1])/D3,(v3[0]-v1[0])/D3,0,0},
									{-(v3[0]-v1[0])/D3,(v3[1]-v1[1])/D3,0,0},
									{0,0,1,0},
									{0,0,0,1}
									};

	 
	//Matrix rotation complete.							
	arma:: dmat Comp = Rz* Rx * Ry * T;
	
	return(Comp);
		
}

int main(int argc, char *argv[]){

	//Three points in 3D view
	arma::dcolvec v1 ={3,4,3,1};
	arma::dcolvec v2 ={7,10,7,1};
	arma::dcolvec v3 ={5,14,3,1};
	
	arma::fcolvec p = {0,10,0,1};
	
	//We obtain the Matrix Comp from v1, v2, v3 to rotate it.
	arma:: dmat Comp = rot_v1v2v3(v1,v2,v3);
	
	//p1 to origin 
	arma::dcolvec v1p = Comp *v1;
	//p2 to positive z axe
	arma:: dcolvec v2p = Comp	*v2;
	//p3 to positive y axe		
	arma:: dcolvec v3p = Comp	*v3;
		
				
	arma::dcolvec pp =Comp *p;
	
	cout<<"V1 PPP: \n"<<v1p <<endl<<"V2 PPP: \n" <<v2p <<endl<<"V3 PPP: \n" << v3 <<"Pp: \n"<<pp<<endl;
}


