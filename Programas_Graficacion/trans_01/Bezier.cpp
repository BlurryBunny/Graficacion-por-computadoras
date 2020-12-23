#include <iostream>
#include <armadillo>
#include <cmath>

using namespace std;
int main(int argc, char *argv[]){


	//columnas de nuestra nueva matriz
	arma::fcolvec v1 ={3,4,0,1};
	arma::fcolvec v2 ={7,10,3,1};
	
	
	//Rotacion en X
	arma::fmat T = {{1,0,0,-3},
									{0,1,0,-4},
									{0,0,1,0},
									{0,0,0,1}};

	arma:: fcolvec v1p = T *v1;
	arma:: fcolvec v2p = T *v2;

	cout<< "Matriz origen :";
	cout<< "\nV1:\n "<< v1 << "\nV2:\n "<< v2 << endl;
	
	cout<< "Matriz prima :";
	cout<< "\nV1p:\n "<< v1p << "\nV2p:\n "<< v2p << endl;
	
	//Rotacion en Y
	arma::fmat Ry = {{0.80,0,0.60,0},
									{0,1,0,0},
									{-0.60,0,0.80,0},
									{0,0,0,1}};

	arma:: fcolvec v1pp = Ry *v1p;
	arma:: fcolvec v2pp = Ry *v2p;
	
	cout<< "Matriz Bi-Prima :";
	cout<< "\nV1pp:\n "<< v1pp << "\nV2pp:\n "<< v2pp << endl;
	
	//Rotacion en Z
	arma::fmat Rz = {{0.64,0.77,0,0},
									{-0.77,0.64,0,0},
									{0,0,1,0},
									{0,0,0,1}};
									
	arma:: fcolvec v1ppp = Rz *v1pp;
	arma:: fcolvec v2ppp = Rz *v2pp;								
	
	cout<< "Matriz BTri-Prima :";
	cout<< "\nV1ppp:\n "<< v1ppp << "\nV2ppp:\n "<< v2ppp << endl;			
	
	cout<< "Matrix resultado";
	cout << v1ppp << endl << v2ppp << endl;
	
	
	//Curvas de rotacion.  Buscando el angulo de rotacion.
	
	//Vertices 
	arma::fcolvec v1_p ={3,4,0};
	arma::fcolvec v2_p ={6,5,0};
	
	//Paso 1
	//T(-3,-4,0)
	
	cout<<"\nCalculo con angulos"<<endl;
	arma::fmat T2 = {{1,0,0,-3},
									{0,1,0,-4},
									{0,0,1,0},
									{0,0,0,1}};
	
	//Paso 2
	double theta = abs (atan (3.0/4.0));
	//Ry(tan -1(3/4)= + 36.87)
	arma::dmat Ry2 = {{cos(theta),sin(theta),0,0},
									{0,1,0,0},
									{-sin(theta),0,cos(theta),0},
									{0,0,0,1}};
	//Paso 3
	double alpha = abs (atan (6.0/5.0));
	//Ry(tan -1(3/4)= - 36.87)
	arma::dmat Rz2 = {{cos(-alpha),-sin(-alpha),0,0},
									{sin(-alpha),cos(-alpha),0,0},
									{0,0,1,0},
									{0,0,0,1}};
	
	arma :: dmat Comp = Rz2 *Ry2 * T2;
	arma:: dcolvec v1p_2 = Comp *v1_p;
	arma:: dcolvec v2p_2 = Comp *v2_p;
	
	cout<< v1p_2 << endl << v2p_2 << endl;
	
	
	
				
}
