#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

// define the dimensionless constants
const double c    = 13;      
const double b      = 1.0e12; // 1.0e-4; //1.0e6;  //100;
      double a     = pow(c,2.0)/100;        // used the relation M = 10 Mp

// derived constants
      double P = 1.0/b;
      double Q = a * P;
      double R = 1.0 / (a * b);
      double phi0 = 1;

double A(double phi){
      return 1.0 + a * pow(phi,2.0);
}
double V(double phi){
      return -0.5 * pow(phi,2.0) + 0.25 * pow(phi,4.0);
}
double dA(double phi){
      return phi;
}
double dV(double phi){
      return -1.0 * phi + pow(phi,3.0);
}
double den_chi(
      double phi,
      double chi,
      double dphi,
      double dchi){
  
      return 0.5 * pow(A(phi) * dchi,2.0) + 0.5 * b * pow(A(phi),4.0) * pow(chi,2.0);   
}
double den_phi(
      double phi,
      double chi,
      double dphi,
      double dchi){

      return 0.5 * pow(dphi,2.0) + V(phi);
}
double hubble(
      double phi,
      double chi,
      double dphi,
      double dchi){

      return c * pow(1.0/3.0 * ( -V(1.0) + den_phi(phi,chi,dphi,dchi) + R * den_chi(phi,chi,dphi,dchi)) ,0.5);   //   
}
double Veff(
      double phi,
      double chi,
      double dphi,
      double dchi){
      return V(phi) + R * A(phi) * den_chi(phi,chi,dphi,dchi);		  // in units of mu^4/lambda
}
// The driving force terms in the EOMs of phi and chi
double F_phi(
      double phi,
      double chi,
      double dphi,
      double dchi){

      return dV(phi) + P * dA(phi) / A(phi) * den_chi(phi,chi,dphi,dchi);
}
double F_chi(
      double phi,
      double chi,
      double dphi,
      double dchi){

      return b * ( pow(A(phi), 2.0) * chi + Q * 2.0 * dA(phi)/A(phi) * dphi *dchi );
}

double f1(
      double phi,
      double chi,
      double dphi,
      double dchi){
      return dphi;		  
}
double f2(
      double phi,
      double chi,
      double dphi,
      double dchi){
      return   dchi;	
}
double f3(
      double phi,
      double chi,
      double dphi,
      double dchi){
      return -3.0 * hubble(phi,chi,dphi,dchi) * dphi - F_phi(phi,chi,dphi,dchi);			  
}
double f4(
      double phi,
      double chi,
      double dphi,
      double dchi){
      return -3.0 * hubble(phi,chi,dphi,dchi) * dchi - F_chi(phi,chi,dphi,dchi);		  
}
double dden_phi(
      double phi,
      double chi,
      double dphi,
      double dchi){
      return -dphi * (3.0* dphi *hubble(phi,chi,dphi,dchi) + F_phi(phi,chi,dphi,dchi)) + dphi * dV(phi); 
}
double dden_chi(
      double phi,
      double chi,
      double dphi,
      double dchi){
      return 0.5 * dA(phi) * dphi * pow(dchi, 2.0) - A(phi) * dchi * (3.0 * dchi * hubble(phi,chi,dphi,dchi) + F_chi(phi,chi,dphi,dchi)) + 1.5 * b * dphi * pow(A(phi)*chi,2.0) + b * pow(A(phi),3.0) * chi * dchi; 
}
double dH(
      double phi,
      double chi,
      double dphi,
      double dchi){
      return 0.5 * pow(c,2.0) /(hubble(phi,chi,dphi,dchi)) * (1.0/3.0 * (dden_phi(phi,chi,dphi,dchi) + R * dden_chi(phi,chi,dphi,dchi)));
}
double epsilon(
      double phi,
      double chi,
      double dphi,
      double dchi){
      return -dH(phi,chi,dphi,dchi) * pow(hubble(phi,chi,dphi,dchi) , -2.0);
}
int main(){
	  /* Test the symmetron model for different choices of initial conditions for phi and chi fields */
      double phi    = 0.01; //0.818973; //0.0243837;   //0.0101014;  //0.00924935; //0.01; //0.0124481;  //0.00924943; //0.01; //0.0001; //0.01;
      double dphi   = 0;   //0.0687869; //0.00214899;  //0.000890267; //0.000814808; //0;    //0.00109709; //0.000814816; //0; //0.0001; //0.01;
      double chi    = 1.0e2; //1.0e6; //0;//(-9.88131e-324;) //1.09402e-05;  //0.00487858; //1.0e5; //-5.92879e-323; //0.00620518; //1.0e5; //1.0/sqrt(2.0) ;
      double dchi   = 0;    //0;//(-5.24751e-317;)  //-17.7571;    //-3854.41;   //0;     //3.89013e-318;  //-322.902 ; //0; //20/sqrt(b); //1/sqrt(2.0)/sqrt(b);
      
      double t      = 0; //51.9948; //11.9988;   //1.9998;   //0.9999;   //0; 
      double efolds = -60; //128.276;  //-7.53769;  //-45.0513; //-48.8034; //-60.0;
      double Hubble = 0;
      double den__phi = 0;
      double den__chi = 0;            // density of phi and chi fields in units of the vacuum energy
      double den_chi_mod = 0;
      double E_vac   = -V(1.0);
      double scale   = 0;            // scale factor
      double Epsilon = 0;
      double h     = 0.00000001;
      double hh    = 10.0 * h;   //40*h;   //10 * h;
      int    nsteps = 1000000000;
      double k1,k2,k3,k4,j1,j2,j3,j4,l1,l2,l3,l4,m1,m2,m3,m4;
      int i,j;
      ofstream out("output.txt");

      cout <<"t \t a \t phi \t chi \t den-phi \t den-chi \t Hubble \t Epsilon \n";
     
      for(j = 1; j < nsteps; j++){
          t       += hh;
          den__phi = den_phi(phi,chi,dphi,dchi) / E_vac;
          den__chi = R * den_chi(phi,chi,dphi,dchi) / E_vac;
          den_chi_mod = den__chi / A(phi);
          Hubble = hubble(phi,chi,dphi,dchi);
          efolds += hh * Hubble;
          scale = exp(efolds);
          Epsilon = epsilon(phi,chi,dphi,dchi);
          out<< t <<'\t'<< scale <<'\t'<< phi <<'\t'<< chi <<'\t'<< den__phi <<'\t'<< den__chi <<'\t'<< den_chi_mod <<'\t'<< Hubble <<'\t' << Epsilon<<'\n';
        
          /* Solve the coupled equations using the 4th-order Runge Kutta method */
          k1 = f1(phi,chi,dphi,dchi);
          j1 = f2(phi,chi,dphi,dchi);
          l1 = f3(phi,chi,dphi,dchi);
          m1 = f4(phi,chi,dphi,dchi);
          k2 = f1(phi + hh/2.0 * k1, chi + hh/2.0 * j1, dphi + hh/2.0 * l1, dchi + hh/2.0 * m1);
          j2 = f2(phi + hh/2.0 * k1, chi + hh/2.0 * j1, dphi + hh/2.0 * l1, dchi + hh/2.0 * m1);
          l2 = f3(phi + hh/2.0 * k1, chi + hh/2.0 * j1, dphi + hh/2.0 * l1, dchi + hh/2.0 * m1);
          m2 = f4(phi + hh/2.0 * k1, chi + hh/2.0 * j1, dphi + hh/2.0 * l1, dchi + hh/2.0 * m1);
          k3 = f1(phi + hh/2.0 * k2, chi + hh/2.0 * j2, dphi + hh/2.0 * l2, dchi + hh/2.0 * m2);
          j3 = f2(phi + hh/2.0 * k2, chi + hh/2.0 * j2, dphi + hh/2.0 * l2, dchi + hh/2.0 * m2);
          l3 = f3(phi + hh/2.0 * k2, chi + hh/2.0 * j2, dphi + hh/2.0 * l2, dchi + hh/2.0 * m2);
          m3 = f4(phi + hh/2.0 * k2, chi + hh/2.0 * j2, dphi + hh/2.0 * l2, dchi + hh/2.0 * m2);
          k4 = f1(phi + hh * k3, chi + hh * j3, dphi + hh * l3, dchi + hh * m3);
          j4 = f2(phi + hh * k3, chi + hh * j3, dphi + hh * l3, dchi + hh * m3);
          l4 = f3(phi + hh * k3, chi + hh * j3, dphi + hh * l3, dchi + hh * m3);
          m4 = f4(phi + hh * k3, chi + hh * j3, dphi + hh * l3, dchi + hh * m3);
          phi = phi + hh/6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
          chi = chi + hh/6.0 * (j1 + 2.0 * j2 + 2.0 * j3 + j4);
          dphi = dphi + hh/6.0 * (l1 + 2.0 * l2 + 2.0 * l3 + l4);
          dchi = dchi + hh/6.0 * (m1 + 2.0 * m2 + 2.0 * m3 + m4);

          if(j%1000 == 0){
			 cout<< t <<'\t'<< phi <<'\t'<< chi <<'\t'<< dphi <<'\t'<< dchi <<'\t'<< efolds <<'\n';
		  }
      }
      
      
}
