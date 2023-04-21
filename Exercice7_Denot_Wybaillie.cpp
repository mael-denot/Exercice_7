#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string.h>
#include "ConfigFile.tpp"
#include <algorithm>

using namespace std;

double norme(vector<double> const& v)
{
  double norme_(0.);
  for(unsigned int i(0); i<v.size(); ++i)
    norme_ += v[i]*v[i];
  return sqrt(norme_);
}

double norme_sqrt(vector<double> const& v)
{
  double norme_(0.);
  for(unsigned int i(0); i<v.size(); ++i)
    norme_ += v[i];
  return sqrt(norme_);
}

void boundary_condition(vector<double> &fnext, vector<double> &fnow, double const& A, double const& omega,\
		double const& t,double const& dt, \
		vector<double> &beta2, string &bc_l, string &bc_r, int &N)
{
      if (bc_l == "fixe"){
        fnext[0] = fnow[0]; // done : Completer la condition au bord gauche fixe
      }else if(bc_l == "libre"){
        fnext[0] = fnext[1]; // done : Completer la condition au bord gauche libre
      }else if(bc_l== "harmonique"){
        fnext[0] = A*sin(omega*(t+dt)); // done : Completer la condition au bord gauche harmonique
      }else if (bc_l =="sortie"){
        fnext[0] = fnow[0] -  norme_sqrt(beta2)*(fnow[0] - fnow[1]) ; // done : Completer la condition au bord gauche "sortie de l'onde"
      }else{
        cerr << "Merci de choisir une condition valide au bord gauche" << endl;
      }
	      
      if (bc_r == "fixe"){
        fnext[N-1] = fnow[N-1]; // done : Completer la condition au bord droit fixe
      }else if(bc_r == "libre"){
        fnext[N-1] = fnext[N-2]; // done : Completer la condition au bord droit libre
      }else if(bc_r== "harmonique"){
        fnext[N-1] = A*sin(omega*(t+dt)); // done : Completer la condition au bord droit harmonique
      }else if (bc_r =="sortie"){
        fnext[N-1] = fnow[N-1] - sqrt(beta2[N-1])*(fnow[N-1] - fnow[N-2]); // done : Completer la condition au bord droit "sortie de l'onde"
      }else{
        cerr << "Merci de choisir une condition valide au bord droit" << endl;
      }
      cout << "beta2 = " << norme_sqrt(beta2) << endl;
}

//
// TODO : Calcul de l'energie de l'onde
//
double energie(vector<double> const& f, double const& dx)
{
  double energie_(0.);
  for(unsigned int i(0); i<f.size()-1; ++i)
    energie_ += f[i]*f[i]*dx;

  return energie_;
}

//
// Surcharge de l'operateur pour ecrire les elements d'un tableau
//
template <class T> ostream& operator<< (ostream& o, vector<T> const& v)
{
  unsigned int len(v.size());
  for(unsigned int i(0); i < (len - 1); ++i)
    o << v[i] << " ";
  if(len > 0)
    o << v[len-1];
  return o;
}

double E(vector<double> const& f, double const& dx) {
  double E(0.);
  for (int i = 0; i < f.size(); i++) {
    E += f[i] * f[i] * dx;
  }
  return E;
}

//
// Main
//
int main(int argc, char* argv[])
{
  const double PI = 3.1415926535897932384626433832795028841971e0;
  const double g  = 9.81;
  double dx;
  double dt;
  double t;
  int stride(0);

  string inputPath("configuration.in"); // Fichier d'input par defaut
  if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice7 config_perso.in")
    inputPath = argv[1];

  ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

  for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice7 config_perso.in input_scan=[valeur]")
    configFile.process(argv[i]);

  // Parametres de simulation :
  double tfin    = configFile.get<double>("tfin");
  int N          = configFile.get<int>("Npoints");
  double CFL     = configFile.get<double>("CFL");
  double A       = configFile.get<double>("A");
  double fmn     = configFile.get<double>("fmn");
  double minit   = configFile.get<double>("minit");
  double omega   = configFile.get<double>("omega");
  double hL      = configFile.get<double>("hL");
  double hR      = configFile.get<double>("hR");
  double h00     = configFile.get<double>("h00");
  double xa      = configFile.get<double>("xa");
  double xb      = configFile.get<double>("xb");
  double xL      = configFile.get<double>("xL");
  double xR      = configFile.get<double>("xR");
  int n_stride(configFile.get<int>("n_stride"));
  string bc_l           = configFile.get<string>("cb_gauche");
  string bc_r           = configFile.get<string>("cb_droit");
  string initialization = configFile.get<string>("initialization");
  bool v_uniform        = configFile.get<bool>("v_uniform");

  vector<double> h0(N) ;
  vector<double> vel2(N) ;
  vector<double> x(N) ;
  vector<double> fpast(N), fnow(N), fnext(N), beta2(N);

  dx = (xR - xL) / (N-1);
  bool ecrire_f = configFile.get<bool>("ecrire_f"); // Exporter f(x,t) ou non
  string schema = configFile.get<string>("schema");
  

  for(int i(0); i<N; ++i){ 
     x[i] = xL + i * dx ;
     h0[i] = 0.0;
     if(v_uniform){
           h0[i]  = h00;
     } else {
           h0[i]  = hL * (xL<=x[i] && x[i]<=xa) + \
		   (0.5*(hL + hR) + 0.5*(hL - hR)*cos(PI *(x[i]-xa)/(xa - xb))) * (xa<x[i] && x[i]<xb) +\
		   hR * (xb<=x[i] && x[i]<=xR);
	   if(i==0) cout << "v is not uniform"<<endl;
     }
     vel2[i]  = g* h0[i];
  }



  auto max_vel2 = std::max_element(vel2.begin(), vel2.end());
  dt = CFL * dx / sqrt(*max_vel2);
  cout << "dt is "<< dt<< " " <<dx<< endl;

  // Fichiers de sortie :
  string output = configFile.get<string>("output");

  ofstream fichier_x((output + "_x.out").c_str());
  fichier_x.precision(15);

  ofstream fichier_v((output + "_v.out").c_str());
  fichier_v.precision(15);

  ofstream fichier_f((output + "_f.out").c_str());
  fichier_f.precision(15);

  ofstream fichier_E((output + "_E.out").c_str());
  fichier_E.precision(15);



  // Initialisation des tableaux du schema numerique :

  //TODO initialize f and beta
  for(int i(0); i<N; ++i)
  {
    fpast[i] = 0.;
    fnow[i]  = 0.;
    beta2[i] = vel2[i]*dt*dt/(dx*dx);
    if(initialization=="cos"){
	fnow[i] = A*cos(2*PI*fmn*x[i]);
	fpast[i] = A*cos(2*PI*fmn*x[i] - omega*dt);
    }else if(initialization=="sin"){
	fnow[i] = 0.;
	fpast[i] = 0.;
    }
  }
  cout<<"beta2[0] is "<<beta2[0]<<endl;


  // Boucle temporelle :
  for(t=0.; t<tfin-.5*dt; t+=dt)
  {
    // Ecriture :
    if(stride%n_stride == 0)
    {
      if(ecrire_f) fichier_f << t << " " << fnow << endl;
      fichier_E << t << " " << energie(fnow,dx) << endl;
    }
    ++stride;


    // Evolution :
    for(int i(1); i<N-1; ++i)
    {
    fnext[i] = beta2[i]*(fnow[i+1] - 2*fnow[i] + fnow[i-1]) - (fpast[i] - 2*fnow[i]); // Done : Compléter le schéma A
      if(schema == "B"){
        //fnext[i] = dt*dt*(vel2[i+1]*(fnow[i+1] - fnow[i])/dx*dx - vel2[i]*(fnow[i] - fnow[i-1])/dx*dx) - (fpast[i] - 2*fnow[i]); // Done : Compléter le schéma B (en ajoutant des termes au A)
        fnext[i] = beta2[i+1]*(fnow[i+1] - fnow[i]) - beta2[i]*(fnow[i] - fnow[i-1]) - (fpast[i] - 2*fnow[i]);
      }else if(schema == "C"){
        fnext[i] = beta2[i+1]*fnow[i+1] - 2*beta2[i]*fnow[i] + beta2[i-1]*fnow[i-1] - (fpast[i] - 2*fnow[i]); // Done : Compléter le schéma C (en ajoutant des termes au A)
      }
    }



    // add boundary conditions
    boundary_condition(fnext, fnow, A, omega, t, dt, beta2, bc_l, bc_r, N);

    // Mise a jour :
    fpast = fnow;
    fnow  = fnext;
  }


  if(ecrire_f) fichier_f << t << " " << fnow << endl;
  fichier_E << t << " " << energie(fnow,dx) << endl;
  fichier_x << x << endl;
  fichier_v << vel2 << endl;

  fichier_f.close();
  fichier_E.close();
  fichier_x.close();
  fichier_v.close();

  cout << "Fin de la simulation." << endl;

  return 0;
}
