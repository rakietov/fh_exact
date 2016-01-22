#include <iostream>
#include <armadillo>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <complex>

#include "basis.h"
#include "hamiltonian.h"
#include "level_stat.h"

using namespace std;
using namespace arma;




int main()
{
    vector<int> vL = {20}; // for fermions: twice the physical lattice size
    vector<int> vN = {5};
    //vector<double> vdelta = {0.0001,0.001,0.01,0.1, 0.2, 0.3,0.5, 0.8,1.0, 1.3, 1.8, 2.3, 2.8, 5.,6.,10.,15.};
    vector<double> vdelta = {0.1};
    vector<double>  vW = {0.};
 
	//vector<double> vWU = {15};
   vector<double>  vWU = {45., 50.,55., 60.,65., 70.,75., 80.,85., 90.,95., 100., 105., 110., 115., 120.};
 


    int L ;
    int N ;
    int zmax = 300;
    int n_of_eig = 50; 
    int first_eig = 12; 
    int last_eig  = 37;
    double delta;
    double W;
    double WU;
    
    
    
    
    //cin >> L >> N >> W >> WU >> zmax >> n_of_eig >> first_eig >> last_eig >> delta;
    
 
   srand (time(NULL));

    
    for(int il = 0; il < vL.size(); ++il)
    {
        for(int in = 0; in < vN.size(); ++in)    
        {
            Tbasis basis1( vL[il], vN[in], 3, 2);

           
            for(int iW = 0; iW < vW.size(); ++iW) 
            {
                for(int iWU = 0; iWU < vWU.size(); ++iWU) 
                {                        
                    string fname = "L_";
                    stringstream ssL, ssN,ssW, ssWU, sszmax;
                    
                    ssL << vL[il];
                    ssN << vN[in];
                    ssW << vW[iW];
                    ssWU << vWU[iWU];

                    
                    fname += ssL.str();
                    fname += string("N_");
                    fname += ssN.str();
                    fname += string("W_");
                    fname += ssW.str();
                    fname += string("WU_");
                    fname += ssWU.str();  
    
                    basis1.write_vectors( fname );
                    
                    for(int z = 0; z < zmax; ++z)
                    {
                    
                        double t = -1.0;
                        vector<double> U = {};
                        vector<double> eps = {};
                        
                        for( int i = 0; i < vL[il]; ++i)
                        {
		                double ui = vWU[iWU] *( rand()/double(RAND_MAX) );

                            U.push_back(ui);
                            // note: for fermions only U[2i], U[2i+1] matter
							//		 and signify interaction on i'th physical site.
							//       if e_i's are nonzero, eps[2i] = eps[2i+1]  
							//		 must be imposed

                            double e_i = vW[iW]*(rand()/double(RAND_MAX) - 0.5);

							if( i % 2 == 1)e_i == eps[i-1];
                            //cout << "e_i: "<<e_i <<endl;            
                            //double e_i = (i-5.)*0.01;
                            eps.push_back(e_i);
                        }
		        for(int i = 0; i < vL[il]; ++i)
			{
				std::cout << U[i] << " ";
			}
			std::cout << std::endl;


						int n_of_eig = basis1.get_n_vec();
                        Thamiltonian ham1( basis1, t, U, eps, n_of_eig );
                        //ham1.is_symmetric();
                        ham1.diagonalize();
                        //cout << "write_eigenvalues \n";
                        ham1.write_eigenvalues( fname );
                        
                        for(int idel = 0; idel < vdelta.size(); ++idel)
                        { 
                            delta = vdelta[idel];
                            Tlevel_stat lvstat1( ham1.get_eigenvalues(), n_of_eig/4., 3.*n_of_eig/4., delta );

                            string tmpstr = fname;
                            
                            stringstream ssdelta;
                            ssdelta << delta;

                            tmpstr += string("d_");
                            tmpstr += ssdelta.str();
                            
                           // lvstat1.write_eigenvalues();
                           
                           // lvstat1.unfold_spectrum();
                           // lvstat1.write_unfolded_spectrum();
                            //cout << "write_spacings \n";
                            
                            //lvstat1.write_unfolded_level_spacings( tmpstr );
                            lvstat1.write_avRBSL( tmpstr );
                        }
                    }
                }
            }
        }
    }


    return 0;
}


























