/*
 * Try to do it in real space;
 */
#include <iostream>
#include <string>
#include <math.h>
#include <cmath>
#include <fstream>
#include <vector>
#include <complex>
#include <algorithm>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#define PI 3.141592653589
using namespace std;

void get_potential(std::string filename, \
                        int &nx,
                        int &ny,
                        int &nz,
                        double (&cell)[3],
                        std::vector<std::vector<std::vector<double> > > &cube);
void readfile(std::string filename, \
              int &nx, \
              int &ny, \
              int &nz, \
              double (&cell)[3], \
              std::vector<double> &a, \
              std::vector<std::vector<std::vector<double> > > &cube );
void MatxVec(std::vector<double> &Res, \
             std::vector<std::vector<double> > Mat, \
             std::vector<double> Vec, unsigned N);

typedef struct {
    string element;
    double x;
    double y;
    double z;
} atom;

int main ( int argc, char **argv )
{
    std::string tmp;
    int nx=0, ny=0, nz=0;
    double cell[3], tmpdbl;
    int kcut;
    double vdW_Si = 3.24277, vdW_N = 2.63995, vdW[2] = {vdW_Si, vdW_N}; // van der Waals radius of silicon atom, in Bohr.
    std::string element[2] = {"si", "n"};
    if ( argc == 1) { cout << "No input file" << endl; } 
    std::string name, system=argv[1];
    std::vector< std::vector< std::vector<double> > > cube;    // Original data are reshaped into a (nx, ny, nz) 3D array
    double dV, ddV;  // dV: first order derivative, ddV: second order derivative
    std::vector<double> coeff(3), Vec(3); // Coefficients calculated. 
    std::vector< std::vector<double> > Mat = {{10,-4,0.5},{-15,7,-1},{6,-3,0.5}};

    /*
     *  We will read solid calculation from file and combine it with single atom results.
     */
    int NX, NY, NZ; // Capital variables are used for parameters of solid system. 
    double CELL[3];
    std::vector<std::vector<std::vector<double> > > CUBE;
    get_potential(system.append(".dat"), NX, NY, NZ, CELL, CUBE);
    double DX=CELL[0]/NX, DY=CELL[1]/NY, DZ=CELL[2]/NZ;
    cout << DX << "  " << DY << "  " << DZ <<endl;
    std::vector<atom> solid;
    atom tmpat;
    system = argv[1];
    ifstream xyzFile(system.append("-bohr.xyz"));
    std::getline(xyzFile, tmp);
    std::getline(xyzFile, tmp);
    while ( xyzFile >> tmpat.element >> tmpat.x >> tmpat.y >> tmpat.z ) {
        solid.push_back(tmpat);
    }
    xyzFile.close();
    int NAT = solid.size(); // Number of atoms in one cell of solid;
    // 暴力扩展晶胞，在原有晶胞周围添加原子，范围由vdW确定；
    // Add mirror atoms out of the cell but within a shell whose length is vdW
    int ie; // index of element, 0 for si, 1 for n
    for ( int n = 0; n < NAT; n++ ) {
        for ( int i = -1; i < 2; i++ ) { // i = -1, 0, 1, along x direction
            for ( int j = -1; j < 2; j++ ) { // j = -1, 0, 1, along y direction
                for ( int k = -1; k < 2; k++ ) { // k = -1, 0, 1, along z direction
                    if ( i*i+j*j+k*k > 0 ) {
                        tmpat.element = solid[n].element;
                        tmpat.x = solid[n].x - CELL[0]*i;
                        tmpat.y = solid[n].y - CELL[1]*j;
                        tmpat.z = solid[n].z - CELL[2]*k;
                        if ( tmpat.element == "si" ) { ie = 0; }
                        else if ( tmpat.element == "n" ) { ie = 1;}
                        if ( tmpat.x > 0-vdW[ie] && tmpat.x < CELL[0]+vdW[ie] \
                          && tmpat.y > 0-vdW[ie] && tmpat.y < CELL[1]+vdW[ie] \
                          && tmpat.z > 0-vdW[ie] && tmpat.z < CELL[2]+vdW[ie] ) {
                            solid.push_back(tmpat);
                        }
                    }
                }
            }
        }
    }
    cout << solid.size() << endl;                
    system = argv[1];
    ofstream bigint(system.append("-big-angstrom.xyz"));
    bigint << solid.size() << endl;
    bigint << "expanded solid" << endl;
    for (unsigned int n = 0; n < solid.size(); n++ ) {
        bigint << solid[n].element\
               << "  " << solid[n].x*0.529177 \
               << "  " << solid[n].y*0.529177 \
               << "  " << solid[n].z*0.529177  << endl;
    }
    bigint.close();
    
    std::vector<std::vector<std::vector<double> > > potential(NX, std::vector<std::vector<double> > \
                                                              (NY, std::vector<double>(NZ,0)));
    /*
     * In the following for loop, we read single atom calculation from file
     * and get radius distribution of electrostatic potential.
     */
    for (unsigned int ie =0; ie < sizeof(element)/sizeof(std::string); ie++ ) {
        cube.erase(cube.begin(), cube.end());
        name = element[ie];
        cout << "We are dealing with " << name << " ..." << endl;
        get_potential(name.append("-atom.dat"), nx, ny, nz, cell, cube );
        // pz means potential along z-direction
        double *pz = new double[nz];
        double *bigpz = new double[nz*10];
        double *newpz = new double[nz*10];
        name = element[ie];
        ofstream original_pz(name.append("-original-pz.dat"));
        for (int  k = 0; k < nz; k++ ) {
            // tmpdbl = 0.0;
            // for ( int j = 0; j < ny; j++ ) {
            //     for ( int ii = 0; ii < nx; ii++ ) {
            //         tmpdbl += cube[ii][j][k];
            //     }
            // }
            // pz[k] = tmpdbl/(nx*ny);
            pz[k] = cube[0][0][k];
        }

        for ( int k = 0; k < nz; k++ ) {
            original_pz << cell[2]/nz*k << "  " << pz[k] << endl;
        }
        original_pz.close();

        double *z = new double[nz];            // coordinate along z direction;
        double *bigz = new double[nz*10];
        for ( int k = 0; k < nz; k++ ){
            z[k] = cell[2]/nz*k;
        }
        for ( int k = 0; k < nz*10; k++ ) {
            bigz[k] = cell[2]/(nz*10)*k;
            if ( k == 0 ) { continue; }
            else if ( bigz[k] > vdW[ie] && bigz[k-1] < vdW[ie] ) {
                kcut = k;
            }
        }

        gsl_interp_accel *acc = gsl_interp_accel_alloc();
        gsl_spline *spline_cubic = gsl_spline_alloc(gsl_interp_cspline, nz+1);
        double *tmpz = new double[nz+1];
        double *tmppz = new double[nz+1];
        for ( int k = 0; k < nz+1; k++ ) {
            if ( k <= nz-1 ) { tmpz[k] = z[k]; tmppz[k] = pz[k]; }
            else { tmpz[k] = tmpz[nz-1] + cell[2]/nz; tmppz[k] = pz[0]; }
        }
        gsl_spline_init(spline_cubic, tmpz, tmppz, nz+1);
        for ( int k = 0; k < nz*10; k++ ) {
            bigpz[k] = gsl_spline_eval(spline_cubic, bigz[k], acc);
        }
        dV = gsl_spline_eval_deriv(spline_cubic, bigz[kcut], acc);
        ddV = gsl_spline_eval_deriv2(spline_cubic, bigz[kcut], acc);
        Vec = {bigpz[kcut],\
               (bigz[kcut])*dV,\
               pow(bigz[kcut],2)*ddV};
        MatxVec(coeff, Mat, Vec, 3);
        cout << coeff[0] << "  " << coeff[1] << "  " << coeff[2] << endl;
        gsl_spline_free(spline_cubic);
        gsl_interp_accel_free(acc);

        name = element[ie];
        ofstream newpzout(name.append("-dis-pz.dat"));
        for ( int k = 0; k < nz*10; k++ ) {
            if ( k > kcut ) {
                newpz[k] = 0.0;
            } else {
                tmpdbl = bigz[k]/bigz[kcut];
                newpz[k] = bigpz[k] - pow(tmpdbl,3)*(coeff[0] + coeff[1]*tmpdbl + coeff[2]*pow(tmpdbl,2));
            }
            newpzout << bigz[k] << "  " << newpz[k] << endl;
        }
        newpzout.close();
        
        gsl_interp_accel *acc2 = gsl_interp_accel_alloc();
        gsl_spline *spline_cubic2 = gsl_spline_alloc(gsl_interp_cspline, nz*10);
        gsl_spline_init(spline_cubic2, bigz, newpz, nz*10);
        double dis;
        for ( int i = 0; i < NX; i++ ) {
            for ( int j = 0; j < NY; j++ ) {
                for ( int k = 0; k < NZ; k++ ) {
                    for (unsigned int n = 0; n < solid.size(); n++ ) {
                        tmpat = solid[n];
                        if ( tmpat.element == element[ie] ) {
                            dis = sqrt(pow((i*DX-tmpat.x),2) \
                                      +pow((j*DY-tmpat.y),2) \
                                      +pow((k*DZ-tmpat.z),2));
                            if ( dis > vdW[ie] ) { potential[i][j][k] += 0.0; }
                            else { potential[i][j][k] += gsl_spline_eval(spline_cubic2, dis, acc2); }
                        }
                    }
                }
            }
        }

        gsl_spline_free(spline_cubic2);
        gsl_interp_accel_free(acc2);
        delete[] tmpz;
        delete[] tmppz;
    }

    system = argv[1];
    ofstream cubeout(system.append("_CUBE.dat"));
    system = argv[1];
    ofstream tot_p(system.append("_total_potential.dat"));
    for ( int i = 0; i < NX; i++ ) {
        for ( int j = 0; j < NY; j++ ) {
            tot_p << i*DX << "  " << j*DY << "  " << potential[i][j][96] << endl;
            cubeout << i*DX << "  " << j*DY << "  " << CUBE[i][j][96] << endl;
        }
    }
    cubeout.close();
    tot_p.close();

    std::vector<std::vector<std::vector<double> > > diff(NX, std::vector<std::vector<double> > \
                                                        (NY, std::vector<double> (NZ,0.0))); // diff = CUBE - potential;
    std::vector<double> ave_diff(NZ, 0.0), ave_CUBE(NZ, 0.0);
    double tmpdbl2;
    system = argv[1];
    ofstream ave_diff_out(system.append("_ave_diff.dat"));
    system = argv[1];
    ofstream ave_CUBE_out(system.append("_ave_CUBE.dat"));
    for ( int k = 0; k < NZ; k++ ) {
        tmpdbl = 0.0;
        tmpdbl2 = 0.0;
        for ( int j = 0; j < NY; j++ ) {
            for ( int i = 0; i < NX; i++ ) {
                diff[i][j][k] = CUBE[i][j][k] - potential[i][j][k];
                tmpdbl += diff[i][j][k];
                tmpdbl2 += CUBE[i][j][k];
            }
        }
        ave_diff[k] = tmpdbl/(NX*NY);
        ave_CUBE[k] = tmpdbl2/(NX*NY);
        ave_diff_out << k*DZ << "  " << ave_diff[k] << endl;
        ave_CUBE_out << k*DZ << "  " << ave_CUBE[k] << endl;
    }




    return 0;
}

void get_potential(std::string filename, \
                   int &nx,\
                   int &ny,\
                   int &nz,\
                   double (&cell)[3],\
                   std::vector<std::vector<std::vector<double> > > &cube)
{
    std::vector<double> row;                        // A temperary vector;
    std::vector< std::vector< double > > mat;       // A temperary matrix;
    std::vector<double> a;                          // 1D array to store original data linearly, dimension ( nx x ny x nz )
    cube.erase(cube.begin(), cube.end());
    readfile( filename, nx, ny, nz, cell, a, cube);  // Read data from file and put them into 1D array a,
                                                    // then a has been reshaped into a 3D array cube
                                                    //
}

void readfile(std::string                                      filename, \
              int                                              &nx,       \
              int                                              &ny,       \
              int                                              &nz,       \
              double                                           (&cell)[3],     \
              std::vector<double>                              &a,        \
              std::vector<std::vector<std::vector<double> > >  &cube      ) 
{
    std::string tmp;
    int nxm, nym, nzm, nat, ntyp, ibrav, i, j, k;
    double tmpdbl, a1, a2, a3, a4, a5;
    std::vector<double> row;  // A temperary vector;
    std::vector< std::vector< double > > mat;     // A temperary matrix;
    ifstream inFile( filename );
    std::getline(inFile, tmp) ;
    inFile >> nx >> ny >> nz >> nxm >> nym >> nzm >> nat >> ntyp;
    std::cout << "nx=" << nx  << endl\
              << "ny=" << ny  << endl\
              << "nz=" << nz  << endl\
              << "nat="<< nat << endl;
    inFile >> ibrav >> cell[0];
    
    std::getline(inFile, tmp);
    if ( ibrav == 1 ){
        cell[2] = cell[1] = cell[0];
        std::getline(inFile, tmp);
    } else if ( ibrav == 0 ) {
        std::getline(inFile, tmp);
        std::getline(inFile, tmp);
        inFile >> tmpdbl;
        inFile >> tmpdbl;
        inFile >> tmpdbl;
        cout << tmpdbl<<endl;
        cell[1] = cell[0];
        cell[2] = cell[0]*tmpdbl;
        std::getline(inFile, tmp); 
        std::getline(inFile, tmp);
    }
    std::cout << "ibrav=" << ibrav   << endl \
              << "a = "   << cell[0] << " bohr" << endl \
              << "b = "   << cell[1] << " bohr" << endl \
              << "c = "   << cell[2] << " bohr" << endl;
    for ( i = 0; i < ntyp; i++ ) {
        std::getline(inFile, tmp);
    }
    for ( i = 0; i < nat; i++ ) {
        std::getline(inFile, tmp);
    }
    while ( inFile >> a1 >> a2 >> a3 >> a4 >> a5 ){
        a.push_back(a1);
        a.push_back(a2);
        a.push_back(a3);
        a.push_back(a4);
        a.push_back(a5);
    }
    inFile.close();

    for ( i = 0; i < nx; i++ ) {
        mat.erase(mat.begin(), mat.end());
        for ( j = 0; j < ny; j++ ) {
            row.erase(row.begin(), row.end());
            for ( k = 0; k < nz; k++ ) {
                row.push_back(0.0);
            }
            mat.push_back(row);
        }
        cube.push_back(mat);
    }
                

    for ( k = 0; k < nz; k++ ){
        for ( j = 0; j < ny; j++ ){
            for ( i = 0; i < nx; i++ ){
                cube[i][j][k] = a[i+1+(j+1-1)*nx+(k+1-1)*nx*ny-1];
            }
        }
    }
    cout << " first 10 elements: " << endl;
    for ( i = 0; i < 10; i++ ) {
        cout << i << "  " << cube[i][0][0] << endl;
    }

}


void MatxVec(std::vector<double> &Res, \
             std::vector<std::vector<double> > Mat, \
             std::vector<double> Vec, unsigned N)
{
    /*
     * Multiply matrix Mat (N by N) and vector Vec (N by 1)
     */
    double tmp;
    for ( unsigned int i = 0; i < N; i++ ) {
        tmp = 0.0;
        for ( unsigned int j = 0; j < N; j++ ) {
            tmp += Mat[i][j]*Vec[j];
        }
        Res[i] = tmp;
    }
}

