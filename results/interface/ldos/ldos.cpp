#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <fstream>
#include <numeric>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include "EnergyFile.h"
#define PI (3.141592653589)
using namespace std;

// Pre-declaration
// Gaussian function;
double gaussian(double x, double mu, double sig);
// Integral of gaussian, i.e. error function.
double int_gaussian ( double x, double mu, double sig );

int main( int argc, char* argv[] )
{
    // Inputs:
    // 1) Energy range of the ldos, in eV, { start, end, step }
    double input[] = { 5.5, 11.0, 0.01};
    // 2) Name of your input file, by default wfc.files
    //    You may want to run the program with a different input file by doing
    //    ./ldos.x filename
    //
    //    Each line of inFile should contains 4 items:
    //    a) [Float] Energy of this wavefunction, in eV
    //    b) [Float] Gaussian broadening you want use, in eV, e.g. 0.1
    //    c) [Float] Norm of the wfc, not used, supposed to be 1.0
    //    d) [String] Where the wave function is. The file is supposed to end with .plavz. No quotes.
    string inFile = "wfc.files";
    if ( argc == 2 )
        inFile = argv[1];
    // 3) The delta parameter used in Anh's
    //    APPLIED PHYSICS LETTERS 102, 241603 (2013)
    //    to determine HOMO LUMO edges.
    //    This parameter is tricky.
    //    Try to tune this parameter and see how the plot looks like.
    const double delta=0.003; // eV

    // Outputs:
    //    You cannot modify those filenames in command-line.
    //    If you want, modifty them below, remember to recompile the code.
    //
    // ldosFile: LDOS, by default cmap.dat
    //    First line: number of date points along z, number of energy steps
    //    From 2nd line, each line contains three numbers
    //    a) Position along z direction, in Bohr
    //    b) Energy, in eV
    //    c) LDOS at this position of this energy
    // HLFile: HOMO and LUMO orbitals
    //    homo
    //    Position[bohr] Energy[eV]
    //    lumo
    //    Position[bohr] Energy[eV]
    // dosFile: density of states, not scaled to one
    //    Energy[eV]  dos
    string ldosFile, HLFile, dosFile;
    ldosFile = "cmap.dat";
    HLFile  = "homo_lumo.dat";
    dosFile = "dos.dat";

    //
    double FL, tmp_dbl;
    std::vector< std::vector<double> > xmap, ymap;
    std::vector<double> tmp_v1, tmp_v2;
    EnergyFile ef( inFile, input );

    // Compute density of states and write it to file 
    std::vector<double> dos( ef.e_axis.size() ); 
    ofstream dos_out;
    dos_out.open( dosFile, ios::out );
    for ( unsigned int iex = 0; iex < ef.e_axis.size(); iex++ ) { 
        dos[iex] = 0.0;
        for ( unsigned int ie = 0; ie <ef.ne; ie++ ) {
            dos[iex] += gaussian( ef.e_axis[iex], ef.e[ie], ef.sigs[ie])*ef.norms[ie];
        }
        dos_out << ef.e_axis[iex] << "  " << dos[iex] << endl;
        if ( iex == 0 ) { tmp_dbl = dos[iex]; FL = ef.e_axis[iex]; }
        else if ( tmp_dbl > dos[iex] ) { tmp_dbl = dos[iex]; FL = ef.e_axis[iex]; }
    }
    dos_out.close();
 
    cout << "Fermi Level: " << FL << endl;
    // 
    // Initialize xmap and cmap
    // 
    for ( unsigned int i = 0; i < ef.x_axis.size(); i++ ) {
        for ( unsigned int j = 0; j < ef.e_axis.size(); j++ ) {
            tmp_v1.push_back(ef.x_axis[i]);
            tmp_v2.push_back(ef.e_axis[j]);
        }
        xmap.push_back(tmp_v1);
        ymap.push_back(tmp_v2);
    }
    
    // 
    // Generate fmap, which is interpolated wfc along x_axis;
    // By generating fmap and keeping it in memory, we reduce the cost of generating cmap
    // because we don't have to do interpolation in every loop.
    // 
    // Cubic spline interpolation was done with GNU scientific library ( GSL ).
    // 
    std::vector< std::vector<double> > cmap (ef.x_axis.size(), std::vector<double> (ef.e_axis.size()));
    std::vector< std::vector<double> > fmap (ef.x_axis.size(), std::vector<double> (ef.ne));
    double *xpt = new double[ef.nx];
    double *ypt = new double[ef.nx];
    for ( unsigned int ie = 0; ie < ef.ne; ie++ ) {
        gsl_interp_accel *acc = gsl_interp_accel_alloc();
        gsl_spline *spline_cubic = gsl_spline_alloc(gsl_interp_cspline, ef.nx);
        xpt = &ef.x[0];
        ypt = &ef.wfcs[ie][0];
        gsl_spline_init(spline_cubic, xpt, ypt, ef.nx);
        for (unsigned int ixx = 0; ixx < ef.x_axis.size(); ixx++ ) {
            fmap[ixx][ie] = gsl_spline_eval(spline_cubic, ef.x_axis[ixx], acc);
        }
        gsl_spline_free(spline_cubic);
        gsl_interp_accel_free(acc);
    }
  
    
    // Generate cmap and write cmap to file to be processed by python matplotlib
    ofstream cmapFile;
    cmapFile.open( ldosFile, ios::out );
    cmapFile << ef.x_axis.size() << " " << ef.e_axis.size() << " " << endl;
    for ( unsigned int ixx = 0; ixx < ef.x_axis.size(); ixx++ ) {
        for ( unsigned int iex = 0; iex < ef.e_axis.size(); iex++ ) {
            cmap[ixx][iex] = 0.0;
            for (unsigned int ie = 0; ie < ef.ne; ie++ ) {
                cmap[ixx][iex] += fmap[ixx][ie] * \
                                  gaussian(ef.e_axis[iex], ef.e[ie], ef.sigs[ie]) * ef.norms[ie];
            }
            cmapFile << ef.x_axis[ixx] << " " << ef.e_axis[iex] << " " << cmap[ixx][iex] << endl;
        }
    }
    cmapFile.close();
    
    // 
    // Find out HOMO and LUMO of the system
    // 
    unsigned tmp_iex;
    double tmp_homoe, tmp_lumoe;
    std::vector<double> homoe( ef.x_axis.size() ), lumoe( ef.x_axis.size() );
    std::vector<double> Integral( ef.x_axis.size() ), tmp_int( ef.e_axis.size() );
    for ( unsigned int ixx = 0; ixx < ef.x_axis.size(); ixx++ ) {
        Integral[ixx] = 0.0;
        for ( unsigned int ie = 0; ie < ef.ne; ie++ ) {
            if ( ef.e[ie] <= FL - 6*ef.sigs[ie] ) {
                Integral[ixx] += fmap[ixx][ie];
            } else if ( ef.e[ie] > FL - 6*ef.sigs[ie] && ef.e[ie] <= FL + 6*ef.sigs[ie] ) {
                Integral[ixx] += fmap[ixx][ie]*ef.norms[ie]*int_gaussian(FL, ef.e[ie], ef.sigs[ie]);
            } else {
                break;
            }
        }
    }
    for ( unsigned int ixx = 0; ixx < ef.x_axis.size(); ixx++ ) { 
        tmp_iex = 0;
        for ( unsigned int iex = 0; iex < ef.e_axis.size(); iex++ ) {
            // Compute CBM
            tmp_int[iex]=0.0;
            if ( ef.e_axis[iex] <= FL ) { continue; }
            else {
                for ( unsigned int ie = 0; ie < ef.ne; ie++ ) {
                    if ( FL-6*ef.sigs[ie] <= ef.e[ie] && ef.e[ie] <= ef.e_axis[iex] + 6*ef.sigs[ie] ) {
                        tmp_int[iex] += fmap[ixx][ie]*ef.norms[ie]*\
                                        ( int_gaussian(ef.e_axis[iex], ef.e[ie], ef.sigs[ie])-\
                                          int_gaussian(FL, ef.e[ie], ef.sigs[ie]) );
                    }
                }
            tmp_int[iex] = tmp_int[iex] - delta*Integral[ixx];
            if ( tmp_iex == 0 ) { tmp_iex = iex; tmp_lumoe = tmp_int[iex]; }
            else if ( fabs(tmp_int[iex]) < fabs(tmp_lumoe) ) { tmp_iex = iex; tmp_lumoe = tmp_int[iex]; }
            }
        }
        lumoe[ixx] = ef.e_axis[tmp_iex];
    }
    for ( unsigned int ixx = 0; ixx < ef.x_axis.size(); ixx++ ) {    
        tmp_iex = 0;
        for ( unsigned int iex = 0; iex < ef.e_axis.size(); iex++ ) {
            // Compute VBM
            tmp_int[iex]=0.0;
            if ( ef.e_axis[iex] >= FL ) { break; }
            else {
                for ( unsigned int ie = 0; ie < ef.ne; ie++ ) {
                    if ( ef.e_axis[iex]-6*ef.sigs[ie] <= ef.e[ie] && ef.e[ie] <= FL + 6*ef.sigs[ie] ) {
                        tmp_int[iex] += fmap[ixx][ie]*ef.norms[ie]*\
                                        ( int_gaussian(FL, ef.e[ie], ef.sigs[ie])-\
                                          int_gaussian(ef.e_axis[iex], ef.e[ie], ef.sigs[ie]) );
                    }
                }
            tmp_int[iex] = tmp_int[iex] - delta*Integral[ixx];
            if ( iex == 0 ) { tmp_iex = iex; tmp_homoe = tmp_int[iex]; }
            else if ( fabs(tmp_int[iex]) < fabs(tmp_homoe) ) { tmp_iex = iex; tmp_homoe = tmp_int[iex]; }
            }
        }
        homoe[ixx] = ef.e_axis[tmp_iex];
    }

    // Output HOMO and LUMO
    ofstream MO;
    MO.open( HLFile, ios::out );
    MO << "homo" << endl;
    for ( unsigned int ixx = 0; ixx < ef.x_axis.size(); ixx++ ) {
        MO << ef.x_axis[ixx] << "  " << homoe[ixx] << endl;
    }
    MO << "lumo" << endl;
    for ( unsigned int ixx = 0; ixx < ef.x_axis.size(); ixx++ ) {
        MO << ef.x_axis[ixx] << "  " << lumoe[ixx] << endl;
    }
    MO.close();

    return 0;
}

double gaussian(double x, double mu, double sig)
{
    return 1.0/(sqrt(2*PI)*sig)*exp(-pow((x-mu)/sig,2.0)/2.0);
}

double int_gaussian ( double x, double mu, double sig )
{
    // Integral of gaussian(t, mu, sig) from -infinity to x;
    //   ~ x
    //  |            Gau(t,mu,sig) dt
    // ~  -infinity
    return 0.5*erfc((mu-x)/(sqrt(2)*sig));
}
