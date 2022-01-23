#include <iostream>
#include <string>
#include <vector>
#include <fstream>
using namespace std;

class EnergyFile
{
    public:
        std::string fname;
        unsigned int ne;                               // number of energy levels readin from file
        unsigned int nx;                               // number of x readin from file
        int times;
        std::vector<double>  x_axis, e_axis;  // energy sample points in the plot, 'case we need more points in plot
        std::vector<double>  x;
        double               xmax, xmin;
        std::vector<double>  e;               // energy read from file
        double               emax, emin;
        std::vector<double>  sigs;
        std::vector<double>  norms;
        std::vector<string>  names;
        std::vector< vector<double> >  wfcs;
        void   readin ( string fname );
        EnergyFile( string fname, double *input );
};

void EnergyFile::readin ( string fname )
{
    double tmp_e, tmp_sig, tmp_norm, tmp_x, tmp_dbl;
    string tmp_name, tmp_str;
    vector<double>  tmp_vec;

    ifstream enfile;
    enfile.open( fname, ios::in );
    this->ne = 0;
    while ( enfile >> tmp_e >> tmp_sig >> tmp_norm >> tmp_name ) {
        this->e.push_back(tmp_e);
        this->sigs.push_back(tmp_sig);
        this->norms.push_back(tmp_norm);
        this->names.push_back(tmp_name);
        this->ne += 1;
    }
    this->emax = e.back(); 
    this->emin = e[0];
    enfile.close();
    
    ifstream wfcFile;
    this->nx = 0;
    for ( int ie = 0; ie < this->ne; ie ++ ) {
        wfcFile.open( this->names[ie], ios::in );
        std::getline( wfcFile, tmp_str );
        tmp_vec.clear( );
        while ( wfcFile >> tmp_x >> tmp_dbl ) {
            if ( ie == 0 ) {
                this->x.push_back(tmp_x);
                this->nx += 1;
            }
            tmp_vec.push_back( tmp_dbl );
        }
        this->wfcs.push_back( tmp_vec );
        wfcFile.close( );
    }
    this->xmax = x.back();
    this->xmin = x[0];
}

// Constructor
EnergyFile::EnergyFile( string fname, double *input )
{
    double tmp_x, tmp_e;
    double dx;
    this->fname = fname;
    this->readin(this->fname);

    cout << "   Input range of energy: [" << input[0]<< ", " << input[1]<< "]" << endl;
    cout << "Possible range of energy: [" << this->emin << ", " << this->emax << "]" << endl;

    if ( input[0] < this->emin ) {
        cout << "Energy requested is too low!" << endl;
        exit(0);
    }else if ( input[1] > this->emax ) {
        cout << "Energy requested is too high!" << endl;
        exit(0);
    }
    // Initialize energy axis 
    tmp_e = input[0];
    while ( tmp_e <= input[1] ) {
        this->e_axis.push_back( tmp_e );
        tmp_e += input[2];
    }
    // Initialize x_axis, x_axis has more values than x read from wfc files.
    this->times = 5;
    dx = ( this->xmax - this->xmin ) / ( this->x.size() * this->times );
    tmp_x = this->xmin;
    while ( tmp_x <= this->xmax ) {
        this->x_axis.push_back( tmp_x );
        tmp_x += dx;
    }
}
