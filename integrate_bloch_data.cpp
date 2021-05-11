#include <iostream>
#include <fstream>
#include <boost/array.hpp>
#include <cmath>
#include <chrono>
#include <vector>
#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;
using namespace std::chrono;

const double gamma_n = -1.83247171e8;
const double gamma_He = -2.03789e8;

typedef vector< double > state_type;

double u_mag(double x1, double x2, double x3)
{
    return sqrt(pow(x1,2)+pow(x2,2)+pow(x3,2));
}

void neutron( const state_type &x , state_type &dxdt , double t )
{
    dxdt[0] = gamma_n * x[4] * x[1];
    dxdt[1] = gamma_n * ( x[3] * x[2] - x[4] * x[0]);
    dxdt[2] = -gamma_n * x[3] * x[1];
}

void helium( const state_type &y , state_type &dydt , double t )
{
    dydt[0] = gamma_He * y[4] * y[1];
    dydt[1] = gamma_He * ( y[3] * y[2] - y[4] * y[0]);
    dydt[2] = -gamma_He * y[3] * y[1];
}

int main(int argc, char **argv)
{
    vector<double> t;
    double u;

    ifstream infile;
    infile.open("Input_Data/Times.txt");
    
    double tmp;
    while(!infile.eof()){
        infile >> tmp;
        t.push_back(tmp);
    }
    t.pop_back();
    infile.close();

    const double dt = t[1]-t[0];
    const int num_steps = t.size();
    const int num_particles = 10;

    double** B0;
    double** B1;
    B0 = new double* [num_steps];
    B1 = new double* [num_steps];
    for(int i=0; i<num_steps; i++){
        B0[i] = new double [num_particles];
        B1[i] = new double [num_particles];
    }

    ifstream infile1;
    infile1.open("Input_Data/B0_data.txt");
    ifstream infile2;
    infile2.open("Input_Data/B1_data.txt");
    for(int i=0; i<num_steps; i++){
        for(int j=0; j<num_particles; j++){
            infile1 >> B0[i][j];
            infile2 >> B1[i][j];
        }
    }
 
    infile1.close();
    infile2.close();

    ofstream outfile;

    runge_kutta_dopri5< state_type > stepper;

    auto start = high_resolution_clock::now();

    for(int j=0; j<num_particles; j++){

        outfile.open("Neutron_Data/Neutron_"+to_string(j)+"_data.txt");
        state_type x = { 0 , 0 , 1.0 , 0 , 0}; // initial conditions

        for(size_t i=0; i<num_steps; i++){

            x[3] = B1[i][j];
            x[4] = B0[i][j];
            stepper.do_step( neutron, x , t[i] , dt);
            u = u_mag(x[0],x[1],x[2]);
            outfile << t[i] << "\t" << x[0] << "\t" << x[1] << "\t" << x[2] << "\t" << x[4] << "\t" << u << "\n";
        }
        outfile.close();
    }

    for(int j=0; j<num_particles; j++){

        outfile.open("Helium_Data/Helium_"+to_string(j)+"_data.txt");
        state_type y = { 0 , 0 , 1.0 , 0 , 0}; // initial conditions

        for(size_t i=0; i<num_steps; i++){

            y[3] = B1[i][j];
            y[4] = B0[i][j];
            stepper.do_step( helium, y , t[i] , dt);
            u = u_mag(y[0],y[1],y[2]);
            outfile << t[i] << "\t" << y[0] << "\t" << y[1] << "\t" << y[2] << "\t" << y[3] << "\t" << u << "\n";
        }
        outfile.close();
    }


    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    cout << "Integration time = " << duration.count() << " milliseconds" << "\n";
    cout << "Number of steps = " << num_steps << "\n";
    cout << "Number of particles = " << num_particles << "\n";

    return 0;
}
