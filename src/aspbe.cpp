




#include <iostream>
#include <fstream>
#include <numeric>
#include <cmath>
#include <vector>
#include <string>
#include <iomanip>
#include <time.h>
#include <boost/numeric/odeint.hpp>
#include "aspbe.hpp"
#include "myode.hpp"

int andx1[MAX_N_C][MAX_N_AGGR];
int andx2[MAX_N_C][MAX_N_AGGR];
int naggr[MAX_N_C];
double nakern[MAX_N_C][MAX_N_C];
double ahw1[MAX_N_C][MAX_N_AGGR];
double ahw2[MAX_N_C][MAX_N_AGGR];
double tcp0;

int main( const int argc , const char **argv )
{
    s_pars my_pars;
    s_idat my_idat;
    s_dprop csdprops;
    s_dprop cldprops;
    int i, n_steps, tat;
    double V_cryst;
    ofstream fdat,frcsd;
    clock_t tStart = clock();

    tcp0=0.0;


 
    set_defaults(my_pars);
    read_parameters(argc,argv,my_pars);



    state_type vd_n(    my_pars.NClasses_max+3 , 0.0 );
    state_type vd_n0(   my_pars.NClasses_max ,   0.0 );
    state_type vd_cld(  my_pars.NClasses_max ,   0.0 );

    vd_n[my_pars.NClasses_max]   = my_pars.ml_protein_0;
    vd_n[my_pars.NClasses_max+1] = 0.0;
    vd_n[my_pars.NClasses_max+2] = my_pars.ml_water_0;
    


    s_grid* my_grid;
    my_grid = new s_grid[my_pars.NClasses_max];    
    tat = makegrid(my_pars, my_grid);



    V_cryst = makeseed(vd_n, my_pars, my_grid);
    for(i=0; i < my_pars.NClasses_max; i++ ) { vd_n0[i] = vd_n[i]; }



    string outfil = my_pars.runid + "-data";
    if(my_pars.printlevel>0) fdat.open (outfil.c_str());

    if(my_pars.printlevel==2) {
      string rcsdfil = my_pars.runid + "-csdn";
      frcsd.open (rcsdfil.c_str());
    }
    
    if(my_pars.printlevel>0) writeleg(V_cryst, my_pars, my_grid, fdat);
    


    n_steps = integrate_const(    
              make_controlled( 1e-06 , 1e-06 , stepper_type() ), 
              population(my_pars,my_grid,my_idat),
              vd_n, 
              0.0 , my_pars.t_end, my_pars.dt_write, streaming_observer( fdat, frcsd, my_grid, my_idat ));



                             csdprops =   ana_csd(        vd_n, my_pars, my_grid);
    if(my_pars.DOCLD==1)     cldprops =   ana_cld(vd_cld, vd_n, my_pars, my_grid);
    if(my_pars.printlevel>0) V_cryst = write_csd(vd_cld, vd_n, vd_n0, my_pars, my_grid);

    if(my_pars.printlevel>0) {
      fdat << "# " << n_steps << " steps of numerical inegration ";
      fdat << "in " << ((double)(clock() - tStart)/CLOCKS_PER_SEC) << "(" << tat << ") seconds. " << endl;
      fdat << "final totel crystal volume = " << V_cryst  <<  " cm3" << endl;
    }
  
    cout << my_pars.runid << "  ";

    cout << setw(12) << scientific << setprecision(4) << my_pars.k_nucl_a;
    cout << setw(12) << scientific << setprecision(4) << my_pars.k_nucl_b;
    cout << setw(12) << scientific << setprecision(4) << my_pars.k_growth_a;
    cout << setw(12) << scientific << setprecision(4) << my_pars.k_growth_b;

    cout << setw(12) << scientific << setprecision(4) << tcp0;
    cout << setw(12) << scientific << setprecision(4) << 1.0e+04 * csdprops.mode;
    cout << setw(12) << scientific << setprecision(4) << 1.0e+04 * csdprops.fwhm;
    cout << setw(12) << scientific << setprecision(4) << csdprops.pvat;
    if(my_pars.DOCLD==1) {
      cout << setw(12) << scientific << setprecision(4) << 1.0e+04 * cldprops.mode;
      cout << setw(12) << scientific << setprecision(4) << 1.0e+04 * cldprops.fwhm;
      cout << setw(12) << scientific << setprecision(4) << cldprops.pvat;
    }
    cout << endl;

    if(my_pars.printlevel==2) frcsd.close();
    if(my_pars.printlevel>0) fdat.close();
/*    free(my_grid);*/

    return 0;
}
