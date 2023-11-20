#include <iostream>
#include <fstream>
#include "aspbe.hpp"

void writeleg(const double V_cryst, const s_pars my_pars, const s_grid* my_grid, std::ostream &fdat) {

    double V_liq, cl_protein, cl_AS, sol, supsat, rGro, rNuc;
  
    fdat << "# parameters:        "  << endl;
    fdat << "# "                     << std::endl;
    fdat << "# confil,f           "  << my_pars.confil <<           "    ...configuration file name" << endl;
    fdat << "# runid,i            "  << my_pars.runid <<            "    ...run-id, used for output files" << endl;
    fdat << "# " << std::endl;
    fdat << "# xtal size        = " << 10000*my_grid[0].avg << " - " 
                                    << 10000*my_grid[my_pars.NClasses_max-1].avg 
                                    << " mum"  << std::endl;
    fdat << "# bin widths       = " << 10000*my_grid[0].width << " - " 
                                    << 10000*my_grid[my_pars.NClasses_max-1].width 
                                    << " mum"  << std::endl;
    fdat << "# " << std::endl;
    fdat << "# name               "             << my_pars.name <<         "    ...protein name" << endl;
    fdat << "# asname             "           << my_pars.asname <<         "    ...anti-solvent name" << endl;
    fdat << "# SizeNucl           "         << my_pars.SizeNucl <<         "    ...average size of nucleus [cm]" << endl;
    fdat << "# " << std::endl;
    fdat << "# ga                 "               << my_pars.k_growth_a <<       "    ...growth kg, G = kg ((c-csat)/csat)^dg [cm/sec]" << endl;
    fdat << "# gb                 "               << my_pars.k_growth_b <<       "    ...growth dg, G = kg ((c-csat)/csat)^dg []" << endl;
    fdat << "# na                 "               << my_pars.k_nucl_a <<         "    ...nucleation kn, G = kn ((c-csat)/csat)^dn [#/cm3/sec]" << endl;
    fdat << "# nb                 "               << my_pars.k_nucl_b <<         "    ...nucleation dn, G = kn ((c-csat)/csat)^dn []" << endl;
    fdat << "# ak                 "               << my_pars.ak <<               "    ...aggregation kernel const" << endl;
    fdat << "# at                 "                << my_pars.at <<              "    ...aggregation type - 0: none, 1: const, 2: sum, 3: product, 4: (Li+Lj)^3, 5: G * (Li+Lj)^3" << endl;
    fdat << "# " << std::endl;
    fdat << "# ds_protein         "       << my_pars.ds_protein <<       "    ...density protein in crystal [g/cm3]" << endl;
    fdat << "# ds_water           "         << my_pars.ds_water  <<      "    ...density water in crystal [g/cm3]" << endl;
    fdat << "# printlevel         "       << my_pars.printlevel  <<      "    ...boolean, 0: print nothing, 1: print CSD of tn, 2: print CSD for all t" << endl;
    fdat << "# dt_0               "             << my_pars.dt_0 <<       "    ...initial time step [sec]" << endl;
    fdat << "# dt_write           "         << my_pars.dt_write <<       "    ...write interval for running output [sec]" << endl;
    fdat << "# t_end              "            << my_pars.t_end <<       "    ...total simulation time [sec]" << endl;
    fdat << "# sizeexp            "          << my_pars.sizeexp <<       "    ...exponent for grid size generation" << endl;
    fdat << "# typecsd0           "         << my_pars.typecsd0 <<       "    ...type of seed CSD (0: none, 1: square, 2: normal dist)" << endl;
    fdat << "# mcsd0              "            << my_pars.mcsd0 <<       "    ...average size of seed CSD" << endl;
    fdat << "# scsd0              "            << my_pars.scsd0 <<       "    ...stdev/halfwidth of seed CSD" << endl;
    fdat << "# ms_protein_0       "     << my_pars.ms_protein_0 <<       "    ...mass of protein in seed" << endl;
    fdat << "# shapefac           "         << my_pars.shapefac <<       "    ...shapefactor ratio (crystal volume)/cube" << endl;
    fdat << "# fr_AS              "            << my_pars.fr_AS <<       "    ...flow rate of antisolvent [cm3/sec]" << endl;
    fdat << "# Temp               "             << my_pars.Temp <<       "    ...Temperature [degC]" << endl;
    fdat << "# pH                 "    << my_pars.pH <<                  "    ...pH []" << endl;
    fdat << "# ml_protein_0       "     << my_pars.ml_protein_0 <<       "    ...mass of protein in solution at t=0 [g]" << endl;
    fdat << "# ml_AS_max          "        << my_pars.ml_AS_max <<       "    ...total volume of antisolvent to add" << endl;
    fdat << "# ml_water_0         "       << my_pars.ml_water_0 <<       "    ...mass of water in solution at t=0 [g]" << endl;
    fdat << "# nclasses           "         << my_pars.NClasses_max <<   "    ...number of grid points" << endl;
    fdat << "# dl_water_res       "     << my_pars.dl_water_res <<       "    ...density of water in AS solution [g/cm3]" << endl;
    fdat << "# dl_AS_res          "        << my_pars.dl_AS_res <<       "    ...density of water in AS solution [g/cm3]" << endl;
    fdat << "# const_supsat       "     << my_pars.const_supsat <<       "    ...constant super saturation (for testing)" << endl;
    fdat << "# " << std::endl;
    fdat << "# " << std::endl;
    
    fdat << setw(16) << "#      t      "; 
    fdat << setw(16) << "ml_protein";     
    fdat << setw(16) << "ml_AS";          
    fdat << setw(16) << "ml_water";       
    fdat << setw(16) << "cl_protein";     
    fdat << setw(16) << "V_c";            
    fdat << setw(16) << "V_liq";          
    fdat << setw(16) << "solubility";     
    fdat << setw(16) << "supsat";         
    fdat << setw(16) << "rGro";           
    fdat << setw(16) << "rNuc";           
    fdat << setw(16) << "cfr";            
    fdat << endl;

    V_liq          = volume_in_ml_from_masses_in_g(my_pars.name, my_pars.asname, my_pars.ml_protein_0, 0.0, my_pars.ml_water_0);
    cl_protein     = 1000.0*my_pars.ml_protein_0/V_liq;
    

    
    cl_AS          = 0.0;


    if(my_pars.name=="lysozyme"  && my_pars.asname=="nacl")    { sol = 346.01309430864*exp(-0.102667383846361*cl_AS);}
    else if(my_pars.name=="ifng" && my_pars.asname=="nh42so4") { sol = 1000000;}
    else                                                       { fprintf(stderr, "solubility not parameterized\n"); exit(1);}
    
    if(sol < cl_protein) supsat = (cl_protein - sol)/sol;
    else                 supsat =  0.0;
    rGro = my_pars.k_growth_a * pow(supsat,my_pars.k_growth_b);                 
    rNuc = my_pars.k_nucl_a   * pow(supsat,my_pars.k_nucl_b);                   
    
    fdat << setw(16) << scientific << setprecision(6) << 0;
    fdat << setw(16) << scientific << setprecision(6) << my_pars.ml_protein_0;  
    fdat << setw(16) << scientific << setprecision(6) << 0.0;                   
    fdat << setw(16) << scientific << setprecision(6) << my_pars.ml_water_0;    
    fdat << setw(16) << scientific << setprecision(6) << cl_protein;            
    fdat << setw(16) << scientific << setprecision(6) << V_cryst;               
    fdat << setw(16) << scientific << setprecision(6) << V_liq;                 
    fdat << setw(16) << scientific << setprecision(6) << sol;                   
    fdat << setw(16) << scientific << setprecision(6) << supsat;                
    fdat << setw(16) << scientific << setprecision(6) << rGro;                  
    fdat << setw(16) << scientific << setprecision(6) << rNuc;                  
    fdat << setw(16) << scientific << setprecision(6) << my_pars.fr_AS;         
    fdat << endl;
   
    fdat << "# ";
    
}
