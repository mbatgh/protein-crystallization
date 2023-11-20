
#include <iostream>
#include <fstream>
#include "aspbe.hpp"

double write_csd( const state_type cld, const state_type n, const state_type n0, const s_pars pars, const s_grid* grid ) {

  double totv, totv0, ntotv, ntotv0, cldtot, ncldtot, vcldtot, nvcldtot;
  int i;
  string csdfil = pars.runid + "-csd1"; 
  ofstream fcsd; 
  
  fcsd.open(csdfil.c_str());

  totv  = 0.0;
  totv0 = 0.0;
  ntotv = 0.0; 
  ntotv0  = 0.0; 
  cldtot  = 0.0; 
  ncldtot  = 0.0; 
  vcldtot  = 0.0; 
  nvcldtot = 0.0; 
  

  for(i=0; i < pars.NClasses_max; i++ ) {
    totv  += n[i]  * grid[i].width * pow(grid[i].avg, 3.0);
    totv0 += n0[i] * grid[i].width * pow(grid[i].avg, 3.0);
    ntotv  += n[i]  * pow(grid[i].avg, 3.0);
    ntotv0 += n0[i] * pow(grid[i].avg, 3.0);
    cldtot += cld[i];
    ncldtot += cld[i] / grid[i].width;
    vcldtot += cld[i] * pow(grid[i].avg,3.0);
    nvcldtot += cld[i] * pow(grid[i].avg,2.0);
  }

  if(totv==0.0) totv=1.0;
  if(totv0==0.0) totv0=1.0;
  if(ntotv==0.0) totv=1.0;
  if(ntotv0==0.0) totv0=1.0;
  if(cldtot==0.0) cldtot=1.0; 
  if(ncldtot==0.0) ncldtot=1.0; 
  if(vcldtot==0.0)  vcldtot=1.0; 
  if(nvcldtot==0.0)  nvcldtot=1.0; 

  fcsd  << "# CSDs at t=0 and t=" << pars.t_end << endl;
  fcsd  << setw(14) << "      #Li[mum]";  
  fcsd  << setw(14) << " Delta-Li[mum]";  
  fcsd  << setw(14) << "   Li^3[mum^3]";  
  fcsd  << setw(14) << "         rho-0";  
  fcsd  << setw(14) << "           rho";  
  fcsd  << setw(14) << "            n0";  
  fcsd  << setw(14) << "             n";  
  fcsd  << setw(14) << "    %Vol-0-abs";  
  fcsd  << setw(14) << "      %Vol-abs";  
  fcsd  << setw(14) << "   %Vol-0-norm";  
  fcsd  << setw(14) << "     %Vol-norm";  
  fcsd  << setw(14) << "    %n-cld-abs";  
  fcsd  << setw(14) << "   %n-cld-norm";  
  fcsd  << setw(14) << "  %vol-cld-abs";  
  fcsd  << setw(14) << " %vol-cld-norm" << endl; 

  for(i=0; i < pars.NClasses_max; i++ ) {
    fcsd  << setw(14)               << setprecision(4) << 10000.0*grid[i].avg;
    fcsd  << setw(14)               << setprecision(4) << 10000.0*grid[i].width;
    fcsd  << setw(14)               << setprecision(4) << pow(10000.0*grid[i].avg,3.0);
    fcsd  << setw(14) << scientific << setprecision(4) << n0[i];
    fcsd  << setw(14) << scientific << setprecision(4) << n[i];
    fcsd  << setw(14) << scientific << setprecision(4) << n0[i] * grid[i].width;
    fcsd  << setw(14) << scientific << setprecision(4) << n[i]  * grid[i].width;
    fcsd  << setw(14) << scientific << setprecision(4) << n0[i] * grid[i].width * pow(grid[i].avg,3.0) / totv0*100.0;
    fcsd  << setw(14) << scientific << setprecision(4) << n[i]  * grid[i].width * pow(grid[i].avg,3.0) / totv*100.0;    
    fcsd  << setw(14) << scientific << setprecision(4) << n0[i] * pow(grid[i].avg,3.0) / ntotv0*100.0;
    fcsd  << setw(14) << scientific << setprecision(4) << n[i]  * pow(grid[i].avg,3.0) / ntotv*100.0;
    fcsd  << setw(14) << scientific << setprecision(4) << cld[i] / cldtot*100.0;
    fcsd  << setw(14) << scientific << setprecision(4) << cld[i] / grid[i].width / ncldtot*100.0;
    fcsd  << setw(14) << scientific << setprecision(4) << cld[i] * pow(grid[i].avg,3.0) / vcldtot * 100.0;
    fcsd  << setw(14) << scientific << setprecision(4) << cld[i] * pow(grid[i].avg,2.0) / nvcldtot * 100.0 << endl;
  }

  fcsd.close();

  
  return(totv/pars.shapefac);
}
