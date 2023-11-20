
#include <iostream>
#include <fstream>
#include <vector>
#include "aspbe.hpp"
#include "spline.hpp"

s_dprop ana_csd(const state_type n, const s_pars my_pars, const s_grid *my_grid) {

  using namespace magnet::math;

  int i, imod, i_bot_left, i_bot_right;
  s_dprop dprop;
  state_type nnv(my_pars.NClasses_max, 0.0 );
  double y,xmax,ymax,vtot,ntot,xleft,xright;
  Spline spline;
  
  dprop.mode = 0.0;
  dprop.fwhm = 0.0;
  dprop.pvat = 0.0;
  
  xmax = 0.0;
  ymax = 0.0;
  imod = 0;

  ntot=0.0;
  for(i=0;i<my_pars.NClasses_max;i++) {
    nnv[i] = n[i]*my_grid[i].width*pow(my_grid[i].avg,3.0);
    ntot += nnv[i];
  }
  
  for(i=0;i<my_pars.NClasses_max;i++) {
    nnv[i] /= ntot/100.0;
  }

  vtot=0.0;
  for(i=0;i<my_pars.NClasses_max;i++) { 
    if(my_grid[i].avg>my_pars.sizethresh)  dprop.pvat+= nnv[i];
    vtot += nnv[i];
  }
  dprop.pvat /= vtot/100.0;

  ymax=0.0;
  for(i=0; i < my_pars.NClasses_max; i++ ) {
    if(nnv[i] > ymax) { 
      ymax = nnv[i]; 
      imod = i;
    }
  }

  if(imod<3) {
    dprop.mode = 0.0;
    dprop.fwhm = 0.0;
    dprop.pvat = 0.0;
    return(dprop);
  }
    
  if(imod>my_pars.NClasses_max-3) {
    dprop.mode = my_grid[my_pars.NClasses_max-1].avg;
    dprop.fwhm = 0.0;
    dprop.pvat = 100.0;
    return(dprop);
  }
  
  xmax = my_grid[imod].avg;

  
  
  dprop.mode = xmax;
  

  
  i=imod;
  do {
    i--;
    y = nnv[i];
  } while(y>ymax/2.0&&i>0);
  i_bot_left = i;
  
  i=imod;
  do {
    i++;
    y = nnv[i];
  } while(y>ymax/2.0&&i<my_pars.NClasses_max-1);
  i_bot_right = i;

  xleft  = my_grid[i_bot_left].avg + 
           (ymax/2.0-nnv[i_bot_left]) / (nnv[i_bot_left]-nnv[i_bot_left+1]) * 
           (my_grid[i_bot_left+1].avg-my_grid[i_bot_left].avg);
  xright = my_grid[i_bot_right].avg - 
           (ymax/2.0-nnv[i_bot_right]) / (nnv[i_bot_right]-nnv[i_bot_right-1]) * 
           (my_grid[i_bot_right].avg-my_grid[i_bot_right-1].avg);  

  dprop.fwhm = xright-xleft;

  return(dprop); 
}
