
#include <vector>
#include <iostream>
#include "aspbe.hpp"

double makeseed(state_type& vd_n, const s_pars pars, const s_grid* grid)
{

    double fac=1.0;
    double totv=0.0;
    double mpis=0.0;
    int i;

    vd_n[0] = 0;
    
    if(pars.typecsd0==2) {



      for(i=0 ; i < pars.NClasses_max; i++ )
      {
        vd_n[i] = 1.0/(sqrt(2.0*M_PI)*pars.scsd0) *  
	  exp( -(grid[i].avg - pars.mcsd0)*(grid[i].avg - pars.mcsd0) / (2.0 * pars.scsd0 * pars.scsd0)) * grid[i].width;

        if(grid[i].avg > pars.mcsd0 + 3.0*pars.scsd0 || grid[i].avg < pars.mcsd0 - 3.0*pars.scsd0 ) vd_n[i]=0.0;
        totv += vd_n[i] * pow(grid[i].avg,3.0);
      }
    } 
    else if(pars.typecsd0==1) {



      if(pars.mcsd0-pars.scsd0<0.0 || pars.mcsd0+pars.scsd0>grid[pars.NClasses_max-1].avg) {
        cout << "error 1 in makegrid" << endl;
        exit(1);
      }

      for(i=0 ; i < pars.NClasses_max; i++ )
      {
        if(grid[i].avg >= pars.mcsd0-pars.scsd0 && grid[i].avg <= pars.mcsd0+pars.scsd0) {
          vd_n[i] = 1.0 * grid[i].width;
        } else {
          vd_n[i] = 0.0;
        }
        totv += vd_n[i] * pow(grid[i].avg, 3.0);
      }

    } else {



        return(0.0);
    }

    if(pars.mcsd0+pars.scsd0 > grid[pars.NClasses_max-1].high) {
      cout << "error 1 in makegrid" << endl;
      cout << "largest size: " << grid[pars.NClasses_max-1].high << endl;
      exit(1);
    }

    if(pars.ms_protein_0==0.0) {
      cout << "error 2 in makegrid" << endl;
      cout << "ms_protein_0 = " << pars.ms_protein_0 << endl;
      exit(1);
    }

    mpis = totv*pars.shapefac*pars.ds_protein;
    fac = pars.ms_protein_0/mpis;



    
    for(i=0 ; i < pars.NClasses_max; i++ ) { vd_n[i] *= fac/grid[i].width; }

    totv=0.0;

    for(i=0 ; i < pars.NClasses_max; i++ ) {
      totv += vd_n[i] * grid[i].width * pow(grid[i].avg, 3.0);
    }




 
    return(totv);
}
