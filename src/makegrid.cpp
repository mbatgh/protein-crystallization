
#include "aspbe.hpp"
#include "time.h"

int makegrid( const s_pars pars, s_grid* grid)
{

    double x,lavg,tmp,maxkern=1.0;
    int i,j,k;

    grid[0].low    = 0.0;
    grid[0].avg    = pars.SizeNucl;
    grid[0].high   = 2.0*pars.SizeNucl;
    grid[0].width  = 2.0*pars.SizeNucl;
    grid[0].volume = pow(pars.SizeNucl,3.0);
   


    for(i=1 ; i < pars.NClasses_max; i++ )
    {
        x = (double)i;
        grid[i].low    = grid[i-1].high;
        grid[i].high   = grid[i].low + (2.0 * pars.SizeNucl) * pow(pars.sizeexp,x);
        grid[i].avg    = 0.5*(grid[i].high+grid[i].low);
        grid[i].width  = grid[i].high - grid[i].low;
        grid[i].volume = pow(grid[i].avg,3.0);

    }
    


    clock_t tStart = clock();

    if(pars.at > 0) {
      
      for(i=0 ; i < pars.NClasses_max-1; i++ )
      {
        naggr[i]=0;
      
        for(k=0 ; k < pars.NClasses_max; k++ )
        {
          for(j=0 ; j < pars.NClasses_max; j++ )
          {
            lavg = pow(grid[k].volume+grid[j].volume, 1.0/3.0);
            if(lavg <= grid[i+1].avg && lavg > grid[i].avg)
            {
              naggr[i]++;
              if(naggr[i]>MAX_N_AGGR) {
                cout << "need to recompile with MAX_N_AGGR > " << naggr[i] << endl;
                exit(1);
	      }
            }
            andx1[i][naggr[i]] = k;
            andx2[i][naggr[i]] = j;
            tmp = (1.0 - pow(lavg/grid[i+1].avg,3.0)) / (1.0 - pow(grid[i].avg/grid[i+1].avg,3.0));
            ahw1[i][naggr[i]] = 0.5*tmp;
            ahw2[i][naggr[i]] = 0.5*(1.0-tmp);
          }
        }
      }
      
      if(pars.at==4 || pars.at==5) {
	
        for(i=0 ; i < pars.NClasses_max; i++ )
        {
          for(j=0 ; j < pars.NClasses_max; j++ )
          {
            nakern[i][j] = pow(grid[i].avg+grid[j].avg, 3.0);
          }
        }
        
        maxkern = pow(2.0*grid[pars.NClasses_max-1].avg, 3.0);
        
      } else if(pars.at==3) {
	
        for(i=0 ; i < pars.NClasses_max; i++ )
        {
          for(j=0 ; j < pars.NClasses_max; j++ )
          {
            nakern[i][j] = grid[i].avg*grid[j].avg;
          }
        }
        
        maxkern = grid[pars.NClasses_max-1].avg*grid[pars.NClasses_max-1].avg;
	
      } else if(pars.at==2) {
	
        for(i=0 ; i < pars.NClasses_max; i++ )
        {
          for(j=0 ; j < pars.NClasses_max; j++ )
          {
            nakern[i][j] = grid[i].avg+grid[j].avg;
          }
        }
        
        maxkern = grid[pars.NClasses_max-1].avg+grid[pars.NClasses_max-1].avg;
        
      } else if(pars.at==1) {
	
	for(i=0 ; i < pars.NClasses_max; i++ )
        {
          for(j=0 ; j < pars.NClasses_max; j++ )
          {
            nakern[i][j] = pars.ak;
          }
        }
        
        maxkern = pars.ak;
      }

      tmp=0.0;

      for(i=0 ; i < pars.NClasses_max; i++ )
      {
        for(j=0 ; j < pars.NClasses_max; j++ )
        {
          nakern[i][j] = nakern[i][j] / maxkern * pars.ak;
        }
      }
    }
   
    return(((int)(clock() - tStart)/CLOCKS_PER_SEC));
}
