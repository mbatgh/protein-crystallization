
#include "aspbe.hpp"

void set_defaults(s_pars& pars) {



    pars.t_end          = 10000.0;                 
    pars.dt_0           = 0.1;                     
    pars.fr_AS          = 1.0;                     
    pars.dt_write       = 1000.0;                  
    pars.runid          = "test";                  
    pars.rseed          = 12345;                   
    pars.sizeexp        = 1.003;                   
    pars.NClasses_max   = 1000;                    
    pars.printlevel     = 1;                       
    pars.const_supsat   = 0;                       
    pars.fr_AS2         = 0.0;                     
    pars.fr_AS3         = 0.0;                     
    pars.t1             = 1.0e+32;                 
    pars.t2             = 1.0e+32;                 
    pars.DOCLD          = 0;                       
    pars.sizethresh     = 0.003;                   
    pars.nrot           = 10000;                   
    pars.pcp0           = 5.0;                     
        


    pars.name           = "lysozyme";              
    pars.asname         = "nacl";                  
    pars.k_growth_a     = 1.0e-07;                 
    pars.k_growth_b     = 2.0;                     
    pars.k_nucl_a       = 1.0e+10;                 
    pars.k_nucl_b       = 6.0;                     
    pars.shapefac       = 0.5;                     
    pars.SizeNucl       = 1.0e-05;                 
    pars.ds_water       = 0.75;                    
    pars.ds_protein     = 1.48;                    



    pars.typecsd0       = 0;                       
    pars.mcsd0          = 0.002;                   
    pars.scsd0          = 0.001;                   
    pars.ms_protein_0   = 0.0;                     



    pars.ml_AS_max      = 100.0;                   
    pars.ml_protein_0   = 10.0;                    
    pars.ml_water_0     = 100.0;                   
    pars.dl_water_res   = 1.0;                     
    pars.dl_AS_res      = 0.1;                     
    pars.ak             = 0.0;                     
}
