#include "aspbe.hpp"

void check_parameters(s_pars my_pars)
{
  if(my_pars.NClasses_max>=MAX_N_C) {
    cout << "MAX_N_C too small - recompile!" << endl;
    exit(1);
  }
  
  if(my_pars.at<0 || my_pars.at>5)  {
    cout << "wrong aggregation type" << endl;
    exit(1);
  }
}
