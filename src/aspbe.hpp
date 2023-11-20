
#ifndef ASPBE_H
#define ASPBE_H

#include <vector>
#include <string>
#include <boost/numeric/odeint.hpp>
#include <eigen3/Eigen/Dense>

using namespace Eigen;
using namespace boost::numeric::odeint;
using namespace std;

const int MAX_N_C = 1001;
const int MAX_N_AGGR = 3000;

extern int    andx1[MAX_N_C][MAX_N_AGGR];
extern int    andx2[MAX_N_C][MAX_N_AGGR];
extern double ahw1[MAX_N_C][MAX_N_AGGR];
extern double ahw2[MAX_N_C][MAX_N_AGGR];
extern int    naggr[MAX_N_C];
extern double nakern[MAX_N_C][MAX_N_C];
extern double tcp0;



struct s_pars
{

  string confil;
  string runid;
  string asname;
  double t_end;
  double dt_0;
  double fr_AS;
  double fr_AS2;
  double fr_AS3;
  double t1;
  double t2;
  double Temp;
  double pH;
  double dt_write;
  double sizeexp;
  double const_supsat;
  double sizethresh;
  double pcp0;
  int    NClasses_max;
  int    printlevel;
  int    DOCLD;

  string name;
  double k_growth_a;
  double k_growth_b;
  double k_nucl_a;
  double k_nucl_b;
  double shapefac;
  double SizeNucl;

  int    typecsd0;
  double mcsd0;
  double scsd0;
  double ms_protein_0;

  double ml_AS_max;
  double ml_protein_0;
  double ml_water_0;
  double ds_protein;
  double ds_water;
  double dl_AS_res;
  double dl_water_res;
  double ak;
  int    at;
  int    rseed;
  int    nrot;
};

struct s_grid
{
  double low;
  double avg;
  double high;
  double width;
  double volume;
};

struct s_idat {
  int    csdoft;
  double cl_protein;
  double V_liq;
  double V_c;
  double solubility;
  double supsat;
  double rGro;
  double rNuc;
  double cfr;
};

struct s_dprop {
  double pvat;
  double mode;
  double fwhm;
};

typedef struct s_idat s_idat;
typedef struct s_pars s_pars;
typedef struct s_grid s_grid;
typedef vector< double > state_type;
typedef struct s_dprop s_dprop;

void    set_defaults(s_pars&);
void    read_parameters( const int, const char**, s_pars&);
int     makegrid(const s_pars, s_grid*);
double  makeseed(state_type&, const s_pars, const s_grid*);
double  write_csd(const state_type, const state_type, const state_type, const s_pars, const s_grid*);
s_dprop ana_csd(const state_type, const s_pars, const s_grid*);
s_dprop ana_cld(state_type&, const state_type, const s_pars, const s_grid*);
void    writeleg(const double , const s_pars, const s_grid*, std::ostream&);
void    randrot(Vector3d*, int);
double  randd(void);
double  getpp30csd(const state_type, const s_pars, const s_grid*);

inline double solubility_in_gpl_from_ASconc_in_gpl(const string prname,
                                                   const string asname,
                                                   const double cl_AS,
                                                   const double Temp,
                                                   const double pH)
{
    if(prname=="lysozyme" && asname=="nacl") {
      
      return(346.01309430864*exp(-0.102667383846361*cl_AS));
      
    }
    else
    {
      cout << "Solubility not parameterized" << endl;
      exit(1);
    }
}


inline double volume_in_ml_from_masses_in_g(const string prname,
                                            const string asname,
                                            const double ml_protein,
                                            const double ml_AS,
                                            const double ml_water)
{
    double mtot = ml_protein + ml_AS + ml_water;
    double rho, wp;

    if(prname=="lysozyme" && asname=="nacl") {
      return((mtot+sqrt(mtot*mtot - 4.0*mtot*(0.65*ml_AS+0.303*ml_protein)))/2.0);
    }
    else if(prname=="ifng" && asname=="nh42so4") {
      wp = 100.0*ml_AS/ml_water;
      rho = 1.0+0.00575*wp;
      return(mtot/rho);
    }
    else
    {
      cout << "Specific volume of not parameterized" << endl;
      exit(1);
    }
}

#endif

