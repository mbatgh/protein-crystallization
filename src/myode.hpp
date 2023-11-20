

typedef runge_kutta_dopri5< state_type > stepper_type;

struct population
{
    const s_pars &m_pars;
    const s_grid *m_grid;
          s_idat &m_idat;

    population( const s_pars& pars, const s_grid* grid, s_idat& idat ) : m_pars (pars), m_grid(grid), m_idat(idat) {}

    void operator()( const state_type &n, state_type &dndt , double t ) const
    {
      double ml_protein,
             ml_AS,
             ml_water;

      double V_c,
             cl_protein,
	     cl_AS,
             V_liq,
             solubility,
	     supsat,
             dV_cryst,
             rGro,
             rNuc,
	     coll,
	     cfac,
	     stmp,
             cfr;
	     
      double rAggrBirth[MAX_N_C];
      double rAggrDeath[MAX_N_C];

      int    m, i, j, k, l;

      m = m_pars.NClasses_max-1;
      
      ml_protein = n[m+1];
      ml_AS      = n[m+2];
      ml_water   = n[m+3];

      V_liq          = volume_in_ml_from_masses_in_g(m_pars.name,m_pars.asname,ml_protein,ml_AS,ml_water);
      cl_protein     = 1000.0*ml_protein/V_liq;
      cl_AS          = 1000.0*ml_AS/V_liq;
     
      if(ml_protein < m_pars.ml_protein_0 * m_pars.pcp0/100.0 && tcp0==0.0) tcp0 = t;
      
      if(m_pars.name=="lysozyme"  && m_pars.asname=="nacl")     {
        
        solubility = 346.01309430864*exp(-0.102667383846361*cl_AS);
        
      }
      else if(m_pars.name=="ifng" && m_pars.asname=="nh42so4")  {
        
        stmp=100.0*ml_AS/ml_water;
        if(stmp<1.0) solubility=1000.0;
        else solubility = 2.0e+25 * pow(stmp,-17.6698);
        
      }
      else { fprintf(stderr, "solubility not parameterized\n"); exit(1);}
     


      if(solubility < cl_protein) supsat = (cl_protein - solubility)/solubility;
      else                        supsat =  0.0;
      
      if(m_pars.const_supsat > 0.0) supsat = m_pars.const_supsat;

      rGro = m_pars.k_growth_a * pow(supsat,m_pars.k_growth_b);                 
      rNuc = m_pars.k_nucl_a   * pow(supsat,m_pars.k_nucl_b);                   



      m = m_pars.NClasses_max-1;

      dndt[0] =   - rGro * (n[0] - 0     )/(m_grid[0].avg - 0               ) + rNuc * V_liq / m_grid[0].width;
      dndt[m] =   - rGro * (n[m] - n[m-1])/(m_grid[m].avg - m_grid[m-1].avg );

      for(i = 1; i < m; i++)
      {
        dndt[i] = -rGro * (n[i] - n[i-1])/(m_grid[i].avg - m_grid[i-1].avg );
      }



      if(m_pars.at>0) {
        if(m_pars.at<5) {

          for(i = 0; i <=m; i++) {
            rAggrBirth[i]=0.0;
            rAggrDeath[i]=0.0;
          }

          for(i = 0; i <=m; i++) {
            for(j = 0; j <=m; j++) { 
              coll = nakern[i][j] * n[i] * n[j] * m_grid[i].width * m_grid[j].width;
              rAggrDeath[i] +=  coll;
            }
          }
      
          for(i = 0; i <=m; i++) {
            for(l=1; l<=naggr[i]; l++) {
              k = andx1[i][l];
              j = andx2[i][l];
              coll = nakern[k][j] * n[k] * m_grid[k].width * n[j] * m_grid[j].width;
              rAggrBirth[i]   +=   ahw1[i][l] * coll;
              rAggrBirth[i+1] +=   ahw2[i][l] * coll;
            }
          }

          for(i = 0; i <= m; i++) {
            dndt[i] += (rAggrBirth[i] - rAggrDeath[i]) / m_grid[i].width;
          }

	} else {

	  cfac = supsat*supsat;

	  for(i = 0; i <=m; i++) {
            rAggrBirth[i]=0.0;
            rAggrDeath[i]=0.0;
          }

          for(i = 0; i <=m; i++) {
            for(j = 0; j <=m; j++) { 
              coll = cfac * nakern[i][j] * n[i] * n[j] * m_grid[i].width * m_grid[j].width;
              rAggrDeath[i] +=  coll / m_grid[i].width;
            }
          }
      
          for(i = 0; i <=m; i++) {
            for(l=1; l<=naggr[i]; l++) {
              k = andx1[i][l];
              j = andx2[i][l];
              coll = cfac * nakern[k][j] * n[k] * m_grid[k].width * n[j] * m_grid[j].width;
              rAggrBirth[i]   +=   ahw1[i][l] * coll / (m_grid[i].width);
              rAggrBirth[i+1] +=   ahw2[i][l] * coll / (m_grid[i+1].width);
            }
          }

          for(i = 0; i <= m; i++) {
            dndt[i] += rAggrBirth[i] - rAggrDeath[i];
          }
	}
      }
      


      dndt[m]   = 0.0;
      dndt[m-1] = 0.0;



      dV_cryst = 0.0;
      V_c = 0.0;
      
      for(i = 0; i <= m; i++)
      {
        dV_cryst += rGro * n[i] * m_grid[i].width * m_grid[i].avg * m_grid[i].avg  * 3.0;
        V_c      +=        n[i] * m_grid[i].width * m_grid[i].avg * m_grid[i].avg * m_grid[i].avg;
      }

      dV_cryst += rNuc * V_liq * pow(m_pars.SizeNucl,3.0);

      dndt[m+1] = -dV_cryst * m_pars.ds_protein * m_pars.shapefac;
      dndt[m+3] = -dV_cryst * m_pars.ds_water   * m_pars.shapefac;
      
      if(ml_AS<m_pars.ml_AS_max) {
        
        cfr = m_pars.fr_AS;
        if(t>m_pars.t1) cfr=m_pars.fr_AS2;
        if(t>m_pars.t2) cfr=m_pars.fr_AS3;



  
        dndt[m+2] =  cfr * m_pars.dl_AS_res;
        dndt[m+3] += cfr * m_pars.dl_water_res;
      }
      else
      {
        cfr = 0.0;
        dndt[m+2] = 0.0;
      }



      m_idat.csdoft     = m_pars.printlevel;
      m_idat.cl_protein = cl_protein;
      m_idat.V_c        = V_c;
      m_idat.V_liq      = V_liq;
      m_idat.solubility = solubility;
      m_idat.supsat     = supsat;
      m_idat.rGro       = rGro;
      m_idat.rNuc       = rNuc;
      m_idat.cfr        = cfr;
    }
};

struct streaming_observer
{
    std::ostream &m_outdat;
    std::ostream &m_outcsd;
    const s_grid *m_grid;
    s_idat       &m_idat;
    streaming_observer( std::ostream &outdat, std::ostream &outcsd, const s_grid* grid, s_idat& idat ) : 
                        m_outdat( outdat ),
                        m_outcsd( outcsd ),
                        m_grid  ( grid   ),
                        m_idat  ( idat   ) {}
    void operator()( const state_type &n , double t ) const
    {
      m_outdat << setw(16) << scientific << setprecision(6) << t;
      m_outdat << setw(16) << scientific << setprecision(6) << n[n.size()-3];         
      m_outdat << setw(16) << scientific << setprecision(6) << n[n.size()-2];         
      m_outdat << setw(16) << scientific << setprecision(6) << n[n.size()-1];         
      m_outdat << setw(16) << scientific << setprecision(6) << m_idat.cl_protein;     
      m_outdat << setw(16) << scientific << setprecision(6) << m_idat.V_c;            
      m_outdat << setw(16) << scientific << setprecision(6) << m_idat.V_liq;          
      m_outdat << setw(16) << scientific << setprecision(6) << m_idat.solubility;     
      m_outdat << setw(16) << scientific << setprecision(6) << m_idat.supsat;         
      m_outdat << setw(16) << scientific << setprecision(6) << m_idat.rGro;           
      m_outdat << setw(16) << scientific << setprecision(6) << m_idat.rNuc;           
      m_outdat << setw(16) << scientific << setprecision(6) << m_idat.cfr;            
      m_outdat << endl;

      if(m_idat.csdoft==2) {
        for(size_t u=0; u<(n.size())-3; u++) {
            m_outcsd << setw(14) << scientific << setprecision(6) << t 
                     << setw(14) << scientific << setprecision(6) << 10000.0*m_grid[u].avg
                     << setw(14) << scientific << setprecision(6) << n[u]*m_grid[u].width << endl;
        }
        m_outcsd << endl;
      }
    }
};

