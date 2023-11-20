
#include "aspbe.hpp"
#include <stdio.h>
#include <iostream>
#include <vector>
#include "spline.hpp"

s_dprop ana_cld(state_type& cld, const state_type n, const s_pars my_pars, const s_grid *my_grid)
{

    using namespace magnet::math;
    
    int i, j;
    int nplane;
    int MISSED;
    int GOTONEXT;
    int isiz;
    double hb, maxx, maxy, crosss, D, N, tE, tL, t, xleft, xright;
    
    Vector3d c[8];
    Vector3d V[6];
    Vector3d nn[6];
    Vector3d P0,P1,v1,v2,dP,dS,cs1,cs2;

    Spline spline;

    int imod, i_bot_left, i_bot_right;
    s_dprop dprop;
    state_type nnv(my_pars.NClasses_max, 0.0 );
    double y,xmax,ymax,vtot,ntot;
    
    dprop.mode = 0.0;
    dprop.fwhm = 0.0;
    dprop.pvat = 0.0;

    xmax = 0.0;
    ymax = 0.0;
    imod = 0;  

    srand(my_pars.rseed);

    for(i=0;i<my_pars.NClasses_max;i++) cld[i]=0.0;
    
    for(i=0;i<my_pars.NClasses_max;i++) {

      if(n[i]>0.0) {
       
        hb=my_grid[i].avg*0.5;
        maxx=sqrt(3.0)*my_grid[i].avg;
        maxy=maxx;

        for(j=1;j<=my_pars.nrot;j++) {
          
          c[0] << -hb,-hb, hb;
          c[1] << -hb, hb, hb;
          c[2] <<  hb,-hb, hb;
          c[3] <<  hb, hb, hb;
          c[4] << -hb,-hb,-hb;
          c[5] << -hb, hb,-hb;
          c[6] <<  hb,-hb,-hb;
          c[7] <<  hb, hb,-hb;  
          
          randrot(c, 8);
          
          v1 = c[0]-c[1];
          v2 = c[3]-c[1];
          nn[0]=v1.cross(v2);
          v1 = c[6]-c[7];
          v2 = c[5]-c[7];
          nn[1]=v1.cross(v2);
          v1 = c[1]-c[5];
          v2 = c[7]-c[5];
          nn[2]=v1.cross(v2);
          v1 = c[2]-c[6];
          v2 = c[4]-c[6];
          nn[3]=v1.cross(v2);
          v1 = c[3]-c[7];
          v2 = c[6]-c[7];
          nn[4]=v1.cross(v2);
          v1 = c[0]-c[4];
          v2 = c[5]-c[4];
          nn[5]=v1.cross(v2);
          
          tE=0.0;
          tL=1.0;
          
          P0 << maxx*(randd()-0.5), maxy*(randd()-0.5), 2.0*my_grid[i].avg;
          dP << 0,0,-4.0*my_grid[i].avg;
          P1 = P0 + dP;
          dS = P1 - P0;
          
          V[0]=c[1];
          V[1]=c[7];
          V[2]=c[5];
          V[3]=c[6];
          V[4]=c[7];
          V[5]=c[4];
          
          MISSED=0;

          for(nplane=0;nplane<6;nplane++) {
            
            GOTONEXT = 0;
            v1 = P0 - V[nplane];
            N = -1.0 * (v1.dot(nn[nplane]));
            D = dS.dot(nn[nplane]);

            if(D==0.0) {
              if(N<0.0) {
                MISSED = 1;
              } else {
                GOTONEXT = 1;
              }
            } else {
              if(GOTONEXT==0) {
                
                t = N/D;
          
                if(D<0.0) {
                  if(t>tE) tE = t;
                  if(tE>tL) MISSED = 1;
                }
            
                if(D>0.0) {
                  if(t<tL) tL=t;
                  if(tL<tE) MISSED = 1;
                } 
                
              }
            }
          }

          if(MISSED==0) {
            
            cs2 = P0 + tE * dS;
            cs1 = P0 + tL * dS;
            v1 =  cs2-cs1;
            crosss = v1.norm();
            isiz=0;
            
            while(my_grid[isiz].avg<crosss && isiz<my_pars.NClasses_max-1) isiz++;
            
            cld[isiz] += n[i];
          } 
        }
      } 
    }
    
  
  
    ntot=0.0;
    for(i=0;i<my_pars.NClasses_max;i++) {
      nnv[i] = cld[i]*pow(my_grid[i].avg,2.0);
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


void randrot(Vector3d* c, int n) {

   int i;
   double phi, Phi, xi;
   double cos_phi, sin_phi, cos_Phi, sin_Phi, eta;
   double u, v, w, x, y, z;
   double a[3][3];

   phi=2.*M_PI*randd();
   Phi=2.*M_PI*randd();
   xi=2.*randd()-1.;

   cos_phi=cos(phi); sin_phi=sin(phi);
   cos_Phi=cos(Phi); sin_Phi=sin(Phi);
   eta=sqrt(1.-xi*xi);

   a[0][0]=cos_Phi*cos_phi-xi*sin_Phi*sin_phi;
   a[0][1]=cos_Phi*sin_phi+xi*sin_Phi*cos_phi;
   a[0][2]=eta*sin_Phi;
   a[1][0]= -sin_Phi*cos_phi-xi*cos_Phi*sin_phi;
   a[1][1]= -sin_Phi*sin_phi+xi*cos_Phi*cos_phi;
   a[1][2]=eta*cos_Phi;
   a[2][0]=eta*sin_phi;
   a[2][1]= -eta*cos_phi;
   a[2][2]=xi;

   for(i=0;i<n;i++)
   {  u=c[i](0); v=c[i](1); w=c[i](2);
      x=a[0][0]*u+a[0][1]*v+a[0][2]*w;
      y=a[1][0]*u+a[1][1]*v+a[1][2]*w;
      z=a[2][0]*u+a[2][1]*v+a[2][2]*w;
      c[i] << x,y,z;
   }
}


double randd(void) {
  return (double)random()/(double)RAND_MAX;
}

