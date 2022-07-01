/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Version Sep/22/2014
   Zhenxing Wang (KU)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(conp,FixConp)

#else

#ifndef LMP_FIX_CONP_H
#define LMP_FIX_CONP_H

#include "fix.h"
#include "pppm.h"
#include "ewald.h"
#include "pair_lj_cut_coul_long.h"
#include "pair.h"
#include "pair_coul_long.h"
#include "input.h"
#include "comm.h"
#include "neighbor.h"
#include "memory.h"
#include "neigh_request.h"


namespace LAMMPS_NS {

class FixConp : public Fix {
 public:
  FixConp(class LAMMPS *, int, char **);
  ~FixConp();
  int setmask();
  void init();
  void setup(int);
  void pre_force(int);
  void force_cal(int);
  void a_cal();
  void a_read();
  void b_cal();
  void equation_solve();
  void update_charge();
  int electrode_check(int);
  void sincos_a(double **);
  void sincos_b();
  void cg();
  void inv();
  void coul_cal(int, double *,int *);
  void pot_wall_wall();
  void pot_wall_wall_2();
  void pot_wall_ions();
  void GQ_H_call();
  void update_charge_2();
  void force_extra_ions_1();
  void force_extra_ions_2();
  void V_cal(); //function that calculates voltage

 private:
  int me,runstage;
  double Btime,Btime1,Btime2;
  double Ctime,Ctime1,Ctime2;
  double Ktime,Ktime1,Ktime2;
  double cgtime,cgtime1,cgtime2;
  FILE *outf,*outa,*a_matrix_fp;
  int a_matrix_f;
  int minimizer;
  double vL,vR;
  int molidL,molidR;
  int maxiter;
  double tolerance;
  double C_penalty;       //RS on 7-1-2022: constant of added penalty function

  double rms(int,double,bigint,double);
  void coeffs();

  double unitk[3];
  double *ug;
  double g_ewald,    gsqmx,volume,slab_volfactor;	//RS on 21-4-2021: changed to eliminate eta
//double g_ewald,eta,gsqmx,volume,slab_volfactor;
  int *kxvecs,*kyvecs,*kzvecs;
  double ***cs,***sn,**csk,**snk;
  int kmax,kmax3d,kmax_created,kcount;
  int kxmax,kymax,kzmax;
  double *sfacrl,*sfacrl_all,*sfacim,*sfacim_all;
  int everynum;
  int elenum,elenum_old,elenum_all;
  double *eleallq;
  double *aaa_all,*bbb_all;
  double *sss_all;
  int *tag2eleall,*eleall2tag,*curr_tag2eleall,*ele2tag;

  double *Psi_w_w;
  double *Psi_w_wi;
  double *Psi_w_i;
  double *Bq_CP4M2;
  double *V_min_Bq_CP4M2;

  //Needed for constant charge, declared as global because might cause segmentation fault
  double *O_1;
  double *O_2;
  double *OAinv_1;
  double *OAinv_2;
  double *Z;

  const double q_L = -0.001;    //charges used for wall particles for the calculation of the potential on the wall particles due to the ions
  const double q_R = 0.001;
  const double conv = 4186.8/(6.02214e23*1.60218e-19);
  const double evscale = 0.069447;
  double *GQ;   //RS on 8-6-2022: for constant charge
  double *H;    //RS on 8-6-2022:for constant charge
  double C_pp;
  double C_pm;
  double C_mp;
  double C_mm;

  PPPM obj_kspace = PPPM(lmp);
  //Ewald obj_kspace = Ewald(lmp);
  PairCoulLong obj_CoulLong = PairCoulLong(lmp);  //declaration of an instance of the class PairCoulLong
};

};


#endif
#endif
