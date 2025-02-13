/* ---------------------------------------------------------------------
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
   Version: Sep/22/2014
   Zhenxing Wang(KU)
------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "stddef.h"
#include "fix_conp.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "force.h"
#include "respa.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "compute.h"

#include "pair.h"
#include "kspace.h"
#include "comm.h"
#include "mpi.h"
#include "math_const.h"
#include "neigh_list.h"
#include "domain.h"
#include "iostream"

#define EWALD_F   1.12837917
#define EWALD_P   0.3275911
#define A1        0.254829592
#define A2       -0.284496736
#define A3        1.421413741
#define A4       -1.453152027
#define A5        1.061405429

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

enum{CONSTANT,EQUAL,ATOM};

extern "C" {
  void dgetrf_(const int *M,const int *N,double *A,const int *lda,int *ipiv,int *info);
  void dgetri_(const int *N,double *A,const int *lda,const int *ipiv,double *work,const int *lwork,int *info);
}
/* ---------------------------------------------------------------------- */
//RS on 15-07-2022: the fix command has been adjusted to
// fix [ID] all conp [Nevery] [penalty] [Molecule-ID 1] [Molecule-ID 2] [Method] [Input arg 1] [Input arg 2] [Minimization method] [Log] [Matrix]
// for [Method] use the word potential (for constant potential simulation) or charge (for constant charge, charging or discharging)
// for simulating constant potential: [Input arg 1] = potential left, [Input arg 2] = potential right (in Volt)
// for simulating constant charge: [Input arg 1] = dummy value, [Input arg 2] = 0
// for simulating charging: [Input arg 1] = potential of source (in Volt), [Input arg 2] = Resistance (in fs*V)/e
// for simulating discharging: [Input arg 1] = 0, [Input arg 2] = Resistance (in fs*V)/e

FixConp::FixConp(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 11) error->all(FLERR,"Illegal fix conp command");
  maxiter = 100;
  tolerance = 0.000001;
  everynum = utils::numeric(FLERR,arg[3],false,lmp);
  C_penalty = utils::numeric(FLERR,arg[4],false,lmp);	//RS: reads in constant of penalty function instead of eta as in old code
  molidL = utils::inumeric(FLERR,arg[5],false,lmp);
  molidR = utils::inumeric(FLERR,arg[6],false,lmp);
  if (strcmp(arg[7],"potential") == 0){
      method = 0;
      vL = evscale*utils::numeric(FLERR,arg[8],false,lmp); //RS: voltage unit in transformed from V to e/Angstrom by multiplying by evscale (1 over Coulomb constant),
      vR = evscale*utils::numeric(FLERR,arg[9],false,lmp);
  }
  else if (strcmp(arg[7],"charge") == 0){
      method = 1;
      Vs = evscale*utils::numeric(FLERR,arg[8],false,lmp);
      R = evscale*utils::numeric(FLERR,arg[9],false,lmp);   //RS: as the unit of voltage is transformed, also the unit of R is transformed
  }
  else {
      error->all(FLERR,"invalid method");
  }
  if (strcmp(arg[10],"cg") == 0) {
    minimizer = 0;
  } else if (strcmp(arg[10],"inv") == 0) {
    minimizer = 1;
  } else error->all(FLERR,"Unknown minimization method");

  outf = fopen(arg[11],"w");
  if (narg == 13) {
    outa = NULL;
    a_matrix_fp = fopen(arg[12],"r");
    if (a_matrix_fp == NULL) error->all(FLERR,"Cannot open A matrix file");
    if (strcmp(arg[12],"org") == 0) {
      a_matrix_f = 1;
    } else if (strcmp(arg[12],"inv") == 0) {
      a_matrix_f = 2;
    } else {
      error->all(FLERR,"Unknown A matrix type");
    }
  } else {
    a_matrix_f = 0;
  }
  elenum = elenum_old = 0;
  csk = snk = NULL;
  aaa_all = NULL;
  bbb_all = NULL;
  Psi_w_w = NULL;
  Psi_w_i = NULL;
  Psi_w_wi = NULL;
  Bq_CP4M2 = NULL;
  V_min_Bq_CP4M2 = NULL;

  //needed for constant charge
  sss_all = O_1 = O_2 = OAinv_1 = OAinv_2 = G = H = NULL;

  tag2eleall = eleall2tag = curr_tag2eleall = ele2tag = NULL;
  Btime = cgtime = Ctime = Ktime = 0;
  runstage = 0; //after operation
                //0:init; 1: a_cal; 2: first sin/cos cal; 3: inv only, aaa inverse
}

/* ---------------------------------------------------------------------- */

FixConp::~FixConp()
{
  fclose(outf);
  memory->destroy3d_offset(cs,-kmax_created);
  memory->destroy3d_offset(sn,-kmax_created);
  delete [] aaa_all;
  delete [] bbb_all;
  delete [] curr_tag2eleall;
  delete [] tag2eleall;
  delete [] eleall2tag;
  delete [] ele2tag;
  delete [] kxvecs;
  delete [] kyvecs;
  delete [] kzvecs;
  delete [] ug;
  delete [] sfacrl;
  delete [] sfacim;
  delete [] sfacrl_all;
  delete [] sfacim_all;
  delete [] Psi_w_w;
  delete [] Psi_w_i;
  delete [] Psi_w_wi;
  delete [] Bq_CP4M2;
  delete [] V_min_Bq_CP4M2;

  delete [] sss_all;
  delete [] O_1;
  delete [] O_2;
  delete [] OAinv_1;
  delete [] OAinv_2;
  delete [] Z;
  delete [] G;
  delete [] H;

}

/* ---------------------------------------------------------------------- */

int FixConp::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixConp::init()
{
  MPI_Comm_rank(world,&me);
}

/* ---------------------------------------------------------------------- */

void FixConp::setup(int vflag)
{
  g_ewald = force->kspace->g_ewald;
  slab_volfactor = force->kspace->slab_volfactor;
  double accuracy = force->kspace->accuracy;

  int i;
  double qsqsum = 0.0;
  for (i = 0; i < atom->nlocal; i++) {
    qsqsum += atom->q[i]*atom->q[i];
  }
  double tmp,q2;
  MPI_Allreduce(&qsqsum,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
  qsqsum = tmp;
  q2 = qsqsum * force->qqrd2e / force->dielectric;

// Copied from ewald.cpp
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  double zprd_slab = zprd*slab_volfactor;
  volume = xprd * yprd * zprd_slab;

  unitk[0] = 2.0*MY_PI/xprd;
  unitk[1] = 2.0*MY_PI/yprd;
  unitk[2] = 2.0*MY_PI/zprd_slab;

  bigint natoms = atom->natoms;
  double err;
  kxmax = 1;
  kymax = 1;
  kzmax = 1;

  err = rms(kxmax,xprd,natoms,q2);
  while (err > accuracy) {
    kxmax++;
    err = rms(kxmax,xprd,natoms,q2);
  }


  err = rms(kymax,yprd,natoms,q2);
  while (err > accuracy) {
    kymax++;
    err = rms(kymax,yprd,natoms,q2);
  }

  err = rms(kzmax,zprd_slab,natoms,q2);
  while (err > accuracy) {
    kzmax++;
    err = rms(kzmax,zprd_slab,natoms,q2);
  }

  //RS: code added start
  //to overrule kmax, are different for different initial charges, which results in different amatrix, the effect is especially large when ions are removed from box
    //kxmax = 6;
    //kymax = 6;
    //kzmax = 12;

    FILE *out_k_values = fopen("k_values", "a");
    fprintf(out_k_values,"%20d %20d %20d\n", kxmax, kymax, kzmax);
    fclose(out_k_values);

    //RS: code added end

  kmax = MAX(kxmax,kymax);
  kmax = MAX(kmax,kzmax);
  kmax3d = 4*kmax*kmax*kmax + 6*kmax*kmax + 3*kmax;

  kxvecs = new int[kmax3d];
  kyvecs = new int[kmax3d];
  kzvecs = new int[kmax3d];
  ug = new double[kmax3d];

  double gsqxmx = unitk[0]*unitk[0]*kxmax*kxmax;
  double gsqymx = unitk[1]*unitk[1]*kymax*kymax;
  double gsqzmx = unitk[2]*unitk[2]*kzmax*kzmax;
  gsqmx = MAX(gsqxmx,gsqymx);
  gsqmx = MAX(gsqmx,gsqzmx);

  gsqmx *= 1.00001;

  coeffs();
  kmax_created = kmax;

//copied from ewald.cpp end

  int nmax = atom->nmax;
  //double evscale = 0.069447; //RS: now declared as global constant, so it can be used in other functions
  //if (method == 0){           //RS: now just done as soon as the input aruments are read
  //    vL *= evscale;
  //    vR *= evscale;
  //}

  memory->create3d_offset(cs,-kmax,kmax,3,nmax,"fixconp:cs");
  memory->create3d_offset(sn,-kmax,kmax,3,nmax,"fixconp:sn");
  sfacrl = new double[kmax3d];
  sfacim = new double[kmax3d];
  sfacrl_all = new double[kmax3d];
  sfacim_all = new double[kmax3d];
  tag2eleall = new int[natoms+1];
  curr_tag2eleall = new int[natoms+1];
    if (runstage == 0) {
    int i;
    int nlocal = atom->nlocal;
    for ( i = 0; i < nlocal; i++) {
      if (electrode_check(i)) ++elenum;
    }
    MPI_Allreduce(&elenum,&elenum_all,1,MPI_INT,MPI_SUM,world);

    eleall2tag = new int[elenum_all];
    aaa_all = new double[elenum_all*elenum_all];
    bbb_all = new double[elenum_all];
    ele2tag = new int[elenum];
    for (i = 0; i < natoms+1; i++) tag2eleall[i] = -1;
    for (i = 0; i < natoms+1; i++) curr_tag2eleall[i] = -1;
    if (minimizer == 0) {
      eleallq = new double [elenum_all];
    }
    if (a_matrix_f == 0) {
      if (me == 0) outa = fopen("amatrix","w");
      a_cal();
    } else {
      a_read();
    }

    Psi_w_w = new double[elenum_all];
    Psi_w_wi = new double[elenum_all];
    Psi_w_i = new double[elenum_all];
    Bq_CP4M2 = new double[elenum_all];
    V_min_Bq_CP4M2 = new double[elenum_all];

    sss_all = new double[elenum_all*elenum_all];
    O_1 = new double[elenum_all*elenum_all];
    O_2 = new double[elenum_all*elenum_all];
    OAinv_1 = new double[elenum_all*elenum_all];
    OAinv_2 = new double[elenum_all*elenum_all];
    Z = new double[elenum_all*elenum_all];
    G = new double[elenum_all];
    H = new double[elenum_all*elenum_all];

    pot_wall_wall();
    //pot_wall_wall_2();

    runstage = 1;
    }

}

/* ---------------------------------------------------------------------- */

void FixConp::pre_force(int vflag)
{
    //RS: for printing start
    double *q = atom->q;
    int *tag = atom-> tag;
    int nlocal = atom->nlocal;
    //RS: for printing end

  if(update->ntimestep % everynum == 0) {
    if (strstr(update->integrate_style,"verlet")) { //not respa
      Btime1 = MPI_Wtime();

      //FILE *out_q_before_b_cal = fopen("q_before_b_cal", "a");
      //for (int i = 0; i < nlocal; i++) {
      //    fprintf(out_q_before_b_cal,"%20d %20f\n", tag[i], q[i]);
      //}
      //fclose(out_q_before_b_cal);
      //b_cal();
      pot_wall_ions();
      //FILE *out_q_after_b_cal = fopen("q_after_b_cal", "a");
      //for (int i = 0; i < nlocal; i++) {
      //    fprintf(out_q_after_b_cal,"%20d %20f\n", tag[i], q[i]);
      //}
      //fclose(out_q_after_b_cal);


      Btime2 = MPI_Wtime();
      Btime += Btime2-Btime1;
      if (update->laststep == update->ntimestep) {
        double Btime_all;
        MPI_Reduce(&Btime,&Btime_all,1,MPI_DOUBLE,MPI_SUM,0,world);
        double Ctime_all;
        MPI_Reduce(&Ctime,&Ctime_all,1,MPI_DOUBLE,MPI_SUM,0,world);
        double Ktime_all;
        MPI_Reduce(&Ktime,&Ktime_all,1,MPI_DOUBLE,MPI_SUM,0,world);
        if (me == 0) {
          Btime = Btime_all/comm->nprocs;
          Ctime = Ctime_all/comm->nprocs;
          Ktime = Ktime_all/comm->nprocs;
          fprintf(outf,"B vector calculation time = %g\n",Btime);
          fprintf(outf,"Coulomb calculation time = %g\n",Ctime);
          fprintf(outf,"Kspace calculation time = %g\n",Ktime);
        }
      }
    }
      //Psi_w_wi = new double[elenum_all];
      //Psi_w_i = new double[elenum_all];

    equation_solve();

      //FILE *out_q_before_update_charge = fopen("q_before_update_charge", "a");
      //for (int i = 0; i < nlocal; i++) {
      //    fprintf(out_q_before_update_charge,"%20d %20f\n", tag[i], q[i]);
      //}
      //fclose(out_q_before_update_charge);

      if (method == 0){
          update_charge();    //for constant potential
      }
      else{
          update_charge_2();    //for constant charge
      }

      //FILE *out_q_after_update_charge = fopen("q_after_update_charge", "a");
      //for (int i = 0; i < nlocal; i++) {
          //fprintf(out_q_after_update_charge,"%20d %20f\n", tag[i], q[i]);
      //}
      //fclose(out_q_after_update_charge);

  }
  //force_cal(vflag);		RS on 22-04-2021: commented out to remove forces due to erfc(eta*rij) terms and the added energy (eta/sqrt(2pi))*sum(Q_i^2), later in the code the function force_cal is also commented out.
  //force_extra_ions_1();       //RS: the extra forces that we calculated initially, but since they add up to zero, we don't include them anymore
  //force_extra_ions_2();
}

/* ---------------------------------------------------------------------- */

int FixConp::electrode_check(int atomid)
{
  int *molid = atom->molecule;
  if (molid[atomid] == molidL) return 1;
  else if (molid[atomid] == molidR) return -1;
  else return 0;
}

/* ----------------------------------------------------------------------*/

void FixConp::b_cal()
{
  Ktime1 = MPI_Wtime();
  int i,j,k;
  int nmax = atom->nmax;
  if (atom->nlocal > nmax) {
    memory->destroy3d_offset(cs,-kmax_created);
    memory->destroy3d_offset(sn,-kmax_created);
    nmax = atom->nmax;
    kmax_created = kmax;
  }
  sincos_b();
  MPI_Allreduce(sfacrl,sfacrl_all,kcount,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(sfacim,sfacim_all,kcount,MPI_DOUBLE,MPI_SUM,world);
  double **x = atom->x;
  double *q = atom->q;
  int *tag = atom->tag;
  int nlocal = atom->nlocal;
  int kx,ky,kz;
  double cypz,sypz,exprl,expim,kspacetmp;
  elenum = 0;
  for (i = 0; i < nlocal; i++) {
    if(electrode_check(i)) elenum++;
  }
  double bbb[elenum];
  j = 0;
  for (i = 0; i < nlocal; i++) {
    if (electrode_check(i) == 1) {
      bbb[j] = vL;
      j++;
    }
    if (electrode_check(i) == -1) {
      bbb[j] = vR;
      j++;
    }
  }
  for (k = 0; k < kcount; k++) {
    kx = kxvecs[k];
    ky = kyvecs[k];
    kz = kzvecs[k];
    j = 0;
    for (i = 0; i < nlocal; i++) {
      if (electrode_check(i)) {
        cypz = cs[ky][1][i]*cs[kz][2][i] - sn[ky][1][i]*sn[kz][2][i];
        sypz = sn[ky][1][i]*cs[kz][2][i] + cs[ky][1][i]*sn[kz][2][i];
        exprl = cs[kx][0][i]*cypz - sn[kx][0][i]*sypz;
        expim = sn[kx][0][i]*cypz + cs[kx][0][i]*sypz;
        bbb[j] -= 2.0*ug[k]*(exprl*sfacrl_all[k]+expim*sfacim_all[k]);
        j++;
      }
    }
  }

  //elenum_list and displs for gathering ele tag list and bbb
  int nprocs = comm->nprocs;
  int elenum_list[nprocs];
  MPI_Allgather(&elenum,1,MPI_INT,elenum_list,1,MPI_INT,world);
  int displs[nprocs];
  displs[0] = 0;
  int displssum = 0;
  for (i = 1; i < nprocs; ++i) {
    displssum += elenum_list[i-1];
    displs[i] = displssum;
  }

  //slabcorrection and create ele tag list in current timestep
  double slabcorrtmp = 0.0;
  double slabcorrtmp_all = 0.0;
  for (i = 0; i < nlocal; i++) {
    if (electrode_check(i) == 0) {
      slabcorrtmp += 4*q[i]*MY_PI*x[i][2]/volume;
    }
  }
  MPI_Allreduce(&slabcorrtmp,&slabcorrtmp_all,1,MPI_DOUBLE,MPI_SUM,world);
  j = 0;
  for (i = 0; i < nlocal; i++) {
    if (electrode_check(i)) {
      bbb[j] -= x[i][2]*slabcorrtmp_all;
      ele2tag[j] = tag[i];
      j++;
    }
  }
  Ktime2 = MPI_Wtime();
  Ktime += Ktime2-Ktime1;

  coul_cal(1,bbb,ele2tag);



  //gather ele tag list
  int ele_taglist_all[elenum_all];
  int tagi;
  MPI_Allgatherv(ele2tag,elenum,MPI_INT,&ele_taglist_all,elenum_list,displs,MPI_INT,world);
  for (i = 0; i < elenum_all; i++) {
    tagi = ele_taglist_all[i];
    curr_tag2eleall[tagi] = i;
  }

  //gather b to bbb_all and sort in the same order as aaa_all
  double bbb_buf[elenum_all];
  MPI_Allgatherv(&bbb,elenum,MPI_DOUBLE,&bbb_buf,elenum_list,displs,MPI_DOUBLE,world);
  int elei;
  for (i = 0; i < elenum_all; i++) {
    tagi = eleall2tag[i];
    elei = curr_tag2eleall[tagi];
    bbb_all[i] = bbb_buf[elei];
  }
}

/*----------------------------------------------------------------------- */
void FixConp::equation_solve()
{
//solve equations
  if (minimizer == 0) {
    cgtime1 = MPI_Wtime();
    cg();
    cgtime2 = MPI_Wtime();
    cgtime += cgtime2-cgtime1;
    if (update->laststep == update->ntimestep) {
      double cgtime_all;
      MPI_Reduce(&cgtime,&cgtime_all,1,MPI_DOUBLE,MPI_SUM,0,world);
      if (me == 0) {
        cgtime = cgtime_all/comm->nprocs;
        if (screen) fprintf(screen,"conjugate gradient solver time = %g\n",cgtime);
        if (logfile) fprintf(logfile,"conjugate gradient solver time = %g\n",cgtime);
      }
    }
  } else if (minimizer == 1) {
    inv();
  }
}

/*----------------------------------------------------------------------- */
void FixConp::a_read()
{
  int i = 0;
  int idx1d;
  if (me == 0) {
    int maxchar = 21*elenum_all+1;
    char line[maxchar];
    char *word;
    while(fgets(line,maxchar,a_matrix_fp) != NULL) {
      word = strtok(line," \t");
      while(word != NULL) {
        if (i < elenum_all) {
          eleall2tag[i] = atoi(word);
        } else {
          idx1d = i-elenum_all;
          aaa_all[idx1d] = atof(word);
        }
        word = strtok(NULL," \t");
        i++;
      }
    }
    fclose(a_matrix_fp);
  }
  MPI_Bcast(eleall2tag,elenum_all,MPI_INT,0,world);
  MPI_Bcast(aaa_all,elenum_all*elenum_all,MPI_DOUBLE,0,world);

  int tagi;
  for (i = 0; i < elenum_all; i++) {
    tagi = eleall2tag[i];
    tag2eleall[tagi] = i;
  }
}

/*----------------------------------------------------------------------- */
void FixConp::a_cal()
{
  double t1,t2;
  t1 = MPI_Wtime();
  Ktime1 = MPI_Wtime();
  if (me == 0) {
    fprintf(outf,"A matrix calculating ...\n");
  }

  double **eleallx = NULL;
  memory->create(eleallx,elenum_all,3,"fixconp:eleallx");

  int nprocs = comm->nprocs;
  int nlocal = atom->nlocal;
  int *tag = atom->tag;
  int i,j,k;
  int elenum_list[nprocs];
  MPI_Allgather(&elenum,1,MPI_INT,elenum_list,1,MPI_INT,world);
  int displs[nprocs];
  displs[0] = 0;
  int displssum = 0;
  for (i = 1; i < nprocs; ++i) {
    displssum += elenum_list[i-1];
    displs[i] = displssum;
  }
  j = 0;
  for (i = 0; i < nlocal; i++) {
    if (electrode_check(i)) {
      ele2tag[j] = tag[i];
      j++;
    }
  }

  //gather tag,x and q
  double **x = atom->x;

  double *elexyzlist = new double[3*elenum];
  double *elexyzlist_all = new double[3*elenum_all];
  j = 0;
  for (i = 0; i < nlocal; i++) {
    if (electrode_check(i)) {
      elexyzlist[j] = x[i][0];
      j++;
      elexyzlist[j] = x[i][1];
      j++;
      elexyzlist[j] = x[i][2];
      j++;
    }
  }
  MPI_Allgatherv(ele2tag,elenum,MPI_INT,eleall2tag,elenum_list,displs,MPI_INT,world);
  int displs2[nprocs];
  int elenum_list2[nprocs];
  for (i = 0; i < nprocs; i++) {
    elenum_list2[i] = elenum_list[i]*3;
    displs2[i] = displs[i]*3;
  }
  MPI_Allgatherv(elexyzlist,elenum*3,MPI_DOUBLE,elexyzlist_all,elenum_list2,displs2,MPI_DOUBLE,world);

  double *aaa = new double[elenum*elenum_all];
  for (i = 0; i < elenum*elenum_all; i++) {
    aaa[i] = 0.0;
  }
  j = 0;
  for (i = 0; i < elenum_all; i++) {
    if (i == 0 && me == 0) fprintf(outa," ");
    if (me == 0) fprintf (outa,"%20d",eleall2tag[i]);
    tag2eleall[eleall2tag[i]] = i;
    eleallx[i][0] = elexyzlist_all[j];
    j++;
    eleallx[i][1] = elexyzlist_all[j];
    j++;
    eleallx[i][2] = elexyzlist_all[j];
    j++;
  }
  if (me == 0) fprintf (outa,"\n");


  memory->create(csk,kcount,elenum_all,"fixconp:csk");
  memory->create(snk,kcount,elenum_all,"fixconp:snk");
  sincos_a(eleallx);
  delete [] elexyzlist;
  delete [] elexyzlist_all;

  int elealli,elei,idx1d;
  double zi;
  double CON_4PIoverV = MY_4PI/volume;
  double CON_s2overPIS = sqrt(2.0)/MY_PIS;
  double CON_2overPIS = 2.0/MY_PIS;
  for (i = 0; i < nlocal; ++i) {
    zi = x[i][2];
    if (electrode_check(i)) {
      elealli = tag2eleall[tag[i]];
      for (k = 0; k < elenum; ++k) {
        if (ele2tag[k] == tag[i]) {
          elei = k;
          break;
        }
      }
      for (j = 0; j < elenum_all; ++j) {
        idx1d = elei*elenum_all+j;
        for (k = 0; k < kcount; ++k) {
          aaa[idx1d] += 2.0*ug[k]*(csk[k][elealli]*csk[k][j]+snk[k][elealli]*snk[k][j]);
        }
        aaa[idx1d] += CON_4PIoverV*zi*eleallx[j][2];
      }
      idx1d = elei*elenum_all+elealli;
      aaa[idx1d] +=                  -CON_2overPIS*g_ewald; //gaussian self correction		//RS on 22-04-2021: removed eta term in self correction term in the A-matrix
    //aaa[idx1d] += CON_s2overPIS*eta-CON_2overPIS*g_ewald; //gaussian self correction
      aaa[idx1d] += C_penalty;
    }
  }

  memory->destroy(eleallx);
  memory->destroy(csk);
  memory->destroy(snk);

  coul_cal(2,aaa,ele2tag);

  int elenum_list3[nprocs];
  int displs3[nprocs];
  for (i = 0; i < nprocs; i++) {
    elenum_list3[i] = elenum_list[i]*elenum_all;
    displs3[i] = displs[i]*elenum_all;
  }
  MPI_Allgatherv(aaa,elenum*elenum_all,MPI_DOUBLE,aaa_all,elenum_list3,displs3,MPI_DOUBLE,world);
  delete [] aaa;
  aaa = NULL;
  for (i = 0; i < elenum_all; ++i) {
    for (j = 0; j < elenum_all; ++j) {
      idx1d = i*elenum_all+j;
      if (j != 0 && me == 0) fprintf(outa," ");
      if (me == 0) fprintf (outa,"%20.12f",aaa_all[idx1d]);
    }
    if (me == 0) fprintf (outa,"\n");
  }
  if(me == 0) fclose(outa);

  t2 = MPI_Wtime();
  double tsum = t2 - t1;
  double tsum_all;
  MPI_Allreduce(&tsum,&tsum_all,1,MPI_DOUBLE,MPI_SUM,world);
  if (me == 0) {
    tsum = tsum_all/nprocs;
    fprintf(outf,"A matrix calculation time  = %g\n",tsum);
  }
  Ktime2 = MPI_Wtime();
  Ktime += Ktime2-Ktime1;
}
/*--------------------------------------------------------------*/

void FixConp::sincos_a(double **eleallx)
{
  int i,m,k,ic;
  int kx,ky,kz;
  double ***csele,***snele;
  memory->create3d_offset(csele,-kmax,kmax,3,elenum_all,"fixconp:csele");
  memory->create3d_offset(snele,-kmax,kmax,3,elenum_all,"fixconp:snele");
  double sqk,cypz,sypz;
  for (ic = 0; ic < 3; ic++) {
    sqk = unitk[ic]*unitk[ic];
    if (sqk <= gsqmx) {
      for (i = 0; i < elenum_all; i++) {
        csele[0][ic][i] = 1.0;
        snele[0][ic][i] = 0.0;
        csele[1][ic][i] = cos(unitk[ic]*eleallx[i][ic]);
        snele[1][ic][i] = sin(unitk[ic]*eleallx[i][ic]);
        csele[-1][ic][i] = csele[1][ic][i];
        snele[-1][ic][i] = -snele[1][ic][i];
      }
    }
  }

  for (m = 2; m <= kmax; m++) {
    for (ic = 0; ic < 3; ic++) {
      sqk = m*unitk[ic] * m*unitk[ic];
      if (sqk <= gsqmx) {
        for (i = 0; i < elenum_all; i++) {
          csele[m][ic][i] = csele[m-1][ic][i]*csele[1][ic][i] -
            snele[m-1][ic][i]*snele[1][ic][i];
          snele[m][ic][i] = snele[m-1][ic][i]*csele[1][ic][i] +
            csele[m-1][ic][i]*snele[1][ic][i];
          csele[-m][ic][i] = csele[m][ic][i];
          snele[-m][ic][i] = -snele[m][ic][i];
        }
      }
    }
  }
  for (k = 0; k < kcount; ++k) {
    kx = kxvecs[k];
    ky = kyvecs[k];
    kz = kzvecs[k];
    for (i = 0; i < elenum_all; ++i) {
      cypz = csele[ky][1][i]*csele[kz][2][i] - snele[ky][1][i]*snele[kz][2][i];
      sypz = snele[ky][1][i]*csele[kz][2][i] + csele[ky][1][i]*snele[kz][2][i];
      csk[k][i] = csele[kx][0][i]*cypz - snele[kx][0][i]*sypz;
      snk[k][i] = snele[kx][0][i]*cypz + csele[kx][0][i]*sypz;
    }
  }
  memory->destroy3d_offset(csele,-kmax_created);
  memory->destroy3d_offset(snele,-kmax_created);
}

/*--------------------------------------------------------------*/
void FixConp::sincos_b()
{
  int i,k,l,m,n,ic;
  double cstr1,sstr1,cstr2,sstr2,cstr3,sstr3,cstr4,sstr4;
  double sqk,clpm,slpm;

  double **x = atom->x;
  double *q = atom->q;
  int nlocal = atom->nlocal;

  n = 0;

  // (k,0,0), (0,l,0), (0,0,m)

  for (ic = 0; ic < 3; ic++) {
    sqk = unitk[ic]*unitk[ic];
    if (sqk <= gsqmx) {
      cstr1 = 0.0;
      sstr1 = 0.0;
      for (i = 0; i < nlocal; i++) {
          cs[0][ic][i] = 1.0;
          sn[0][ic][i] = 0.0;
          cs[1][ic][i] = cos(unitk[ic]*x[i][ic]);
          sn[1][ic][i] = sin(unitk[ic]*x[i][ic]);
          cs[-1][ic][i] = cs[1][ic][i];
          sn[-1][ic][i] = -sn[1][ic][i];
        if (electrode_check(i) == 0) {
          cstr1 += q[i]*cs[1][ic][i];
          sstr1 += q[i]*sn[1][ic][i];
        }
      }
      sfacrl[n] = cstr1;
      sfacim[n++] = sstr1;
    }
  }
  for (m = 2; m <= kmax; m++) {
    for (ic = 0; ic < 3; ic++) {
      sqk = m*unitk[ic] * m*unitk[ic];
      if (sqk <= gsqmx) {
        cstr1 = 0.0;
        sstr1 = 0.0;
        for (i = 0; i < nlocal; i++) {
            cs[m][ic][i] = cs[m-1][ic][i]*cs[1][ic][i] -
              sn[m-1][ic][i]*sn[1][ic][i];
            sn[m][ic][i] = sn[m-1][ic][i]*cs[1][ic][i] +
              cs[m-1][ic][i]*sn[1][ic][i];
            cs[-m][ic][i] = cs[m][ic][i];
            sn[-m][ic][i] = -sn[m][ic][i];
          if (electrode_check(i) == 0) {
            cstr1 += q[i]*cs[m][ic][i];
            sstr1 += q[i]*sn[m][ic][i];
          }
        }
        sfacrl[n] = cstr1;
        sfacim[n++] = sstr1;
      }
    }
  }

  // 1 = (k,l,0), 2 = (k,-l,0)
  for (k = 1; k <= kxmax; k++) {
    for (l = 1; l <= kymax; l++) {
      sqk = (k*unitk[0] * k*unitk[0]) + (l*unitk[1] * l*unitk[1]);
      if (sqk <= gsqmx) {
        cstr1 = 0.0;
        sstr1 = 0.0;
        cstr2 = 0.0;
        sstr2 = 0.0;
        for (i = 0; i < nlocal; i++) {
          if (electrode_check(i) == 0) {
            cstr1 += q[i]*(cs[k][0][i]*cs[l][1][i] - sn[k][0][i]*sn[l][1][i]);
            sstr1 += q[i]*(sn[k][0][i]*cs[l][1][i] + cs[k][0][i]*sn[l][1][i]);
            cstr2 += q[i]*(cs[k][0][i]*cs[l][1][i] + sn[k][0][i]*sn[l][1][i]);
            sstr2 += q[i]*(sn[k][0][i]*cs[l][1][i] - cs[k][0][i]*sn[l][1][i]);
          }
        }
        sfacrl[n] = cstr1;
        sfacim[n++] = sstr1;
        sfacrl[n] = cstr2;
        sfacim[n++] = sstr2;
      }
    }
  }

  // 1 = (0,l,m), 2 = (0,l,-m)

  for (l = 1; l <= kymax; l++) {
    for (m = 1; m <= kzmax; m++) {
      sqk = (l*unitk[1] * l*unitk[1]) + (m*unitk[2] * m*unitk[2]);
      if (sqk <= gsqmx) {
        cstr1 = 0.0;
        sstr1 = 0.0;
        cstr2 = 0.0;
        sstr2 = 0.0;
        for (i = 0; i < nlocal; i++) {
          if (electrode_check(i) == 0) {
            cstr1 += q[i]*(cs[l][1][i]*cs[m][2][i] - sn[l][1][i]*sn[m][2][i]);
            sstr1 += q[i]*(sn[l][1][i]*cs[m][2][i] + cs[l][1][i]*sn[m][2][i]);
            cstr2 += q[i]*(cs[l][1][i]*cs[m][2][i] + sn[l][1][i]*sn[m][2][i]);
            sstr2 += q[i]*(sn[l][1][i]*cs[m][2][i] - cs[l][1][i]*sn[m][2][i]);
          }
        }
        sfacrl[n] = cstr1;
        sfacim[n++] = sstr1;
        sfacrl[n] = cstr2;
        sfacim[n++] = sstr2;
      }
    }
  }

  // 1 = (k,0,m), 2 = (k,0,-m)

  for (k = 1; k <= kxmax; k++) {
    for (m = 1; m <= kzmax; m++) {
      sqk = (k*unitk[0] * k*unitk[0]) + (m*unitk[2] * m*unitk[2]);
      if (sqk <= gsqmx) {
        cstr1 = 0.0;
        sstr1 = 0.0;
        cstr2 = 0.0;
        sstr2 = 0.0;
        for (i = 0; i < nlocal; i++) {
          if (electrode_check(i) == 0) {
            cstr1 += q[i]*(cs[k][0][i]*cs[m][2][i] - sn[k][0][i]*sn[m][2][i]);
            sstr1 += q[i]*(sn[k][0][i]*cs[m][2][i] + cs[k][0][i]*sn[m][2][i]);
            cstr2 += q[i]*(cs[k][0][i]*cs[m][2][i] + sn[k][0][i]*sn[m][2][i]);
            sstr2 += q[i]*(sn[k][0][i]*cs[m][2][i] - cs[k][0][i]*sn[m][2][i]);
          }
        }
        sfacrl[n] = cstr1;
        sfacim[n++] = sstr1;
        sfacrl[n] = cstr2;
        sfacim[n++] = sstr2;
      }
    }
  }

  // 1 = (k,l,m), 2 = (k,-l,m), 3 = (k,l,-m), 4 = (k,-l,-m)

  for (k = 1; k <= kxmax; k++) {
    for (l = 1; l <= kymax; l++) {
      for (m = 1; m <= kzmax; m++) {
        sqk = (k*unitk[0] * k*unitk[0]) + (l*unitk[1] * l*unitk[1]) +
          (m*unitk[2] * m*unitk[2]);
        if (sqk <= gsqmx) {
          cstr1 = 0.0;
          sstr1 = 0.0;
          cstr2 = 0.0;
          sstr2 = 0.0;
          cstr3 = 0.0;
          sstr3 = 0.0;
          cstr4 = 0.0;
          sstr4 = 0.0;
          for (i = 0; i < nlocal; i++) {
            if (electrode_check(i) == 0) {
              clpm = cs[l][1][i]*cs[m][2][i] - sn[l][1][i]*sn[m][2][i];
              slpm = sn[l][1][i]*cs[m][2][i] + cs[l][1][i]*sn[m][2][i];
              cstr1 += q[i]*(cs[k][0][i]*clpm - sn[k][0][i]*slpm);
              sstr1 += q[i]*(sn[k][0][i]*clpm + cs[k][0][i]*slpm);

              clpm = cs[l][1][i]*cs[m][2][i] + sn[l][1][i]*sn[m][2][i];
              slpm = -sn[l][1][i]*cs[m][2][i] + cs[l][1][i]*sn[m][2][i];
              cstr2 += q[i]*(cs[k][0][i]*clpm - sn[k][0][i]*slpm);
              sstr2 += q[i]*(sn[k][0][i]*clpm + cs[k][0][i]*slpm);

              clpm = cs[l][1][i]*cs[m][2][i] + sn[l][1][i]*sn[m][2][i];
              slpm = sn[l][1][i]*cs[m][2][i] - cs[l][1][i]*sn[m][2][i];
              cstr3 += q[i]*(cs[k][0][i]*clpm - sn[k][0][i]*slpm);
              sstr3 += q[i]*(sn[k][0][i]*clpm + cs[k][0][i]*slpm);

              clpm = cs[l][1][i]*cs[m][2][i] - sn[l][1][i]*sn[m][2][i];
              slpm = -sn[l][1][i]*cs[m][2][i] - cs[l][1][i]*sn[m][2][i];
              cstr4 += q[i]*(cs[k][0][i]*clpm - sn[k][0][i]*slpm);
              sstr4 += q[i]*(sn[k][0][i]*clpm + cs[k][0][i]*slpm);
            }
          }
          sfacrl[n] = cstr1;
          sfacim[n++] = sstr1;
          sfacrl[n] = cstr2;
          sfacim[n++] = sstr2;
          sfacrl[n] = cstr3;
          sfacim[n++] = sstr3;
          sfacrl[n] = cstr4;
          sfacim[n++] = sstr4;
        }
      }
    }
  }
  if (runstage == 1) runstage = 2;
}

/* ---------------------------------------------------------------------- */
void FixConp::cg()
{
  int iter,i,j,idx1d;
  double d,beta,ptap,lresnorm,resnorm,netcharge,tmp;
  double res[elenum_all],p[elenum_all],ap[elenum_all];
  for (i = 0; i < elenum_all; i++) eleallq[i] = 0.0;
  lresnorm = 0.0;
  for (i = 0; i < elenum_all; ++i) {
    res[i] = bbb_all[i];
    p[i] = bbb_all[i];
    for (j = 0; j < elenum_all; ++j) {
      idx1d= i*elenum_all+j;
      tmp = aaa_all[idx1d]*eleallq[j];
      res[i] -= tmp;
      p[i] -= tmp;
    }
    lresnorm += res[i]*res[i];
  }
  for (iter = 1; iter < maxiter; ++iter) {
    d = 0.0;
    for (i = 0; i < elenum_all; ++i) {
      ap[i] = 0.0;
      for (j = 0; j < elenum_all; ++j) {
        idx1d= i*elenum_all+j;
        ap[i] += aaa_all[idx1d]*p[j];
      }
    }
    ptap = 0.0;
    for (i = 0; i < elenum_all; ++i) {
      ptap += ap[i]*p[i];
    }
    d = lresnorm/ptap;
    resnorm = 0.0;
    for (i = 0; i <elenum_all; ++i) {
      eleallq[i] = eleallq[i]+d*p[i];
      res[i] = res[i]-d*ap[i];
      resnorm += res[i]*res[i];
    }
    if (resnorm/elenum_all < tolerance) {
      netcharge = 0.0;
      for (i = 0; i < elenum_all; ++i) netcharge += eleallq[i];
      if (me == 0) {
        fprintf(outf,"***** Converged at iteration %d. res = %g netcharge = %g\n",
            iter,resnorm,netcharge);
      }
      break;
    }
    beta = resnorm/lresnorm;
    for (i = 0; i < elenum_all; ++i) {
      p[i] = res[i]+beta*p[i];
    }
    lresnorm = resnorm;
    if (me == 0) {
      fprintf(outf,"Iteration %d: res = %g\n",iter,lresnorm);
    }
  }
}
/* ---------------------------------------------------------------------- */
void FixConp::inv()
{
  int i,j,k,idx1d;
  if (runstage == 2 && a_matrix_f < 2) {
    int m = elenum_all;
    int n = elenum_all;
    int lda = elenum_all;
    int *ipiv = new int[elenum_all+1];
    int lwork = elenum_all*elenum_all;
    double *work = new double[lwork];
    int info;
    int infosum;
    //RS on 22-04-2021: on an earlier date neutrality was enforced below by transforming the inverse of the A-marix to the S-matrix
    /* RS: added code start */
    int ij;
    double AinvE [elenum_all] = {0};
    double EtAinvE = 0 ;
    double EtAinv [elenum_all] = {0};
    /* RS: added code end */

    dgetrf_(&m,&n,aaa_all,&lda,ipiv,&info);
    infosum = info;
    dgetri_(&n,aaa_all,&lda,ipiv,work,&lwork,&info);
    infosum += info;
    delete [] ipiv;
    ipiv = NULL;
    delete [] work;
    work = NULL;

    if (infosum != 0) error->all(FLERR,"Inversion failed!");

      /* RS: added code start */
      for (i = 0; i < elenum_all; i++) {
          for (j = 0; j < elenum_all; j++) {
              ij = i * elenum_all + j;
              AinvE[i] += aaa_all[ij];
          }
      }

      for (i = 0; i < elenum_all; i++) {
          EtAinvE += AinvE[i];
      }

      for (j = 0; j < elenum_all; j++) {
          for (i = 0; i < elenum_all; i++) {
              ij = i * elenum_all + j;
              EtAinv[j] += aaa_all[ij];
          }
      }

      for (i = 0; i < elenum_all; i++) {
          for (j = 0; j < elenum_all; j++) {
              ij = i * elenum_all + j;
              sss_all[ij] = aaa_all[ij] - ((AinvE[i] * EtAinv[j])/EtAinvE);
          }
      }
      /* RS: added code end */

    if (me == 0) {
      FILE *outinva = fopen("inv_a_matrix","w");
      for (i = 0; i < elenum_all; i++) {
        if(i == 0) fprintf (outinva," ");
        fprintf (outinva,"%12d",eleall2tag[i]);
      }
      fprintf (outinva,"\n");
      for (k = 0; k < elenum_all*elenum_all; k++) {
        if (k%elenum_all != 0) {
          fprintf (outinva," ");
        }
        fprintf(outinva,"%20.10f",aaa_all[k]);
        if ((k+1)%elenum_all == 0) {
          fprintf(outinva,"\n");
        }
      }
      fclose(outinva);
    }
  }
  if (runstage == 2) {
      if (method == 1){     //RS: if constant charge
          G_H_cal();
          Q_cal();
          V_cal();
          FILE *out_V_cal = fopen("V_cal", "a");
          fprintf(out_V_cal, "%20s %20s %20s %20s\n", "Step", "Q_R [e]", "V_L [V]", "V_R [V]");
          fprintf(out_V_cal, "%20ld %20.10f %20.10f %20.10f\n", update->ntimestep-1, Q_sum, V_m/evscale, V_p/evscale);
          fclose(out_V_cal);
      }
      runstage = 3;
  }

}
/* ---------------------------------------------------------------------- */
void FixConp::update_charge()
{
  int i,j,idx1d;
  int elealli,tagi;
  double eleallq_i;
  int *tag = atom->tag;
  int nall = atom->nlocal+atom->nghost;
  double *q = atom->q;
  double **x = atom->x;
  //FILE *out_Bq_CPM = fopen("Bq_CPM", "a");

    //FILE *out_q_before_in_update_charge = fopen("q_before_in_update_charge", "a");
    //for (int i = 0; i < nall; i++) {
    //    fprintf(out_q_before_in_update_charge,"%20d %20f\n", tag[i], q[i]);
    //}
    //fclose(out_q_before_in_update_charge);

    //FILE *out_bbb_CP4M2 = fopen("bbb_CP4M2", "a");
    //FILE *out_aaa_all = fopen("aaa_all", "a");
  for (i = 0; i < nall; ++i) {
    if (electrode_check(i)) {
      tagi = tag[i];
      elealli = tag2eleall[tagi];
      if (minimizer == 0) {
        q[i] = eleallq[elealli];
      } else if (minimizer == 1) {
        eleallq_i = 0.0;
        for (j = 0; j < elenum_all; j++) {
          idx1d=elealli*elenum_all+j;
          //eleallq_i += aaa_all[idx1d]*bbb_all[j];
          eleallq_i += sss_all[idx1d]*V_min_Bq_CP4M2[j]; //RS: changed on 27-08-2021 use B by PPPM
          //fprintf(out_Bq_CPM, "%i %i %20.10f %20.10f\n", tagi, elealli, bbb_all[j], bbb_CP4M2[j]);
            //fprintf(out_bbb_CP4M2, "%i %20.10f\n", tagi, bbb_CP4M2[j]);
            //fprintf(out_aaa_all, "%20.10f\n", aaa_all[idx1d]);
        }
        q[i] = eleallq_i;
      }
    }
  }
    //fclose(out_Bq_CPM);
    //fclose(out_bbb_CP4M2);
    //fclose(out_aaa_all);

    //FILE *out_q_after_in_update_charge = fopen("q_after_in_update_charge", "a");
    //for (int i = 0; i < nall; i++) {
    //    fprintf(out_q_after_in_update_charge,"%20d %20f\n", tag[i], q[i]);
    //}
    //fclose(out_q_after_in_update_charge);
}
/* ---------------------------------------------------------------------- */
//void FixConp::force_cal(int vflag)		//RS on 22-04-2021: commented out force_cal function, in which foces due to erfc(eta*rij) and an extra energy term due to eta were calculated
//{
//  int i;
//  if (force->kspace->energy) {
//    double eleqsqsum = 0.0;
//    int nlocal = atom->nlocal;
//    for (i = 0; i < nlocal; i++) {
//      if (electrode_check(i)) {
//        eleqsqsum += atom->q[i]*atom->q[i];
//      }
//    }
//    double tmp;
//    MPI_Allreduce(&eleqsqsum,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
//    eleqsqsum = tmp;
//    double scale = 1.0;
//    double qscale = force->qqrd2e*scale;
//    force->kspace->energy += qscale*eta*eleqsqsum/(sqrt(2)*MY_PIS);
//  }
//  coul_cal(0,NULL,NULL);

//}
/* ---------------------------------------------------------------------- */
void FixConp::coul_cal(int coulcalflag,double *m,int *ele2tag)
{
  Ctime1 = MPI_Wtime();
  //coulcalflag = 2: a_cal; 1: b_cal; 0: force_cal
  int i,j,k,ii,jj,jnum,itype,jtype,idx1d;
  int checksum,elei,elej,elealli,eleallj;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz;
  double r,r2inv,rsq,grij,etarij,expm2,t,erfc,dudq;
  double forcecoul,ecoul,prefactor,fpair;

  int inum = force->pair->list->inum;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  int *atomtype = atom->type;
  int *tag = atom->tag;
  int *ilist = force->pair->list->ilist;
  int *jlist;
  int *numneigh = force->pair->list->numneigh;
  int **firstneigh = force->pair->list->firstneigh;

  double qqrd2e = force->qqrd2e;
  double **cutsq = force->pair->cutsq;
  int itmp;
  double *p_cut_coul = (double *) force->pair->extract("cut_coul",itmp);
  double cut_coulsq = (*p_cut_coul)*(*p_cut_coul);
  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = atomtype[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      checksum = abs(electrode_check(i))+abs(electrode_check(j));
      if (checksum == 1 || checksum == 2) {
        if (coulcalflag == 0 || checksum == coulcalflag) {
          delx = xtmp - x[j][0];
          dely = ytmp - x[j][1];
          delz = ztmp - x[j][2];
          rsq = delx*delx + dely*dely + delz*delz;
          jtype = atomtype[j];
          if (rsq < cutsq[itype][jtype]) {
            r2inv = 1.0/rsq;
            if (rsq < cut_coulsq) {
              dudq =0.0;
              r = sqrt(rsq);
              if (coulcalflag != 0) {
                grij = g_ewald * r;
                expm2 = exp(-grij*grij);
                t = 1.0 / (1.0 + EWALD_P*grij);
                erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
                dudq = erfc/r;
              }
              //if (checksum == 1) etarij = eta*r;			//RS on 22-04-2021: this line, the next 4 lines and a line further in the code (dudq = -erfc/r) are commented out to remove erfc(eta*rij) terms from A and B matrix)
              //else if (checksum == 2) etarij = eta*r/sqrt(2);
              //expm2 = exp(-etarij*etarij);
              //t = 1.0 / (1.0+EWALD_P*etarij);
              //erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;

              //if (coulcalflag == 0) {					//RS on 22-04-2021: in force_cal, coulcalflag =0, since force_cal is removed also this part is removed
               // prefactor = qqrd2e*qtmp*q[j]/r;
                //forcecoul = -prefactor*(erfc+EWALD_F*etarij*expm2);
                //fpair = forcecoul*r2inv;
                //f[i][0] += delx*forcecoul;
                //f[i][1] += dely*forcecoul;
                //f[i][2] += delz*forcecoul;
                //if (newton_pair || j < nlocal) {
                  //f[j][0] -= delx*forcecoul;
                  //f[j][1] -= dely*forcecoul;
                  //f[j][2] -= delz*forcecoul;
                //}
                //ecoul = -prefactor*erfc;
                //force->pair->ev_tally(i,j,nlocal,newton_pair,0,ecoul,fpair,delx,dely,delz); //evdwl=0
              //} else {
	      if (coulcalflag != 0) {	//RS on 22-04-2021: added this line since "if (coulcalflag == 0) {" was removed earlier
                //dudq -= erfc/r;	//RS on 22-04-2021: commented out to eliminate erfc(eta*rij) terms in A and B matrix
                elei = -1;
                elej = -1;
                for (k = 0; k < elenum; ++k) {
                  if (i < nlocal) {
                    if (ele2tag[k] == tag[i]) {
                      elei = k;
                      if (coulcalflag == 1) {
                        m[k] -= q[j]*dudq;
                        break;
                      }
                    }
                  }
                  if (j < nlocal) {
                    if (ele2tag[k] == tag[j]) {
                      elej = k;
                      if (coulcalflag == 1) {
                        m[k] -= q[i]*dudq;
                        break;
                      }
                    }
                  }
                }
                if (coulcalflag == 2 && checksum == 2) {
                  elealli = tag2eleall[tag[i]];
                  eleallj = tag2eleall[tag[j]];
                  if (elei != -1) {
                    idx1d = elei*elenum_all+eleallj;
                    m[idx1d] += dudq;
                  }
                  if (elej != -1) {
                    idx1d = elej*elenum_all+elealli;
                    m[idx1d] += dudq;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  Ctime2 = MPI_Wtime();
  Ctime += Ctime2-Ctime1;
}

/* ---------------------------------------------------------------------- */
double FixConp::rms(int km, double prd, bigint natoms, double q2)
{
  double value = 2.0*q2*g_ewald/prd *
    sqrt(1.0/(MY_PI*km*natoms)) *
    exp(-MY_PI*MY_PI*km*km/(g_ewald*g_ewald*prd*prd));
  return value;
}

/* ---------------------------------------------------------------------- */
void FixConp::coeffs()
{
  int k,l,m;
  double sqk;

  double g_ewald_sq_inv = 1.0 / (g_ewald*g_ewald);
  double preu = 4.0*MY_PI/volume;

  kcount = 0;

  // (k,0,0), (0,l,0), (0,0,m)

  for (m = 1; m <= kmax; m++) {
    sqk = (m*unitk[0]) * (m*unitk[0]);
    if (sqk <= gsqmx) {
      kxvecs[kcount] = m;
      kyvecs[kcount] = 0;
      kzvecs[kcount] = 0;
      ug[kcount] = preu*exp(-0.25*sqk*g_ewald_sq_inv)/sqk;
      kcount++;
    }
    sqk = (m*unitk[1]) * (m*unitk[1]);
    if (sqk <= gsqmx) {
      kxvecs[kcount] = 0;
      kyvecs[kcount] = m;
      kzvecs[kcount] = 0;
      ug[kcount] = preu*exp(-0.25*sqk*g_ewald_sq_inv)/sqk;
      kcount++;
    }
    sqk = (m*unitk[2]) * (m*unitk[2]);
    if (sqk <= gsqmx) {
      kxvecs[kcount] = 0;
      kyvecs[kcount] = 0;
      kzvecs[kcount] = m;
      ug[kcount] = preu*exp(-0.25*sqk*g_ewald_sq_inv)/sqk;
      kcount++;
    }
  }

  // 1 = (k,l,0), 2 = (k,-l,0)

  for (k = 1; k <= kxmax; k++) {
    for (l = 1; l <= kymax; l++) {
      sqk = (unitk[0]*k) * (unitk[0]*k) + (unitk[1]*l) * (unitk[1]*l);
      if (sqk <= gsqmx) {
        kxvecs[kcount] = k;
        kyvecs[kcount] = l;
        kzvecs[kcount] = 0;
        ug[kcount] = preu*exp(-0.25*sqk*g_ewald_sq_inv)/sqk;
        kcount++;

        kxvecs[kcount] = k;
        kyvecs[kcount] = -l;
        kzvecs[kcount] = 0;
        ug[kcount] = preu*exp(-0.25*sqk*g_ewald_sq_inv)/sqk;
        kcount++;;
      }
    }
  }

  // 1 = (0,l,m), 2 = (0,l,-m)

  for (l = 1; l <= kymax; l++) {
    for (m = 1; m <= kzmax; m++) {
      sqk = (unitk[1]*l) * (unitk[1]*l) + (unitk[2]*m) * (unitk[2]*m);
      if (sqk <= gsqmx) {
        kxvecs[kcount] = 0;
        kyvecs[kcount] = l;
        kzvecs[kcount] = m;
        ug[kcount] = preu*exp(-0.25*sqk*g_ewald_sq_inv)/sqk;
        kcount++;

        kxvecs[kcount] = 0;
        kyvecs[kcount] = l;
        kzvecs[kcount] = -m;
        ug[kcount] = preu*exp(-0.25*sqk*g_ewald_sq_inv)/sqk;
        kcount++;
      }
    }
  }

  // 1 = (k,0,m), 2 = (k,0,-m)

  for (k = 1; k <= kxmax; k++) {
    for (m = 1; m <= kzmax; m++) {
      sqk = (unitk[0]*k) * (unitk[0]*k) + (unitk[2]*m) * (unitk[2]*m);
      if (sqk <= gsqmx) {
        kxvecs[kcount] = k;
        kyvecs[kcount] = 0;
        kzvecs[kcount] = m;
        ug[kcount] = preu*exp(-0.25*sqk*g_ewald_sq_inv)/sqk;
        kcount++;

        kxvecs[kcount] = k;
        kyvecs[kcount] = 0;
        kzvecs[kcount] = -m;
        ug[kcount] = preu*exp(-0.25*sqk*g_ewald_sq_inv)/sqk;
        kcount++;
      }
    }
  }

  // 1 = (k,l,m), 2 = (k,-l,m), 3 = (k,l,-m), 4 = (k,-l,-m)

  for (k = 1; k <= kxmax; k++) {
    for (l = 1; l <= kymax; l++) {
      for (m = 1; m <= kzmax; m++) {
        sqk = (unitk[0]*k) * (unitk[0]*k) + (unitk[1]*l) * (unitk[1]*l) +
          (unitk[2]*m) * (unitk[2]*m);
        if (sqk <= gsqmx) {
          kxvecs[kcount] = k;
          kyvecs[kcount] = l;
          kzvecs[kcount] = m;
          ug[kcount] = preu*exp(-0.25*sqk*g_ewald_sq_inv)/sqk;
          kcount++;

          kxvecs[kcount] = k;
          kyvecs[kcount] = -l;
          kzvecs[kcount] = m;
          ug[kcount] = preu*exp(-0.25*sqk*g_ewald_sq_inv)/sqk;
          kcount++;

          kxvecs[kcount] = k;
          kyvecs[kcount] = l;
          kzvecs[kcount] = -m;
          ug[kcount] = preu*exp(-0.25*sqk*g_ewald_sq_inv)/sqk;
          kcount++;

          kxvecs[kcount] = k;
          kyvecs[kcount] = -l;
          kzvecs[kcount] = -m;
          ug[kcount] = preu*exp(-0.25*sqk*g_ewald_sq_inv)/sqk;
          kcount++;
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixConp::pot_wall_wall()
{
    //if (me == 0) utils::logmesg(lmp,"pot_wall_wall() start ...\n");

    int* tag = atom->tag;
    int nlocal = atom -> nlocal;
    double* q = atom->q;
    double** f = atom->f;
    int nall = atom->nlocal + atom->nghost;

    //set settings of the kspace object declared in the CPM code to those of the the kspace object used by Lammps
    obj_kspace.accuracy_relative = force->kspace->accuracy_relative;
    obj_kspace.slabflag = force->kspace->slabflag;
    obj_kspace.slab_volfactor = force->kspace->slab_volfactor;

    obj_kspace.init();
    obj_kspace.setup();

    //determine the size of q with the charges, this is larger than the number of particles and the charges are repeating
    // THE SIZE OF Q IS ACTUALLY NALL

    //int size_q = 0;
    //while (tag[size_q]!=0){
    //    size_q++;
    //}

    //FILE *out_size_q_and_nall_0 = fopen("size_q_and_n_all_0", "a");
    //fprintf(out_size_q_and_nall_0,"%20s %20d\n", "size_q", size_q);
    //fprintf(out_size_q_and_nall_0,"%20s %20d\n\n", "nall", nall);
    //fclose(out_size_q_and_nall_0);

    //FILE *out_charges_0 = fopen("out_charges_0", "a");
    //fprintf(out_charges_0,"%20s %20s\n", "tag[i]", "q[i]");
    //for (int i = 0; i < size_q+1000; i++) {
    //    fprintf(out_charges_0,"%20d %20f\n", tag[i], q[i]);
    //}
    //fprintf(out_charges_0,"\n");
    //fclose(out_charges_0);

    //make array with charges used for the calculation of the potential on the wall particles due to the wall particles
    //note the charges of the wall particles are set to the constant q_L and q_R, which are not the real charges in the simulation
    //doing so is convenient because we are interested in the potentials on the wall particles due to the ions, which are
    //calculated by subtracting the potentials on the wall particles due to the wall particles from the potentials on the
    //wall particles due to both the ions and wall particles, i.e., the charges of the wall particles are not relevant
    //this approach allows us to calculate the potentials on the wall particles due to the wall particles only once
    double q_for_cal[nall];
    for (int i = 0; i < nlocal; i++) {
        if (electrode_check(i) == 1) {
            q_for_cal[tag[i]] = q_L;            //the elements in the array are sorted by the global index of the atoms
        } else if (electrode_check(i) == -1){
            q_for_cal[tag[i]] = q_R;
        } else {
            q_for_cal[tag[i]] = 0;
        }
    }

    //make array that has the same size of the array q
    double q_for_cal_2[nall];
    for (int i = 0; i < nall; i++) {
        q_for_cal_2[i] = q_for_cal[tag[i]];     //the elements are sorted by the local index of the processor

    }

    atom->q = &q_for_cal_2[0];  //set the charges to the values required for the calculation

    //Make pointer to a pointer for the forces (atom->f) point to addresses with zero values.
    //At the end of the function, atom->f is reset to point at the "real" addresses again.
    //This is done because the compute() function changes the forces, which should not happen since we only want to
    //calculate the potentials.
    double **f_for_cal;
    f_for_cal = new double*[nlocal];

    for (int i = 0; i < nlocal; i++) {
        f_for_cal[i] = new double[3] ;
    }

    for (int i = 0; i < nlocal; i++) {
        f_for_cal[i][0] = 0;
        f_for_cal[i][1] = 0;
        f_for_cal[i][2] = 0;
    }

    atom->f = f_for_cal;


    obj_kspace.compute(3,1);    //calculates the kspace potentials, used eflag = 3 and vflag =1, as these were the values for the simulations for which we calculated per atom potentials

    //settings for the calculation of the short range potentials
    char *arg_coeff[2] = { "*", "*"};       //to use "pair_coeff * *"

    //extract the global Coulombic cutoff value specified in the pair_style command in the Lammps input, transform it
    //into a string and pass it in an array of strings, where the number of arrays is 1. This is done to make it
    //compatible with the input type of obj_CoulLong.settings
    int itmp;
    double cut_coul = *(double *) force->pair->extract("cut_coul",itmp);
    char string[16];
    sprintf(string,"%f",cut_coul);
    char *arg_settings[1] = {string};

    obj_CoulLong.coeff(2,arg_coeff);                //sets pair_coeff arguments
    obj_CoulLong.settings(1,arg_settings);          //sets second argument to global Coulombic cutoff
    obj_CoulLong.list = force->pair->list;          //since obj_CoulLong.list was not pointing to the same address as force->pair->list

    obj_CoulLong.init();                            //init() sets some internal flags based on user settings

    obj_CoulLong.setup();                           //setup() is called before each dynamic run. Compute forces and counters are also initialised
    obj_CoulLong.compute(3,1);                      //To calculate the short range potentials

    //Adding up the kspace and short range potentials on the wall particles due to the wall particles and store them in
    // the array Psi_w_w.
    int cntr = 0;
    int tag_Psi_w_w[elenum_all];
    for (int i = 0; i < nlocal; i++) {
            if ( (electrode_check(i) == 1) || (electrode_check(i) == -1)) {
                Psi_w_w[cntr] = *(obj_kspace.eatom +i) + *(obj_CoulLong.eatom +i);
                tag_Psi_w_w[cntr] = tag[i];
                cntr++;
            }
    }

    //printing Psi_w_w, also check whether tag_Psi_w_w and eleall2tag match
    //FILE *out_Psi_w_w = fopen("Psi_w_w", "a");
    //for (int i = 0; i < elenum_all; i++) {
    //    fprintf(out_Psi_w_w,"%20d %20d %20d %20f\n", i, tag_Psi_w_w[i], eleall2tag[i], Psi_w_w[i]);
    //}
    //fclose(out_Psi_w_w);

    //reset the values set for the calculation
    atom->q = q;
    atom->f = f;

    for (int i = 0; i < nlocal; i++) {
        delete [] f_for_cal[i];
    }

    delete [] f_for_cal;

   // if (me == 0) utils::logmesg(lmp,"pot_wall_wall() end...\n");

}

void FixConp::pot_wall_ions()
{
    //This function calculates the potential on the wall particles due to both the ions and wall particles
    //At the end of the function the potential on the wall particles due to the ions are calculated by subtracting the
    //potentials on the wall particles due to the wall particles from the potentials on the wall particles due to both
    //the ions and wall particles.

    //if (me == 0) utils::logmesg(lmp,"pot_wall_ions() start...\n");

    if (runstage == 1) runstage = 2;

    int *tag = atom->tag;
    int nlocal = atom -> nlocal;
    double *q = atom->q;
    double **f = atom->f;
    int nall = atom->nlocal + atom->nghost;

    //determine the size of q
    //int size_q = 0;
    //while (tag[size_q]!=0){
    //    size_q++;
    //}

    //FILE *out_size_q_and_nall = fopen("size_q_and_n_all", "a");
    //fprintf(out_size_q_and_nall,"%20s %20d\n", "size_q", size_q);
    //fprintf(out_size_q_and_nall,"%20s %20d\n\n", "nall", nall);
    //fclose(out_size_q_and_nall);

    //FILE *out_charges = fopen("out_charges", "a");
    //fprintf(out_charges,"%20s %20s\n", "tag[i]", "q[i]");
    //for (int i = 0; i < size_q+1000; i++) {
    //    fprintf(out_charges,"%20d %20f\n", tag[i], q[i]);
    //}
    //fprintf(out_charges,"\n");
    //fclose(out_charges);

    double q_for_cal[nall];
    for (int i = 0; i < nlocal; i++) {
        if (electrode_check(i) == 1) {
            q_for_cal[tag[i]] = q_L;
        }
        else if (electrode_check(i) == -1){
            q_for_cal[tag[i]] = q_R;
        }
        else{
            q_for_cal[tag[i]] = q[i];       //the charges of the ions are set to their real values
        }
    }

    double q_for_cal_2[nall];
    for (int i = 0; i < nall; i++) {
        q_for_cal_2[i] = q_for_cal[tag[i]];

    }

    atom->q = &q_for_cal_2[0];

   /* double** f_for_cal;
    f_for_cal = new double[nlocal][3];
    for (int i = 0; i < nlocal; i++) {
        f_for_cal[i][0] = 0;
        f_for_cal[i][1] = 0;
        f_for_cal[i][2] = 0;
    }

    //double* p_f_for_cal = &f_for_cal[0][0];

    //atom->f = &p_f_for_cal;

    atom->f = f_for_cal;       //error: cannot convert ‘double (*)[nlocal][3]’ to ‘double**’ in assignment

    //atom->f = &f_for_cal[0];    //error: cannot convert ‘double (*)[nlocal][3]’ to ‘double**’ in assignment


    //for(int i = 0; i < nlocal; i++){  //gives Segmentation fault (core dumped)
    //    atom->f[i] = f_for_cal[i];
    //}
   */

    double **f_for_cal;
    f_for_cal = new double*[nlocal];

    for (int i = 0; i < nlocal; i++) {
        f_for_cal[i] = new double[3] ;
    }

    for (int i = 0; i < nlocal; i++) {
            f_for_cal[i][0] = 0;
            f_for_cal[i][1] = 0;
            f_for_cal[i][2] = 0;
    }


    atom->f = f_for_cal;        //makes atom->f point to dummy array

    obj_kspace.compute(3,1);    //used eflag = 3 and vflag =1, as these were the values for the simulations for which we calculated per atom potentials
    //information from Lammps mailing list:     eflag != 0 means: compute energy contributions in this step
    //                                          vflag != 0 means: compute virial contributions in this step //
    //                                          the exact value indicates whether per atom or total contributions are supposed to be computed.

    obj_CoulLong.compute(3,1);

    /// \todo This a hack to be able to write restart file at end of the run, since nrequest became 2
    for (int i = 0; i < neighbor->nrequest; i++) {
        delete neighbor->requests[i];
        neighbor->requests[i] = nullptr;
    }
    neighbor->nrequest = 0;

    int cntr = 0;
    int tag_Psi_w_wi[elenum_all];
    int i_saved[elenum];            //save i, so the function electrode_check() can be used to determine if it concerns the left or right wall
    for (int i = 0; i < nlocal; i++) {
        if ( (electrode_check(i) == 1) || (electrode_check(i) == -1)) {
            Psi_w_wi[cntr] = *(obj_kspace.eatom +i) + *(obj_CoulLong.eatom +i);
            tag_Psi_w_wi[cntr] = tag[i];
            i_saved[cntr] = i;
            cntr++;

        }
    }

    //FILE *out_Psi_w_wi = fopen("Psi_w_wi", "a");
    //for (int i = 0; i < elenum_all; i++) {
    //    fprintf (out_Psi_w_wi,"%20d %20d %20d %20f\n", i, tag_Psi_w_wi[i], eleall2tag[i], Psi_w_wi[i]);
    //}
    //fprintf(out_Psi_w_wi,"\n");
    //fclose(out_Psi_w_wi);


    for (int i = 0; i < elenum_all; i++) {
        Psi_w_i[i] = Psi_w_wi[i] - Psi_w_w[i];
        if (electrode_check(i_saved[i]) == 1){
            Bq_CP4M2[i] = (evscale*conv*Psi_w_i[i])/(0.5*q_L);
            if (method == 0){
                V_min_Bq_CP4M2[i] = vL - Bq_CP4M2[i];   //potential on the left wall particles, fo constant potential
            }
        }
        else {
            Bq_CP4M2[i] = (evscale*conv*Psi_w_i[i])/(0.5*q_R);
            if (method == 0){
                V_min_Bq_CP4M2[i] = vR - Bq_CP4M2[i];   //potential on the right wall particles, for constant potential
            }
        }
    }

    /*

   FILE *out_Psi_w_i = fopen("Psi_w_i", "a");
    for (int i = 0; i < elenum; i++) {
        fprintf (out_Psi_w_i,"%20d %20.10f %20.10f\n", eleall2tag[i], Psi_w_i[i], Bq_CP4M2[i]);
    }
    fprintf(out_Psi_w_i,"\n");
    fclose(out_Psi_w_i);

*/




    //reset the values set for the calculation
    atom->q = q;

    atom->f = f;    //make atom->f point at the real force array

    for (int i = 0; i < nlocal; i++) {
            delete [] f_for_cal[i];
    }

    delete [] f_for_cal;

    //for (int i = 0; i < nlocal; i++) {
    //    f[i][0] = 0;
    //    f[i][1] = 0;
     //   f[i][2] = 0;
    //}

    //if (me == 0) utils::logmesg(lmp,"pot_wall_ions() end...\n");

}

void FixConp::G_H_cal() {
    int nlocal = atom->nlocal;
    double *q = atom->q;
    int Eplus[elenum_all] = {0};
    int Emin[elenum_all] = {0};
    //double Q_sum = 0.0;

    double AinvEplus[elenum_all]= {0.0};
    double AinvEmin[elenum_all] = {0.0};
    double D_pp = 0.0;
    double D_mm = 0.0;
    double D_pm = 0.0;
    double D_mp = 0.0;
    double DEdivDD_1[elenum_all] = {0.0};
    double DEdivDD_2[elenum_all] = {0.0};

    /*
    double O_1[elenum_all*elenum_all] = {0.0};
    double O_2[elenum_all*elenum_all] = {0.0};
    double OAinv_1[elenum_all*elenum_all] = {0.0};
    double OAinv_2[elenum_all*elenum_all]; //= {0.0};
    double Z[elenum_all*elenum_all] = {0.0};
     */

    int ij = 0;
    int ik = 0;
    int kj = 0;


    //FILE *out_Q_sum=fopen("Q_sum", "a");
    //fprintf(out_Q_sum, "%20.10d", nlocal);
    int j = 0;
    for (int i = 0; i < nlocal; i++) {
        if (electrode_check(j) == 1) {
            Emin[j] = 1;
            j++;
        }
        if (electrode_check(i) == -1) {
            Eplus[j] = 1;
            j++;
            //Q_sum += q[i];                  //Always apply positive potential on second (right) electrode
        }
    }
    //fprintf(out_Q_sum, "%20.10f", Q_sum);
    //fclose(out_Q_sum);


    for (int i = 0; i < elenum_all; i++) {
        for (int j = 0; j < elenum_all; j++) {
            ij = i * elenum_all + j;
            AinvEplus[i] += aaa_all[ij]*Eplus[j];
            AinvEmin[i] += aaa_all[ij]*Emin[j];
        }
    }

    for (int i = 0; i < elenum_all; i++) {
        D_pp += Eplus[i]*AinvEplus[i];
        D_pm += Eplus[i]*AinvEmin[i];
        D_mp += Emin[i]*AinvEplus[i];
        D_mm += Emin[i]*AinvEmin[i];
    }

    for (int i = 0; i < elenum_all; i++) {
        DEdivDD_1[i] = (D_mm*Eplus[i] - D_mp*Emin[i])/(D_pp*D_mm - D_pm*D_mp);
        DEdivDD_2[i] = (-D_pm*Eplus[i] + D_pp*Emin[i])/(D_pp*D_mm - D_pm*D_mp);
    }


    for (int i = 0; i < elenum_all; i++) {
        //GQ[i] = 0.0;                            //initialise to 0
        G[i] = 0.0;
        for (int j = 0; j < elenum_all; j++) {
            ij = i * elenum_all + j;
            //GQ[i] += aaa_all[ij]*(DEdivDD_1[j] - DEdivDD_2[j])*Q_sum;
            G[i] += aaa_all[ij]*(DEdivDD_1[j] - DEdivDD_2[j]);
        }
    }

    for (int i = 0; i < elenum_all; i++){
        for (int j = 0; j < elenum_all; j++){
            ij = i * elenum_all + j;
            O_1[ij] = DEdivDD_1[i]*Eplus[j];
            O_2[ij] = DEdivDD_2[i]*Emin[j];
        }
    }

    for (int i = 0; i < elenum_all; i++){
        for (int j = 0; j < elenum_all; j++){
            ij = i*elenum_all + j;
            OAinv_1[ij] = 0.0;
            OAinv_2[ij] = 0.0;
            for (int k = 0; k < elenum_all; k++){
                ik = i*elenum_all + k;
                kj = k*elenum_all + j;
                OAinv_1[ij] += O_1[ik]*aaa_all[kj];
                OAinv_2[ij] += O_2[ik]*aaa_all[kj];
            }
        }
    }


    for (int i = 0; i < elenum_all; i++){
        for (int j = 0; j < elenum_all; j++){
            ij = i*elenum_all + j;
            Z[ij] = -OAinv_1[ij] - OAinv_2[ij];
            if (i == j) {
                Z[ij] += 1;
            }
        }
    }


  for (int i = 0; i < elenum_all; i++){
      for (int j = 0; j < elenum_all; j++){
          ij = i*elenum_all + j;
          H[ij] = 0.0;
          for (int k = 0; k < elenum_all; k++){
              ik = i*elenum_all + k;
              kj = k*elenum_all + j;
              H[ij] += -aaa_all[ik]*Z[kj];
          }
      }
  }

  double one_over_DD = 1.0/(D_pp*D_mm - D_pm*D_mp);
  C_pp = one_over_DD*D_mm;
  C_pm = one_over_DD*(-D_pm);
  C_mp = one_over_DD*(-D_mp);
  C_mm = one_over_DD*(D_pp);


  /*
    for (int i = 0; i < elenum_all; i++) {
        AQ_1[i] = 0.0;                            //initialise to 0
        for (int j = 0; j < elenum_all; j++) {
            ij = i * elenum_all + j;
            AQ_1[i] += (DEdivDD_1[j] - DEdivDD_2[j])*Q_sum;       //first (constant) term of AQ, the second term is calculated in update_charge_2
        }
    }


      FILE *out_Ainv=fopen("Ainv", "a");
      for (int i = 0; i < elenum_all; i++) {
          for (int j = 0; j < elenum_all; j++){
              ij = i * elenum_all + j;
              fprintf(out_Ainv, "%20.10f", aaa_all[ij]);
          }
          fprintf(out_Ainv, "\n");
      }
      fclose(out_Ainv);

      FILE *out_AinvEpm =fopen("AinvEpm", "a");
      //fprintf(out_AinvEpm, "%20s %20s\n", "AinvEmin", "AinvEplus");
      for (int i = 0; i < elenum_all; i++) {
          fprintf(out_AinvEpm, "%20.10f %20.10f\n", AinvEmin[i], AinvEplus[i]);
      }
      fclose(out_AinvEpm);


      FILE *out_Emin_Eplus = fopen("Emin_Eplus", "a");
      fprintf(out_Emin_Eplus, "%20s %20d\n", "j=", j);
      fprintf(out_Emin_Eplus, "%20s %20f\n", "Q=", Q_sum);
      fprintf(out_Emin_Eplus, "%20s %20f\n", "Dpp=", D_pp);
      fprintf(out_Emin_Eplus, "%20s %20f\n", "Dpm=", D_pm);
      fprintf(out_Emin_Eplus, "%20s %20f\n", "Dmp=", D_mp);
      fprintf(out_Emin_Eplus, "%20s %20f\n", "Dmm=", D_mm);
      for (int i = 0; i < elenum_all; i++) {
          fprintf(out_Emin_Eplus, "%20d %20d %20.10f %20.10f\n", Emin[i], Eplus[i], DEdivDD_1[i], DEdivDD_2[i]);
      }
      fclose(out_Emin_Eplus);


   FILE *out_GQ=fopen("GQ", "a");
   for (int i = 0; i < elenum_all; i++) {
       fprintf(out_GQ, "%20.10f\n", GQ[i]);
   }
   fclose(out_GQ);


    FILE *out_OAinv_1=fopen("OAinv_1", "a");
    for (int i = 0; i < elenum_all; i++) {
        for (int j = 0; j <elenum_all; j++){
            ij = i * elenum_all + j;
            fprintf(out_OAinv_1, "%20.10f", OAinv_1[ij]);
        }
        fprintf(out_OAinv_1, "\n");
    }
    fclose(out_OAinv_1);


    FILE *out_OAinv_2=fopen("O_2", "a");
    for (int i = 0; i < elenum_all; i++) {
        for (int j = 0; j <elenum_all; j++){
            ij = i * elenum_all + j;
            fprintf(out_O_2, "%20.10f", O_2[ij]);
        }
        fprintf(out_O_2, "\n");
    }
    fclose(out_O_2);


    FILE *out_OAinv_2=fopen("OAinv_2", "a");
    for (int i = 0; i < elenum_all; i++) {
        for (int j = 0; j <elenum_all; j++){
            ij = i * elenum_all + j;
            fprintf(out_OAinv_2, "%20.10f", OAinv_2[ij]);
        }
        fprintf(out_OAinv_2, "\n");
    }
    fclose(out_OAinv_2);

    */

    FILE *out_Z=fopen("Z", "a");
    for (int i = 0; i < elenum_all; i++) {
        for (int j = 0; j <elenum_all; j++){
            ij = i * elenum_all + j;
            fprintf(out_Z, "%20.10f", Z[ij]);
        }
        fprintf(out_Z, "\n");
    }
    fclose(out_Z);

    FILE *out_H=fopen("H", "a");
    for (int i = 0; i < elenum_all; i++) {
        for (int j = 0; j <elenum_all; j++){
            ij = i * elenum_all + j;
            fprintf(out_H, "%20.10f", H[ij]);
        }
        fprintf(out_H, "\n");
    }
    fclose(out_H);


}

void FixConp::update_charge_2()
{
    int i,j,idx1d;
    int elealli,tagi;
    double eleallq_i;
    int *tag = atom->tag;
    int nall = atom->nlocal+atom->nghost;
    double *q = atom->q;
    double **x = atom->x;

    Q_cal();
    V_cal();                    //calculates potentials V_m and V_p

    bigint curr_step = update->ntimestep - 1;
    int print_every = 100;
    int mod = curr_step % print_every;
    bigint first_step = update->firststep;

    if (mod == 0) {
        FILE *out_V_cal = fopen("V_cal", "a");
        fprintf(out_V_cal, "%20ld %20.10f %20.10f %20.10f\n", curr_step, Q_sum, V_m/evscale, V_p/evscale);
        fclose(out_V_cal);
    }
    if (R != 0.0){
        double dt = update->dt;     //timestep, have to use 1 for everynum in constant potential command
        double I = (Vs-(V_p - V_m))/R;
        double dQ = I*dt;
        Q_sum += dQ;
    }

    //FILE *out_Q_cal = fopen("Q_cal", "a");
    for (i = 0; i < nall; ++i) {
        if (electrode_check(i)) {
            tagi = tag[i];
            elealli = tag2eleall[tagi];
            if (minimizer == 0) {
                q[i] = eleallq[elealli];
            } else if (minimizer == 1) {
                eleallq_i = G[elealli]*Q_sum;
                for (j = 0; j < elenum_all; j++) {
                    idx1d=elealli*elenum_all+j;
                    eleallq_i += H[idx1d]*Bq_CP4M2[j];
                    //if (electrode_check(i) == -1){
                    //    Ep_dot_Ainv_Bq += aaa_all[idx1d]*Bq_CP4M2[j];
                    //}
                    //else {
                    //    Em_dot_Ainv_Bq += aaa_all[idx1d]*Bq_CP4M2[j];
                    //}
                }
                //fprintf(out_Q_cal, "%20.10d %20.10f\n", tagi, eleallq_i);
                q[i] = eleallq_i;
            }
            //if (electrode_check(i) == -1){
            //    Q_sum += q[i];
            //}

        }
    }

    //Q_cal();
    //V_cal();
    //fclose(out_Q_cal);

    //V_p = C_pp*(Q_sum + Ep_dot_Ainv_Bq) + C_pm*(-Q_sum + Em_dot_Ainv_Bq);
    //V_m = C_mp*(Q_sum + Ep_dot_Ainv_Bq) + C_mm*(-Q_sum + Em_dot_Ainv_Bq);

    //FILE *out_V_cal = fopen("V_cal", "a");
    //fprintf(out_V_cal, "%20.10f %20.10f\n", V_m, V_p);
    //fclose(out_V_cal);


}

//The below function calculates the forces on the ions if the wall charges are set to zero
void FixConp::force_extra_ions_1()
{
    int *tag = atom->tag;
    int nlocal = atom -> nlocal;
    double *q = atom->q;
    double **f = atom->f;
    int nall = atom->nlocal + atom->nghost;

    double q_for_cal[nall] = {0.0};    //dummy array for the charges of the particles
    for (int i = 0; i < nall; i++) {
        if (electrode_check(i)) {
            q_for_cal[i] = 0.0;        //charges of the wall particles are set to 0
        }
        else{
            q_for_cal[i] = q[i];       //the charges of the ions are set to their real values
        }
    }

    atom->q = &q_for_cal[0];           //q now points the values of q_for_cal

    //declaration and intitialisation of dummy array for forces
    double **f_extra;
    f_extra = new double*[nlocal];      //we are using "newton off” as instructed by Wang et al., forces are then tallied to the owned atoms.
    for (int i = 0; i < nlocal; i++) {
        f_extra[i] = new double[3] ;
    }

    for (int i = 0; i < nlocal; i++) {
        f_extra[i][0] = 0;
        f_extra[i][1] = 0;
        f_extra[i][2] = 0;
    }


    atom->f = f_extra;        //makes atom->f point to the dummy array

    //check force before computation
    FILE *out_F_before = fopen("F_before", "a");
    for (int i = 0; i < nlocal; i++) {
       fprintf (out_F_before,"%20d %20.10f %20.10f %20.10f\n", tag[i], f_extra[i][0], f_extra[i][1], f_extra[i][2]);
    }
    fprintf(out_F_before,"\n");
    fclose(out_F_before);


    obj_kspace.compute(3,1);    //calculation of long range electrostatic interaction

    obj_CoulLong.compute(3,1);  //calculation of short range electrostatic interaction

    //print forces
    FILE *out_F_after_0 = fopen("F_after_0", "a");
    for (int i = 0; i < nlocal; i++) {
        fprintf (out_F_after_0,"%20d %20.10f %20.10f %20.10f\n", tag[i], f_extra[i][0], f_extra[i][1], f_extra[i][2]);
    }
    fprintf(out_F_after_0,"\n");
    fclose(out_F_after_0);


    /// \todo This a hack to be able to write restart file at end of the run, since nrequest became 2
    neighbor->nrequest = 0;

    //subtracting f_extra from forces on ions
    for (int i = 0; i < nlocal; i++) {
        if (electrode_check(i)==0){
            f[i][0] = f[i][0] - f_extra[i][0];
            f[i][1] = f[i][1] - f_extra[i][1];
            f[i][2] = f[i][2] - f_extra[i][2];
        }
    }

    /*

    FILE *out_F_after_0 = fopen("F_after_0", "a");
    for (int i = 0; i < nlocal; i++) {
        fprintf (out_F_after_0,"%20d %20.10f %20.10f %20.10f\n", tag[i], f_extra[i][0], f_extra[i][1], f_extra[i][2]);
    }
    fprintf(out_F_after_0,"\n");
    fclose(out_F_after_0);

    FILE *out_eatom_after_0 = fopen("eatom_after_0", "a");
    for (int i = 0; i < nlocal; i++) {
        fprintf(out_eatom_after_0,"%20d %20.10f %20.10f\n", tag[i], *(obj_kspace.eatom +i), *(obj_CoulLong.eatom +i));
    }
    fclose(out_eatom_after_0);
     */

    //reset the values set for the calculation
    atom->q = q;    //makes atom->q point at the real charge array

    atom->f = f;    //makes atom->f point at the real force array

    //deleting values of dummy array for force
    for (int i = 0; i < nlocal; i++) {
        delete [] f_extra[i];
    }

    delete [] f_extra;

}

// the below function calculates the forces on the ions if the charges of the wall particles are set to
// transpose(H)*(AQ + Bq) or -transpose(Z)*(Q + inv(A)*Bq)
void FixConp::force_extra_ions_2() {

    int *tag = atom->tag;
    int nlocal = atom -> nlocal;
    double *q = atom->q;
    double **f = atom->f;
    int nall = atom->nlocal + atom->nghost;

    //setting the values of the wall charges

    //declaration and initialisation of needed variables
    double q_for_cal[nall] = {0.0};
    int elealli;
    int ij;
    int ji;
    double AinvBq[elenum_all] = {0.0};
    double Q_plus_AinvBq[elenum_all] = {0.0};


    FILE *out_Bq=fopen("Bq_inf_extra", "a");
    for (int i = 0; i < elenum_all; i++) {
        fprintf(out_Bq, "%20.10f\n", Bq_CP4M2[i]);
    }
    fclose(out_Bq);

    //calculation of A^{-1}Bq
    for (int i = 0; i < elenum_all; i++){
        for (int j = 0; j < elenum_all; j++){
            ij = i*elenum_all + j;
            AinvBq[i] += aaa_all[ij]*Bq_CP4M2[j];
        }
    }

    /*
    FILE *out_AinvBq=fopen("AinvBq", "a");
    for (int i = 0; i < elenum_all; i++) {
        fprintf(out_AinvBq, "%20.10f\n", AinvBq[i]);
    }
    fclose(out_AinvBq);
     */

    //calculation of Q + A^{-1}Bq
    for (int i = 0; i < nall; i++){
        if (electrode_check(i)) {
            elealli = tag2eleall[tag[i]];
            Q_plus_AinvBq[elealli] = q[i] + AinvBq[elealli];
        }
    }

    /*
    FILE *out_Q_plus_AinvBq=fopen("Q_plus_AinvBq", "a");
    for (int i = 0; i < elenum_all; i++) {
        fprintf(out_Q_plus_AinvBq, "%20.10f\n", Q_plus_AinvBq[i]);
    }
    fclose(out_Q_plus_AinvBq);


    FILE *out_q_for_cal_0=fopen("q_for_cal_0", "a");
    for (int i = 0; i < nall ; i++) {
        fprintf(out_q_for_cal_0, "%20d %20.10f %20.10f\n", tag[i], q[i], q_for_cal[i]);
    }
    fclose(out_q_for_cal_0);
     */

    //calculation of -Z*[Q + A^{-1}Bq] and setting q_for_cal to it
    for (int i = 0; i < nall; i++) {
        if (electrode_check(i)) {
            elealli = tag2eleall[tag[i]];
            for (int j = 0; j < elenum_all; j++){
                ji = j*elenum_all+elealli;
                q_for_cal[i] += -Z[ji]*Q_plus_AinvBq[j];     //charges of the wall particles are set to -transpose(Z)*(Q + inv(A)*Bq)
            }
        }
        else{
            q_for_cal[i] = q[i];       //the charges of the ions are set to their real values
        }
    }

    atom->q = &q_for_cal[0];            // q now points the values of q_for_cal

    //dummy array for forces
    double **f_extra;
    f_extra = new double*[nlocal];

    for (int i = 0; i < nlocal; i++) {
        f_extra[i] = new double[3] ;
    }

    for (int i = 0; i < nlocal; i++) {
        f_extra[i][0] = 0;
        f_extra[i][1] = 0;
        f_extra[i][2] = 0;
    }

    atom->f = f_extra;        //makes atom->f point to dummy array


    FILE *out_q_for_cal=fopen("q_for_cal", "a");
    for (int i = 0; i < nall ; i++) {
        fprintf(out_q_for_cal, "%20d %20.10f %20.10f\n", tag[i], q[i], q_for_cal[i]);
    }
    fprintf(out_q_for_cal,"break\n");
    fclose(out_q_for_cal);

    FILE *out_F_before_2 = fopen("F_before_2", "a");
    for (int i = 0; i < nlocal; i++) {
        fprintf (out_F_before_2,"%20d %20.10f %20.10f %20.10f\n", tag[i], f_extra[i][0], f_extra[i][1], f_extra[i][2]);
    }
    fprintf(out_F_before_2,"\n");
    fclose(out_F_before_2);

    obj_kspace.compute(3,1);
    obj_CoulLong.compute(3,1);

    FILE *out_F_after_H_AQ_plus_Bq = fopen("F_after_H_AQ_plus_Bq", "a");
    for (int i = 0; i < nlocal; i++) {
        fprintf (out_F_after_H_AQ_plus_Bq,"%20d %20.10f %20.10f %20.10f\n", tag[i], f_extra[i][0], f_extra[i][1], f_extra[i][2]);
    }
    fprintf(out_F_after_H_AQ_plus_Bq,"\n");
    fclose(out_F_after_H_AQ_plus_Bq);

    //added f_extra to forces on ions
    for (int i = 0; i < nlocal; i++) {
        if (electrode_check(i)==0){
            f[i][0] = f[i][0] + f_extra[i][0];
            f[i][1] = f[i][1] + f_extra[i][1];
            f[i][2] = f[i][2] + f_extra[i][2];
        }
    }

    //reset the values set for the calculation
    atom->q = q;    //makes atom->q point at the real charge array

    atom->f = f;    //makes atom->f point at the real force array

    //deleting values of dummy array for force
    for (int i = 0; i < nlocal; i++) {
        delete [] f_extra[i];
    }

    delete [] f_extra;

}


void FixConp::V_cal() {

    //double Q_sum = 0;
    //double V_p = 0;
    //double V_m = 0;
    double Ep_dot_Ainv_Bq = 0;
    double Em_dot_Ainv_Bq = 0;
    double AinvBq[elenum_all] = {0.0};
    double nlocal = atom->nlocal;
    int *tag = atom->tag;
    double *q = atom->q;


    //calculation of A^{-1}Bq
    int ij = 0;
    for (int i = 0; i < elenum_all; i++){
        for (int j = 0; j < elenum_all; j++){
            ij = i*elenum_all + j;
            AinvBq[i] += aaa_all[ij]*Bq_CP4M2[j];
        }
    }

    int elealli = 0;
    for (int i = 0; i < nlocal; i++){
        if (electrode_check(i) == 1) {
            elealli = tag2eleall[tag[i]];
            Em_dot_Ainv_Bq += AinvBq[elealli];
        }
        if(electrode_check(i) == -1){
            elealli = tag2eleall[tag[i]];
            Ep_dot_Ainv_Bq += AinvBq[elealli];
            //Q_sum += q[i];

        }
    }

    V_p = C_pp*(Q_sum + Ep_dot_Ainv_Bq) + C_pm*(-Q_sum + Em_dot_Ainv_Bq);
    V_m = C_mp*(Q_sum + Ep_dot_Ainv_Bq) + C_mm*(-Q_sum + Em_dot_Ainv_Bq);


    //FILE *out_V_cal = fopen("V_cal", "a");
    //fprintf(out_V_cal, "%20.10f %20.10f %20.10f %20.10f %20.10f %20.10f %20.10f %20.10f %20.10f\n", C_mm, C_mp, C_pm, C_pp, Ep_dot_Ainv_Bq, Em_dot_Ainv_Bq, Q_sum, V_m, V_p);
    //fclose(out_V_cal);

}

void FixConp::Q_cal() {
    Q_sum = 0;
    int nlocal = atom->nlocal;
    double *q = atom->q;
    for (int i = 0; i < nlocal; i++){
        if(electrode_check(i) == -1){
            Q_sum += q[i];
        }
    }
}



