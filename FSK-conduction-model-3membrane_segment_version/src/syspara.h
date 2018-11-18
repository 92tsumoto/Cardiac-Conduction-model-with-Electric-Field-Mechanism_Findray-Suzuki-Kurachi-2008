//#ifndef __SYSPARA_H_INCLUDE 
//#define __SYSPARA_H_INCLUDE

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "mkl.h"

#define NN 24
#define BUF 200
#define NUM 20
#define MEDIA_SITE 64
#define MEDIA_PATCH 3
#define MAT_SIZE MEDIA_SITE*MEDIA_PATCH
#define beats 10

#define R 8314.472
#define F 96485.33771638995
#define T 310

#define VNMAX 400*5+1
#define NAIMAX 50*5+1
#define KIMAX 200*5+1
#define dvm 5

struct varstruct {

	int datas;
	int line_wid[NUM];
	
	int n;
	double Rmyo,Gj,D,Istim;
	double coef,dist;
	double Ri,Rg,Rd,Rj;
	double gamma[MEDIA_SITE-1],delta[MEDIA_SITE-1];
	double p1,p2,p3,q1,q2;
	double ri,rd,rj;
	double w1,w2,w3,w4;
	double s1,s2[MEDIA_PATCH];
	double scale1,Srate;
	double RTonF;

// Na channel localize switch
	int na_switch;
	double narate,jmrate1,lmrate1,jmrate2,lmrate2;
	//double vshift;

// state variables
	double x0[NUM][NN];
	//double u[NN][MEDIA_SITE][MEDIA_PATCH];
	double ***u;
	double inner_v[MEDIA_SITE][MEDIA_PATCH];
	double inner_current[MEDIA_SITE][MEDIA_PATCH];
	double disk_membrane[MEDIA_SITE][2];
	double dvdt[MEDIA_SITE][MEDIA_PATCH];

// quantifies of the cleft potential
	double cleft_width,cleft_R;
	double cleft_axial[MEDIA_SITE-1][2],cleft_potential[MEDIA_SITE-1];

// Cell Geometry
	double RGC;
	double length,a;
	double vcell[MEDIA_PATCH],ageo[MEDIA_PATCH];
	double acap[MEDIA_PATCH],vmyo[MEDIA_PATCH];
	double vmito[MEDIA_PATCH],vsr[MEDIA_PATCH];
	double vnsr[MEDIA_PATCH],vjsr[MEDIA_PATCH];
	double vcleft[MEDIA_PATCH];

// Ion Valences 
	double zna,zk,zca;

// Ion concentrations
	double nao,ko,cao;

// Tablization of exp functions
	// Fast sodium current
		double *Tam,*Tbm,*Tah,*Tbh,*Taj,*Tbj,*TEna;
	
	// L-type calcium current
		double *TexpCa,*Tdss,*Ttaud;              
		double *Tfss1,*Tfss2,*Ttauf1,*Ttauf2;
		double *Tfcass1,*Tfcass2,*Tfcatau1,*Tfcatau2;

	// T-type calcium current
		double *Tbss,*Ttaub,*Tgss,*Ttaug;

	// Rapidly activating potassium current
		double *Txrss,*Ttauxr,*Tr,*TEk;
	
	// Slowly activating potassium current
		double *Txs1ss,*Ttauxs1;
	
	// Inward rectificating potassium current
		double *T1Ki,*T2Ki,*T3Ki,*T4Ki;
		double *T5Ki,*T6Ki,*T7Ki,*T8Ki;

	// Plateau Potassium current
		double *Tkp;
	
	// Transient outward current
		double *Trvdv,*Tazdv,*Tbzdv,*Taydv,*Tbydv;

	// using Inak and Inaca
		double *Texp0,*Texp1,*Texp2,*Texp3,*Tnak;

	// using Irel and Irel_ol
		double *Tc,*Tol;

// Change rates for K channel conductances
	double ikrf[MEDIA_SITE][MEDIA_PATCH],iksf[MEDIA_SITE][MEDIA_PATCH];
	double ikxf[MEDIA_SITE][MEDIA_PATCH],itof[MEDIA_SITE][MEDIA_PATCH];

// fiber configuration switch
	int apd_switch;

// Fast sodium current
	double Ena[MEDIA_SITE][MEDIA_PATCH],ina[MEDIA_SITE][MEDIA_PATCH];
	double gna_local[MEDIA_SITE][MEDIA_PATCH],gna_all,gna,am,bm,ah,bh,aj,bj;
	double *total_ina;

// L-type calcium current
	double ilcatot,ibarca[MEDIA_SITE][MEDIA_PATCH],ibarna[MEDIA_SITE][MEDIA_PATCH],ibark[MEDIA_SITE][MEDIA_PATCH];
	double ilcarat[MEDIA_SITE][MEDIA_PATCH];
	double ilca[MEDIA_SITE][MEDIA_PATCH],ilcana[MEDIA_SITE][MEDIA_PATCH],ilcak[MEDIA_SITE][MEDIA_PATCH];
	double dss,taud;                  // Steady-state value of activation gate d and time const
	double fss1,fss2,tauf1,tauf2;     // Steady-state value of inactivation gate f and time const
	double fcass1,fcass2,fcatau1,fcatau2; // Ca dependant inactivation gate
	double kmca,pca,gacai,gacao;
	double pna,ganai,ganao;           // Permiability of membrane to Na (cm/s)
	double pk,gaki,gako;              // Permiability of membrane to K (cm/s)
	double ratgca,fcarat,aCDI,bCDI;

// T-type calcium current
	double gcat,bss,taub,gss,taug;
	double Eca[MEDIA_SITE][MEDIA_PATCH],icat[MEDIA_SITE][MEDIA_PATCH];

// Rapidly activating potassium current
	double gkr,xrss,tauxr,r,gkr_max;
	double Ekr[MEDIA_SITE][MEDIA_PATCH],ikr[MEDIA_SITE][MEDIA_PATCH];
	
// Slowly activating potassium current
	double gks,gks_max,xs1ss,tauxs1,xs2ss,tauxs2;
	double prnak,iks[MEDIA_SITE][MEDIA_PATCH],Eks[MEDIA_SITE][MEDIA_PATCH];

// Potassium current IK1 (time-dependent)
	double gki,aki,bki,kin,c1_ki,c2_ki,c3_ki,c4_ki;
	double Eki[MEDIA_SITE][MEDIA_PATCH],iki[MEDIA_SITE][MEDIA_PATCH];

// Plateau Potassium current
	double gkp,gkp_max,kp;
	double Ekp[MEDIA_SITE][MEDIA_PATCH],ikp[MEDIA_SITE][MEDIA_PATCH];

// Ultra rapid Potassium current
	double gkur,ikur[MEDIA_SITE][MEDIA_PATCH];

// Na-activated potassium channel
	double gkna,pona,pov,kdkna;
	double Ekna[MEDIA_SITE][MEDIA_PATCH],ikna[MEDIA_SITE][MEDIA_PATCH];

// ATP-sensitive potassium channel
	double gkatp,gkbaratp,patp,natp;
	double nicholsarea,atpi,hatp,katp;
	double Ekatp[MEDIA_SITE][MEDIA_PATCH],ikatp[MEDIA_SITE][MEDIA_PATCH];

// Ito Transient Outward Current
// (Dumaine et al. Circ Res 1999;85:803-809)
	double gitodv,gitodv_max,rvdv;
	double azdv,bzdv,aydv,bydv;
	double Ekdv[MEDIA_SITE][MEDIA_PATCH],ito[MEDIA_SITE][MEDIA_PATCH];

// Sodium-Calcium Exchanger V-S
	double c1,c2,gammas,inaca[MEDIA_SITE][MEDIA_PATCH];

// Sodium-Potassium Pump
	double fnak,sigma,ibarnak,kmnai,kmko;
	double inak[MEDIA_SITE][MEDIA_PATCH];
	
// Nonspecific Ca-activated Current 
	double ibarnsna,ibarnsk,pnsca,kmnsca;
	double insna[MEDIA_SITE][MEDIA_PATCH],insk[MEDIA_SITE][MEDIA_PATCH];

// Sarcolemmal Ca Pump
	double ibarpca,kmpca,ipca[MEDIA_SITE][MEDIA_PATCH];

// Ca Background Current
	double gcab,Ecan[MEDIA_SITE][MEDIA_PATCH],icab[MEDIA_SITE][MEDIA_PATCH];

// Na Background Current
	double gnab,Enan[MEDIA_SITE][MEDIA_PATCH],inab[MEDIA_SITE][MEDIA_PATCH];

// Total Ion currents 
	double INa_total[MEDIA_SITE][MEDIA_PATCH];
	double IK_total[MEDIA_SITE][MEDIA_PATCH];
	double ICa_total[MEDIA_SITE][MEDIA_PATCH];
	
// Difference total Ion current 
	double dICa_total[MEDIA_SITE][MEDIA_PATCH];
	
// NSR Ca Ion Concentration Changes
	double iup[MEDIA_SITE][MEDIA_PATCH],ileak[MEDIA_SITE][MEDIA_PATCH];      // Ca uptake from myo. to NSR (mM/ms)
	double kmup,iupbar,nsrbar;

// JSR Ca Ion Concentration Changes
	double Irel_cicr[MEDIA_SITE][MEDIA_PATCH],Irel_jsr_overload[MEDIA_SITE][MEDIA_PATCH];
	double dICa_total_new[MEDIA_SITE][MEDIA_PATCH],dCa_ion[MEDIA_SITE],ICa_total_old[MEDIA_SITE][MEDIA_PATCH];
	double boolien[MEDIA_SITE][MEDIA_PATCH],gmaxrel;
	double jsr_new,grelbarjsrol,tauon,tauoff;
	double csqnbar,kmcsqn,t_cicr[MEDIA_SITE][MEDIA_PATCH],t_overload[MEDIA_SITE][MEDIA_PATCH];
	double bjsr[MEDIA_SITE][MEDIA_PATCH],cjsr[MEDIA_SITE][MEDIA_PATCH],csqn[MEDIA_SITE][MEDIA_PATCH],csqnth,swspontan;
	double djsr,jsr;

// test variable
	double dt;

// Translocation of Ca Ions from NSR to JSR
	double tautr,itr[MEDIA_SITE][MEDIA_PATCH];

// Ca concentration
	double cmdnbar,trpnbar,kmtrpn,kmcmdn;
	double trpn[MEDIA_SITE][MEDIA_PATCH],cmdn[MEDIA_SITE][MEDIA_PATCH];
	double Ca_total[MEDIA_SITE][MEDIA_PATCH],gpig;
	double bmyo[MEDIA_SITE][MEDIA_PATCH],cmyo[MEDIA_SITE][MEDIA_PATCH],dmyo[MEDIA_SITE][MEDIA_PATCH];

// Extracellular ion concentrations
	double Na_bulk,K_bulk,Ca_bulk;
	double tau_diff;

// Base Currnt Stimulus
	double Istim_base;

// Sttimulus parameters
	double BCL;  // Base cycle length = stimulus period
	int beat; // Number of stimulus

// debug variable
	double ca_pre[MEDIA_SITE],dca_now;

    int m;
    int l;

    double tsign[NUM];
    double tend[NUM];

    int write;
    int half;
	int out_data;

// linear algebra
	double *AT,*AT2;
// with Pardiso rutine
	double *b, *solv;
	long long *ia, *ja, *stok_ia, *stok_ja;
	int par_k;
	long long mtype,nrhs,iparm[64],maxfct,mnum,phase;
	long long error, msglvl;
	double ddum;
	long long idum;

/* Internal solver memory pointer pt, */
/* 32-bit: int pt[64]; 64-bit: long int pt[64] */
/* or void *pt[64] should be OK on both architectures */
	void *pt[64];

// using Eular
//double k1[NN][MEDIA_SITE][MEDIA_PATCH];
	double ***k1;

} var;

// for main
void make_ExPTable();
void eular(int n, double h, double t);
void function(int site, int patch, double t);
void input_para(FILE *);
void static_paras(FILE *);
void mem();
void close_mem();

// for dataout
void vm_data(FILE *, double time);
void ina_data(FILE *, double time);
void ik_data(FILE *, FILE *, FILE *, FILE *, FILE *, double time);
void ica_data(FILE *, FILE *, FILE *, FILE *, double time);
void pump_data(FILE *, FILE *, double time);
void cai_data(FILE *, double time);
void intra_v_data(FILE *, double time);
//void itotal_data(FILE *, double time);
//void intra_i_data(FILE *, double time);
//void cleft_v_data(FILE *, FILE *, double time);
//void dvdt_data(FILE *, double time);

// for calculation of Ion currents
void comp_ina(int site, int patch);
void comp_icat(int site, int patch);
void comp_ical(int site, int patch);
void comp_ikr(int site, int patch);
void comp_iki(int site, int patch);
void comp_iks(int site, int patch);
void comp_ikp(int site, int patch);
void comp_ikur(int site, int patch);
void comp_ikna(int site, int patch);
void comp_ikatp(int site, int patch);
void comp_ito(int site, int patch);
void comp_inaca(int site, int patch);
void comp_inak(int site, int patch);
void comp_insca(int site, int patch);
void comp_ipca(int site, int patch);
void comp_icab(int site, int patch);
void comp_inab(int site, int patch);
void comp_iconcent (int site, int patch);
void comp_iconcent2 (int site, int patch);
void conc_nsr(int site, int patch);
void conc_jsr(int site, int patch);
void conc_itr(int site, int patch);
void conc_cai(int site, int patch);
//void conc_cleft(int site, int patch);

MKL_INT init_pardiso();
void linear_coeff();
void linear_vec();
void linear_vec3();
void linear_solve();
void release_of_memory();
//void dgesv_(int *N, int *NRHS, double AT[], int *LDA, int pivot[], double b[], int *LDB, int *ok);

//void main(int argc, char **argv);

