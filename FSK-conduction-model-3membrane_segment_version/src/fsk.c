#include <string.h>
#include "syspara.h"

int mode = 1;
int P = 8;
FILE *fopen(), *fpin, *fp0, *fp1, *fp2, *fp3, *fp4, *fp5;
FILE *fp6, *fp7, *fp8, *fp9, *fp10, *fp11, *fp12, *fp13, *fp14, *fp15;
FILE *fp16, *fp17, *fp18, *fp19, *fp20, *fp21, *fp22, *fp23, *fp24, *fp25;

typedef double Number;
typedef long long Lint;

main(int argc, char **argv)
{
	int i,k,m,l,z;
	int ii=0;
	unsigned int count=0;
	double t = 0.0;
	double time;
	double h,R_all,V_all;
	//double v_old[MEDIA_SITE][MEDIA_PATCH],v_old2[MEDIA_SITE][MEDIA_PATCH],dvdt_new[MEDIA_SITE][MEDIA_PATCH];
	double **v_old,**v_old2,**dvdt_new;
	double d1,d2;
	double v_thre;
	double cut;
	char *tmpname;
	char cmd[BUFSIZ];
	double tend;

/* Action Potential Duration and Max. Info */
	double ***vmax; // Max. Voltage (mV)
	double ***dvdtmax; // Max. dv/dt (mV/ms)
	double ***apd; // Action Potential Duration
	double ***toneapd; // Time of dv/dt Max.
	double ***ttwoapd; // Time of 90% Repolarization
	double ***rmbp; // Resting Membrane Potential
	double ***nair; // Intracellular Na At Rest
	double ***cair; // Intracellular Ca At Rest
	double ***kir; // Intracellular K At Rest
	double ***caimax; // Peak Intracellular Ca
	double ***inamax; // Peak INa
	double **total_inamax; // Peak total INa
	double *total_ina; // Peak total INa within a cell
	double v_thre_time[beats][MEDIA_SITE][MEDIA_PATCH] ; // time when AP passes throght the threshold value
	double v_thre_time2[beats][MEDIA_SITE][MEDIA_PATCH] ; // time when AP passes throght the threshold value
	int thre_flag[MEDIA_SITE][MEDIA_PATCH] ; // time flag when AP passes throght the threshold value
	int thre_flag2[MEDIA_SITE][MEDIA_PATCH] ; // time flag when AP passes throght the threshold value


	var.b=(Number *)calloc(MAT_SIZE,sizeof(Number));
	var.solv=(Number *)calloc(MAT_SIZE,sizeof(Number));
	var.AT=(Number *)calloc(MAT_SIZE*MAT_SIZE,sizeof(Number));

	tmpname = "temp";

	sprintf(cmd, "/usr/bin/cpp -P %s > %s", argv[1],tmpname);
	if(system(cmd) == -1){
		fprintf(stderr,"cannot open %s\n",argv[1]);
		exit(1);
	}
	if((fpin=fopen(tmpname,"r"))==NULL){
		fprintf(stderr,"cannot open %s\n",argv[1]);
		exit(1);
	}
	if ((fp1 = fopen("status.out","w")) == NULL){
		printf("Can't open File\n");
		exit(1);
	}
	if ((fp2 = fopen("ndata.out","w")) == NULL){
		printf("Can't open File\n");
		exit(1);
	}
	if ((fp3 = fopen("ndata_final.out","w")) == NULL){
		printf("Can't open File\n");
		exit(1);
	}
	if ((fp4 = fopen("check_ndata_final.out","w")) == NULL){
		printf("Can't open File\n");
		exit(1);
	}
	if ((fp5 = fopen("act_time.out","w")) == NULL){
		printf("Can't open File\n");
		exit(1);
	}

	if ((fp6 = fopen("vm_data.out","w")) == NULL){
		printf("Can't open File\n");
		exit(1);
	}
	
	if ((fp22 = fopen("initdat.out","w")) == NULL){
		printf("Can't open File\n");
		exit(1);
	}

	// input initial parameters
	printf("start parameter input\n");
	input_para(fpin);
	printf("end parameter input\n");

	if (var.write){
		if ((fp0 = fopen(argv[2],"w"))==NULL){
			fprintf(stderr, "%s cannot open.\n",argv[2]);
			exit(-1);
		}
	}

	for (ii = 0; ii < var.datas; ii++){
		long j;
		time = 0.0;
		tend = var.tend[ii];

		// Eular's step size
		h = 1.0 / var.m;
		// AP threshold value
		v_thre = -30.0;
		// a variable for data output	
		cut = var.m/10.0;

		h *= var.tsign[ii];

	// initialized memory
		mem();
        v_old = (double**)malloc(sizeof(Number)*MEDIA_SITE);
        v_old2 = (double**)malloc(sizeof(Number)*MEDIA_SITE);
        dvdt_new = (double**)malloc(sizeof(Number)*MEDIA_SITE);
		if(v_old==NULL || v_old2==NULL || dvdt_new==NULL) exit(1);
        for(i=0;i<MEDIA_SITE;i++){
                v_old[i] = (Number *)malloc(sizeof(Number)*MEDIA_PATCH);
                v_old2[i] = (Number *)malloc(sizeof(Number)*MEDIA_PATCH);
                dvdt_new[i] = (Number *)malloc(sizeof(Number)*MEDIA_PATCH);
		if(v_old[i]==NULL || v_old2[i]==NULL || dvdt_new[i]==NULL) exit(1);
        }

        vmax = (double***)calloc(beats, sizeof(double**));
        dvdtmax = (double***)calloc(beats, sizeof(double**));
        apd = (double***)calloc(beats, sizeof(double**));
        toneapd = (double***)calloc(beats, sizeof(double**));
        ttwoapd = (double***)calloc(beats, sizeof(double**));
        rmbp = (double***)calloc(beats, sizeof(double**));
        nair = (double***)calloc(beats, sizeof(double**));
        cair = (double***)calloc(beats, sizeof(double**));
        kir = (double***)calloc(beats, sizeof(double**));
        caimax = (double***)calloc(beats, sizeof(double**));
        inamax = (double***)calloc(beats, sizeof(double**));
        total_inamax = (double**)calloc(beats, sizeof(double*));
        total_ina = (double*)calloc(MEDIA_SITE, sizeof(double));
        var.total_ina = (double*)calloc(MEDIA_SITE, sizeof(double));
        if (vmax==NULL || dvdtmax==NULL || apd==NULL || toneapd==NULL
                || ttwoapd==NULL || rmbp==NULL || nair==NULL || cair==NULL
                || kir==NULL || caimax==NULL || inamax==NULL || total_inamax==NULL 
				|| total_ina==NULL || var.total_ina==NULL ) exit(1);
        for (i=0;i<beats;i++){
                vmax[i] = (double**)calloc(MEDIA_SITE, sizeof(double*));
                dvdtmax[i] = (double**)calloc(MEDIA_SITE, sizeof(double*));
                apd[i] = (double**)calloc(MEDIA_SITE, sizeof(double*));
                toneapd[i] = (double**)calloc(MEDIA_SITE, sizeof(double*));
                ttwoapd[i] = (double**)calloc(MEDIA_SITE, sizeof(double*));
                rmbp[i] = (double**)calloc(MEDIA_SITE, sizeof(double*));
                nair[i] = (double**)calloc(MEDIA_SITE, sizeof(double*));
                cair[i] = (double**)calloc(MEDIA_SITE, sizeof(double*));
                kir[i] = (double**)calloc(MEDIA_SITE, sizeof(double*));
                caimax[i] = (double**)calloc(MEDIA_SITE, sizeof(double*));
                inamax[i] = (double**)calloc(MEDIA_SITE, sizeof(double*));
                total_inamax[i] = (double*)calloc(MEDIA_SITE, sizeof(double));
                if (vmax[i]==NULL || dvdtmax[i]==NULL || apd[i]==NULL || toneapd[i]==NULL
                        || ttwoapd[i]==NULL || rmbp[i]==NULL || nair[i]==NULL || cair[i]==NULL
                        || kir[i]==NULL || caimax[i]==NULL || inamax[i]==NULL || total_inamax[i]==NULL ) exit(1);
                for(k=0;k<MEDIA_SITE;k++){
                        vmax[i][k] = (double*)calloc(MEDIA_PATCH, sizeof(double));
                        dvdtmax[i][k] = (double*)calloc(MEDIA_PATCH, sizeof(double));
                        apd[i][k] = (double*)calloc(MEDIA_PATCH, sizeof(double));
                        toneapd[i][k] = (double*)calloc(MEDIA_PATCH, sizeof(double));
                        ttwoapd[i][k] = (double*)calloc(MEDIA_PATCH, sizeof(double));
                        rmbp[i][k] = (double*)calloc(MEDIA_PATCH, sizeof(double));
                        nair[i][k] = (double*)calloc(MEDIA_PATCH, sizeof(double));
                        cair[i][k] = (double*)calloc(MEDIA_PATCH, sizeof(double));
                        kir[i][k] = (double*)calloc(MEDIA_PATCH, sizeof(double));
                        caimax[i][k] = (double*)calloc(MEDIA_PATCH, sizeof(double));
                        inamax[i][k] = (double*)calloc(MEDIA_PATCH, sizeof(double));
                        if(vmax[i][k]==NULL || dvdtmax[i][k]==NULL || apd[i][k]==NULL || toneapd[i][k]==NULL
                                || ttwoapd[i][k]==NULL || rmbp[i][k]==NULL || nair[i][k]==NULL || cair[i][k]==NULL
                                || kir[i][k]==NULL || caimax[i][k]==NULL || inamax[i][k]==NULL ) exit(1);
                }
        }

		printf("finished memory initialization\n");

	// initial values input.
		for (i = 0; i < NN; i++){ 
			for (m = 0; m < MEDIA_SITE; m++){ 
				for (l = 0; l < MEDIA_PATCH; l++){ 
					var.u[i][m][l] = var.x0[ii][i];
					var.inner_v[m][l] = var.u[0][m][l];
				}
			}
		}


	// input static parameters
	static_paras(fp1);
	printf("finished input of static parameters\n");

	if(var.out_data){
		if ((fp7 = fopen("ina_data.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		} 
		if ((fp8 = fopen("ito_data.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
		if ((fp9 = fopen("ikr_data.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
		if ((fp10 = fopen("iks_data.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
		if ((fp11 = fopen("iki_data.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
		if ((fp12 = fopen("ikp_data.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
		if ((fp13 = fopen("ical_data.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
		if ((fp14 = fopen("icat_data.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
		if ((fp15 = fopen("icab_data.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
		if ((fp16 = fopen("ipca_data.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		} 
		if ((fp17 = fopen("inak_data.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
		if ((fp18 = fopen("incx_data.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
		if ((fp19 = fopen("cai_data.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
		if ((fp20 = fopen("vi_data.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
	/*	if ((fp23 = fopen("total_inamax.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
		if ((fp24 = fopen("whole_ina.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		} */
	}	

	// Resistance parameters;
		d1 = var.length/(MEDIA_PATCH-2);
		d2 = var.cleft_width;

	// surface area of Junctional membran unit or cross-section area of myocyte	
		var.s1 = M_PI*var.a*var.a;

	// surface area of junctional and non-Junctional membran unit	
		R_all = 0; V_all = 0.0;
		for(i=0;i<MEDIA_PATCH;i++){
			if (i==0 || i== MEDIA_PATCH-1){
				//var.s2[i] = var.Srate*(M_PI*var.a*var.a);
				var.s2[i] = M_PI*var.a*var.a;
			} else {
				var.s2[i] = 2.0*M_PI*var.a*d1;
			}
			R_all+=var.s2[i];
			printf("s2[%d]=%e\n",i,var.s2[i]);
			fprintf(fp1,"s2[%d]=%e\n",i,var.s2[i]);
		}		
		for(i=0;i<MEDIA_PATCH;i++){
			V_all+=var.vcell[i];
			printf("vcell[%d]=%e\n",i,var.vcell[i]);
			fprintf(fp1,"vcell[%d]=%e\n",i,var.vcell[i]);
		}		
		printf("total_cell_volume=%e:%e\n",V_all,1000*M_PI*var.a*var.length*var.a);
		fprintf(fp1,"total_cell_volume=%e:%e\n",V_all,1000*M_PI*var.a*var.length*var.a);
		printf("total_Surface=%e:%e\n",R_all,2*M_PI*var.a*var.length+2*M_PI*var.a*var.a);
		fprintf(fp1,"total_Surface=%e:%e\n",R_all,2*M_PI*var.a*var.length+2*M_PI*var.a*var.a);
	
		var.Ri = (var.Rmyo*d1)/var.s1;
		printf("Rmyo=%e, Ri=%e\n",var.Rmyo,var.Ri);
		fprintf(fp1,"Rmyo=%e, Ri=%e\n",var.Rmyo,var.Ri);
		
		var.Rg = var.Gj*var.length/var.s1;
		printf("Rg=%e\n",var.Rg);	fprintf(fp1,"Rg=%e\n",var.Rg);

		var.Rj = var.cleft_R;
		//var.Rd = (var.Rmyo*(var.Rmyo/(8*M_PI*var.cleft_R)))/var.s1;
		var.Rd = (150.0*(150.0/(8*M_PI*var.cleft_R)))/var.s1;
		printf("Rd=%e, Rj=%e, cw=%e\n",var.Rd,var.Rj,150.0/(8*M_PI*var.Rj));
		fprintf(fp1,"Rd=%e, Rj=%e\n",var.Rd,var.Rj);
	
		var.p1 = var.Rj+0.5*(var.Ri+var.Rd);
		var.p2 = var.Rj+0.5*var.Rd;
		var.p3 = var.Rj;
		
		var.q1 = 1.0+var.Ri/var.Rg/2.0;
		var.q2 = var.Ri/var.Rg/2.0;

		printf("p1=%e, p2=%e, p3=%e\n",var.p1,var.p2,var.p3);
		printf("q1=%e, q2=%e\n",var.q1,var.q2);
		printf("Rj+0.5(Ri+Rd)=%e, Rj+0.5Rd=%e,Rj=%e\n",var.p1,var.p2,var.p3);
		printf("1+Ri/Rg/2=%e, Ri/Rg/2=%e\n",var.q1,var.q2);
		printf("radius = %e\n",var.a);
		printf("unit length = %e\n",var.length);
		printf("unit surface area 2*Pi*a*dx = %e\n",2.0*M_PI*var.a*var.length);
		printf("Cross section area Pi*a*a = %e\n",M_PI*var.a*var.a);
		printf("intracellular resistance = %e\n",var.Rmyo*var.length/(M_PI*var.a*var.a));
		printf("Ri = %e, (Gi = %e)\n",var.Rmyo*d1/(M_PI*var.a*var.a),M_PI*var.a*var.a/(var.Rmyo*var.length));
		printf("Rg = %e, (Gj = %e)\n",var.Gj*var.length/(M_PI*var.a*var.a),(M_PI*var.a*var.a)/var.length/var.Gj);
		printf("Ena=%f, Ek=%f\n",log(var.u[21][0][1]/var.u[16][0][1])*var.RTonF,log(var.u[22][0][1]/var.u[17][0][1])*var.RTonF);

		fprintf(fp1,"radius = %e\n",var.a);
		fprintf(fp1,"unit length = %e\n",var.length);
		fprintf(fp1,"unit surface area 2*Pi*a*dx = %e\n",2.0*M_PI*var.a*var.length);
		fprintf(fp1,"Cross section area Pi*a*a = %e\n",M_PI*var.a*var.a);
		fprintf(fp1,"intracellular resistance = %e\n",var.Rmyo*var.length/(M_PI*var.a*var.a));
		fprintf(fp1,"Ri = %e, (Gi = %e)\n",var.Rmyo*d1/(M_PI*var.a*var.a),M_PI*var.a*var.a/(var.Rmyo*var.length));
		fprintf(fp1,"Rg = %e, (Gg = %e)\n",var.Gj*var.length/(M_PI*var.a*var.a),(M_PI*var.a*var.a)/var.length*var.Gj);
		fprintf(fp1,"Rj+0.5(Ri+Rd)=%e, Rj+0.5Rd=%e,Rj=%e\n",var.p1,var.p2,var.p3);
		fprintf(fp1,"1+Ri/Rg/2=%e, Ri/Rg/2=%e\n",var.q1,var.q2);
		fprintf(fp1,"Ena=%f, Ek=%f\n",log(var.u[21][0][1]/var.u[16][0][1])*var.RTonF,log(var.u[22][0][1]/var.u[17][0][1])*var.RTonF);

	// for solving algebra equation.
		var.par_k = 9*(MEDIA_SITE-2)+6*2;
		printf("var.par_k = %d\n",var.par_k);
		printf("R matrix size = %d X %d\n",MAT_SIZE,MAT_SIZE);
		printf("k=%d, Mem_size=%d MB\n",var.par_k,((var.par_k)*(int)sizeof(Number))/1024/1024);
	// setting of each element of Resistance Matrix
		var.ia=(Lint *)calloc(var.par_k,sizeof(Lint));
		if(var.ia==NULL){
			printf(" ia null \n");
			exit(1);
		}
		var.ja=(Lint *)calloc(var.par_k,sizeof(Lint));
		if(var.ja==NULL){
			printf(" ja null \n");
			exit(1);
		}
		var.AT=(Number *)calloc(var.par_k,sizeof(Number));
		if(var.AT==NULL){
			printf(" AT null \n");
			exit(1);
		}
		printf("coeff input start\n");
		linear_coeff(); // setting of each element of Resistance Matrix
		printf("coeff input end\n");
		var.mtype = 11; /* Real unsymmetric matrix */
		var.nrhs = 1; /* Number of right hand sides. */
		init_pardiso();
		printf("pardiso initiation ok\n");

		var.b=(Number *)calloc(MAT_SIZE,sizeof(Number));
		if(var.b==NULL){
			printf(" B matrix is null \n");
			exit(1);
		}

		var.solv=(Number *)calloc(MAT_SIZE,sizeof(Number));
		if(var.solv==NULL){
			printf(" solution matrix is null \n");
			exit(1);
		}

	// Tablize exp functions.       
		printf("start tablization\n");
		make_ExpTable();
		printf("finished tablization\n");

	// Initialization time
		time -= h;
		var.dt = h;
		var.beat = 0;
	
		for(var.beat=0; var.beat < beats; var.beat++){
			if(var.beat==beats-1){
				var.l = 2.0*var.BCL;
			} else {
				var.l = var.BCL;
			}


			for (j = 0; j< (var.m * var.BCL ); j++){
				t = h*j;
				time += h;

				for(k=0;k<MEDIA_SITE;k++){
					for(m=0;m<MEDIA_PATCH;m++){
						v_old[k][m] = var.u[0][k][m];
						v_old2[k][m] = var.inner_v[k][m];
					}
				}
					
			/*	for(k=0;k<MEDIA_SITE;k++){
					for(m=0;m<MEDIA_PATCH;m++){
						if (m==0 && k >= 1){
							v_old2[k][m] = var.inner_v[k][0];
						} else if(m==MEDIA_PATCH-1 && k < MEDIA_SITE-1 ){
							v_old2[k][m] = var.inner_v[k][1];
						} else {
							v_old2[k][m] = var.u[0][k][m];
						}
					}
				} */
					
				if ( time-(var.BCL*var.beat+10) >= 0.0 && time-(var.BCL*var.beat+10) < h ){

					for(k=0;k<MEDIA_SITE;k++){
						for(m=0;m<MEDIA_PATCH;m++){
							var.boolien[k][m] = 0;
							rmbp[var.beat][k][m] = var.u[0][k][m];
							nair[var.beat][k][m] = var.u[16][k][m];
							kir[var.beat][k][m]  = var.u[17][k][m];
							cair[var.beat][k][m] = var.u[18][k][m];
							thre_flag[k][m] = 0;
							thre_flag2[k][m] = 0;
						}
					}
				}

				if (time-(var.BCL*(double)var.beat+10.0) >= 0.0 && time-(var.BCL*(double)var.beat+10.0) < 1.0){
					var.Istim = var.Istim_base/(var.s2[(int)floor(MEDIA_PATCH/2.0)]/(2*M_PI*var.a*var.a+2*M_PI*var.a*var.length));
				} else {
					var.Istim = 0.0;
				}

				if (fabs(time) > tend &&  tend != 0.0) break;

				eular(NN,h,t);

				for (k=0;k<MEDIA_SITE;k++){
					for (m=0;m<MEDIA_PATCH;m++){

						dvdt_new[k][m] = (var.u[0][k][m]-v_old[k][m])/h; // LRd -> dvdtnew
						var.dvdt[k][m] = dvdt_new[k][m];

						if(var.beat>=0){
							if (var.u[0][k][m] > vmax[var.beat][k][m] ){
								vmax[var.beat][k][m] = var.u[0][k][m];
							}
							if (var.u[18][k][m] > caimax[var.beat][k][m] ){
								caimax[var.beat][k][m] = var.u[18][k][m];
							}
							if (var.ina[k][m] < inamax[var.beat][k][m] ){
								inamax[var.beat][k][m] = var.ina[k][m];
							}
							if (dvdt_new[k][m] > dvdtmax[var.beat][k][m] ){
								dvdtmax[var.beat][k][m] = dvdt_new[k][m];
								toneapd[var.beat][k][m] = time;
							}
							if (var.u[0][k][m] >= (vmax[var.beat][k][m] -0.9*(vmax[var.beat][k][m] -rmbp[var.beat][k][m] ))){
								ttwoapd[var.beat][k][m] = time;
							}
							
							// chack for activation time
							/* if (m==0 && k >= 1 ){
								if( (v_old2[k][m] - v_thre)<0 && (v_old2[k][m] - v_thre)*(var.inner_v[k][0] - v_thre)<0){
									v_thre_time[var.beat][k][m] = time;
									//printf("pass thre unit[%d][%d] = %lf\n",k,m,v_thre_time[var.beat][k][m]);
								}
							} else if (m==MEDIA_PATCH-1 && k < MEDIA_SITE-1 ){
								if( (v_old2[k][m] - v_thre)<0 && (v_old2[k][m] - v_thre)*(var.inner_v[k][1] - v_thre)<0){
									v_thre_time[var.beat][k][m] = time;
									//printf("pass thre unit[%d][%d] = %lf\n",k,m,v_thre_time[var.beat][k][m]);
								}
							} else{
								if( (v_old2[k][m] - v_thre)<0 && (v_old2[k][m] - v_thre)*(var.u[0][k][m] - v_thre)<0){
									v_thre_time[var.beat][k][m] = time;
									///printf("pass thre unit[%d][%d] = %lf\n",k,m,v_thre_time[var.beat][k][m]);
								}
							}*/
							
							if( (v_old2[k][m] - v_thre)<0 && (v_old2[k][m] - v_thre)*(var.inner_v[k][m] - v_thre)<0){
								if(thre_flag[k][m] != 1){
									v_thre_time[var.beat][k][m] = time;
									//printf("pass thre unit[%d][%d] = %lf\n",k,m,v_thre_time[var.beat][k][m]);
									thre_flag[k][m] = 1;
								}
							}

							if( (v_old[k][m] - v_thre)<0 && (v_old[k][m] - v_thre)*(var.u[0][k][m] - v_thre)<0){
								if(thre_flag2[k][m] != 1){
									v_thre_time2[var.beat][k][m] = time;
									if(m==1){
										printf("pass thre unit[%d][%d] = %lf\n",k,m,v_thre_time2[var.beat][k][m]);
										fprintf(fp5,"%d %d %d %lf\n",var.beat,k,m,v_thre_time2[var.beat][k][m]);
									}
									thre_flag2[k][m] = 1;
								}
							}
							// check end for activation time

						}
					}
				}

				for (k=0;k<MEDIA_SITE;k++){
					for (m=0;m<MEDIA_PATCH;m++){
						if(var.csqn[k][m] >= var.csqnth && var.t_overload[k][m] > 50.0){
							//printf("reset t_overload at %lf\n",var.t_overload[k][m]);
							var.grelbarjsrol = 4;
							var.t_overload[k][m] = 0;
						}
					}
				}
				
				for (k=0;k<MEDIA_SITE;k++){
					for (m=0;m<MEDIA_PATCH;m++){
						total_ina[k]+=var.ina[k][m]*var.s2[m];
					}
					if (var.total_ina[k] < total_inamax[var.beat][k] ){
						total_inamax[var.beat][k] = var.total_ina[k];
					}
					var.total_ina[k]=total_ina[k];
					total_ina[k]=0;
				}

				if (time>= (beats-5)*var.BCL){
					if(j%(int)cut == 0){ // To reduce the recording data.
						vm_data(fp6,time);
						// option
						if(var.out_data){
							ina_data(fp7,time);
							ik_data(fp8,fp9,fp10,fp11,fp12,time);
							ica_data(fp13,fp14,fp15,fp16,time);
							pump_data(fp17,fp18,time);
							cai_data(fp19,time);
							intra_v_data(fp20,time);
						}

					}
				}
				
				for (k=0;k<MEDIA_SITE;k++){
					for (m=0;m<MEDIA_PATCH;m++){
						//dvdt[k][m] = dvdt_new[k][m];
						var.ICa_total_old[k][m] = var.ICa_total[k][m];
						var.dICa_total[k][m] = var.dICa_total_new[k][m];
					}
				}
				

			} // for-loop end; j

		// reset of total_ina
			for (k=0;k<MEDIA_SITE;k++){
				var.total_ina[k]=0.0;
			}

			fprintf(fp22,"beats=%d\n",var.beat+1);
			for(m=0;m<MEDIA_SITE;m++){
				for(k=0;k<MEDIA_PATCH;k++){
					for(i=0;i<NN;i++){
						fprintf(fp22,"%16.15e\n",var.u[i][m][k]);
					}
				}
			}
			fprintf(fp22,"\n");

		} // for-loop end; var.beats

		//printf("data out\n");
		
		for(z=0;z<beats;z++){

			for(m=0;m<MEDIA_SITE;m++){
				for(k=0;k<MEDIA_PATCH;k++){
					apd[z][m][k] = ttwoapd[z][m][k] - toneapd[z][m][k] ;
				}
			}

			for(m=0;m<MEDIA_SITE;m++){
				for(k=0;k<MEDIA_PATCH;k++){
					fprintf(fp2,"%i\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%e\n",
						z+1,vmax[z][m][k],dvdtmax[z][m][k],v_thre_time[z][m][k],v_thre_time2[z][m][k],apd[z][m][k],toneapd[z][m][k],ttwoapd[z][m][k],
							nair[z][m][k],kir[z][m][k],cair[z][m][k],caimax[z][m][k],rmbp[z][m][k],inamax[z][m][k]);
				}
			}

		/*	for(m=0;m<MEDIA_SITE;m++){
				fprintf(fp23,"%16.15e\n",total_inamax[z][m]/R_all);
			}*/	
		}
					
		for(m=0;m<MEDIA_SITE;m++){
			for(k=0;k<MEDIA_PATCH;k++){
				fprintf(fp4,"%d-%d-%d\n",beats-1,m,k);
				fprintf(fp3,"%i\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%e\n",
					beats,vmax[beats-1][m][k],dvdtmax[beats-1][m][k],v_thre_time[beats-1][m][k],v_thre_time2[beats-1][m][k],apd[beats-1][m][k],
						toneapd[beats-1][m][k],ttwoapd[beats-1][m][k],nair[beats-1][m][k],kir[beats-1][m][k],cair[beats-1][m][k],
							caimax[beats-1][m][k],rmbp[beats-1][m][k],inamax[beats-1][m][k]);
			}
		}


/* -------------------------------------------------------------------- */
/* .. Termination and release of memory. */
/* -------------------------------------------------------------------- */
		release_of_memory();

		fclose(fp1); fclose(fp2); fclose(fp3);
		fclose(fp4); fclose(fp5);
		fclose(fp6); fclose(fp22);
		
		if(var.out_data){
			fclose(fp7); fclose(fp8); fclose(fp9); fclose(fp10); fclose(fp11); fclose(fp12);
			fclose(fp13); fclose(fp14); fclose(fp15); fclose(fp16);
			fclose(fp17); fclose(fp18); fclose(fp19); fclose(fp20); 
		}

		for(i=0;i<beats;i++){
			for(j=0;j<MEDIA_SITE;j++){
				free(vmax[i][j]); free(dvdtmax[i][j]); free(apd[i][j]); free(toneapd[i][j]); free(ttwoapd[i][j]);
				free(rmbp[i][j]); free(nair[i][j]); free(cair[i][j]); free(kir[i][j]); free(caimax[i][j]);
				free(inamax[i][j]);
			}
			free(vmax[i]); free(dvdtmax[i]); free(apd[i]); free(toneapd[i]); free(ttwoapd[i]);
			free(rmbp[i]); free(nair[i]); free(cair[i]); free(kir[i]); free(caimax[i]); free(inamax[i]);
		}
		free(vmax); free(dvdtmax); free(apd); free(toneapd); free(ttwoapd);
		free(rmbp); free(nair); free(cair); free(kir); free(caimax); free(inamax);

		free(total_ina);
		for(i=0;i<beats;i++){
			free(total_inamax[i]);
		}
		free(total_inamax);

		close_mem();

	}
}

