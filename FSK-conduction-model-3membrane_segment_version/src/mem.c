#include "syspara.h"

typedef double Number;
typedef long long Lint;

void mem()
{
	int i,k,m,l,z;

// initialized memory for state variables and arguments for Eular method 
        var.u=(double***)malloc(sizeof(double**)*NN);
        var.k1=(double***)malloc(sizeof(double**)*NN);
        if(var.u==NULL || var.k1==NULL) exit(1);
        for(i=0;i<NN;i++){
                var.u[i] = (double**)malloc(sizeof(double*)*MEDIA_SITE);
                var.k1[i] = (double**)malloc(sizeof(double*)*MEDIA_SITE);
                if(var.u[i]==NULL || var.k1[i]==NULL) exit(1);
                for(k=0;k<MEDIA_SITE;k++){
                        var.u[i][k] = (double*)malloc(sizeof(double)*MEDIA_PATCH);
                        var.k1[i][k] = (double*)malloc(sizeof(double)*MEDIA_PATCH);
                        if(var.u[i][k]==NULL || var.k1[i][k]==NULL) exit(1);
                }
        }

// initialized tablization memorys for Exp functions

	var.Tam=(Number *)calloc(VNMAX,sizeof(Number));
	var.Tbm=(Number *)calloc(VNMAX,sizeof(Number));
	var.Tah=(Number *)calloc(VNMAX,sizeof(Number));
	var.Tbh=(Number *)calloc(VNMAX,sizeof(Number));
	var.Taj=(Number *)calloc(VNMAX,sizeof(Number));
	var.Tbj=(Number *)calloc(VNMAX,sizeof(Number));
	var.Txs1ss=(Number *)calloc(VNMAX,sizeof(Number));
	var.Ttauxs1=(Number *)calloc(VNMAX,sizeof(Number));
	var.Txrss=(Number *)calloc(VNMAX,sizeof(Number));
	var.Ttauxr=(Number *)calloc(VNMAX,sizeof(Number));
	var.Tr=(Number *)calloc(VNMAX,sizeof(Number));
	var.T1Ki=(Number *)calloc(VNMAX,sizeof(Number));
	var.T2Ki=(Number *)calloc(VNMAX,sizeof(Number));
	var.T3Ki=(Number *)calloc(VNMAX,sizeof(Number));
	var.T4Ki=(Number *)calloc(VNMAX,sizeof(Number));
	var.Tkp=(Number *)calloc(VNMAX,sizeof(Number));
	var.Trvdv=(Number *)calloc(VNMAX,sizeof(Number));
	var.Tazdv=(Number *)calloc(VNMAX,sizeof(Number));
	var.Tbzdv=(Number *)calloc(VNMAX,sizeof(Number));
	var.Taydv=(Number *)calloc(VNMAX,sizeof(Number));
	var.Tbydv=(Number *)calloc(VNMAX,sizeof(Number));
	var.Tbss=(Number *)calloc(VNMAX,sizeof(Number));
	var.Ttaub=(Number *)calloc(VNMAX,sizeof(Number));
	var.Tgss=(Number *)calloc(VNMAX,sizeof(Number));
	var.Ttaug=(Number *)calloc(VNMAX,sizeof(Number));
	var.TexpCa=(Number *)calloc(VNMAX,sizeof(Number));
	var.Tdss=(Number *)calloc(VNMAX,sizeof(Number));
	var.Ttaud=(Number *)calloc(VNMAX,sizeof(Number));
	var.Tfss1=(Number *)calloc(VNMAX,sizeof(Number));
	var.Tfss2=(Number *)calloc(VNMAX,sizeof(Number));
	var.Ttauf1=(Number *)calloc(VNMAX,sizeof(Number));
	var.Ttauf2=(Number *)calloc(VNMAX,sizeof(Number));
	var.Texp0=(Number *)calloc(VNMAX,sizeof(Number));
	var.Texp1=(Number *)calloc(VNMAX,sizeof(Number));
	var.Texp2=(Number *)calloc(VNMAX,sizeof(Number));
	var.Texp3=(Number *)calloc(VNMAX,sizeof(Number));
	if(var.Tam==NULL || var.Tah==NULL || var.Taj==NULL || var.Tbm==NULL || var.Tbh==NULL || var.Tbj==NULL
		|| var.Txs1ss==NULL || var.Ttauxs1==NULL || var.Txrss==NULL || var.Ttauxr==NULL || var.Tr==NULL
		|| var.Tkp==NULL || var.Tbss==NULL || var.Ttaub==NULL || var.Tgss==NULL || var.Ttaug==NULL
		|| var.TexpCa==NULL || var.Tdss==NULL || var.Ttaud==NULL || var.Tfss1==NULL || var.Tfss2==NULL 
		|| var.Ttauf1==NULL || var.Ttauf2==NULL || var.Trvdv==NULL || var.Tazdv==NULL || var.Tbzdv==NULL
		|| var.Taydv==NULL || var.Tbydv==NULL || var.Texp0==NULL || var.Texp1==NULL || var.Texp2==NULL 
		|| var.Texp3==NULL || var.T1Ki==NULL || var.T2Ki==NULL || var.T3Ki==NULL || var.T4Ki==NULL ) exit(1);
	var.TEna=(Number *)calloc(NAIMAX,sizeof(Number));
	var.TEk=(Number *)calloc(KIMAX,sizeof(Number));
	var.Tnak=(Number *)calloc(NAIMAX,sizeof(Number));
	var.T5Ki=(Number *)calloc(KIMAX,sizeof(Number));
	var.T6Ki=(Number *)calloc(KIMAX,sizeof(Number));
	var.T7Ki=(Number *)calloc(KIMAX,sizeof(Number));
	var.T8Ki=(Number *)calloc(KIMAX,sizeof(Number));
	if(var.TEna==NULL || var.TEk==NULL || var.Tnak==NULL 
		|| var.T5Ki==NULL || var.T6Ki==NULL || var.T7Ki==NULL || var.T8Ki==NULL ) exit(1);


}
		
/* -------------------------------------------------------------------- */
/* .. Termination and release of memory. */
/* -------------------------------------------------------------------- */
void close_mem()
{

	int i,j;

		free(var.b);
		free(var.solv);
		free(var.AT);
		free(var.ia);
		free(var.ja);
		free(var.total_ina);

		for(i=0;i<NN;i++){
			for(j=0;j<MEDIA_SITE;j++){
				free(var.u[i][j]);
				free(var.k1[i][j]);
			}
			free(var.u[i]);
			free(var.k1[i]);
		}
		free(var.u); free(var.k1);

		free(var.Tam);free(var.Tah);free(var.Taj);free(var.Tbm);free(var.Tbh);free(var.Tbj);
		free(var.Txrss);free(var.Ttauxr);free(var.Tr);free(var.Txs1ss);free(var.Ttauxs1);
		free(var.T1Ki);free(var.T2Ki);free(var.T3Ki);free(var.T4Ki);
		free(var.Tkp);free(var.Trvdv);free(var.Tazdv);free(var.Tbzdv);free(var.Taydv);free(var.Tbydv);
		free(var.Tbss);free(var.Ttaub);free(var.Tgss);free(var.Ttaug);
		free(var.TexpCa);free(var.Tdss);free(var.Ttaud);free(var.Tfss1);free(var.Ttauf1);free(var.Tfss2);free(var.Ttauf2);
		free(var.Texp0);free(var.Texp1);free(var.Texp2);free(var.Texp3);
		free(var.TEna);free(var.TEk);free(var.Tnak);free(var.Tc);free(var.Tol);
		free(var.T5Ki);free(var.T6Ki);free(var.T7Ki);free(var.T8Ki);

}
