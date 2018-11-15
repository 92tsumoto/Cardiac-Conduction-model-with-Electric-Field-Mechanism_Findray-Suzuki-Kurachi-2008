#include <string.h>
#include "syspara.h"

void vm_data(FILE *fp2, double time)
{

	int i,j;

	fprintf(fp2,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		for (j=0;j<MEDIA_PATCH;j++){
			fprintf(fp2,"%lf ",var.u[0][i][j]);
		}
	}
	fprintf(fp2,"\n");

} 

void intra_v_data(FILE *fp10, double time)
{

	int i,j;

	fprintf(fp10,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		for (j=0;j<MEDIA_PATCH;j++){
			fprintf(fp10,"%lf ",var.inner_v[i][j]);
		}
	}
	fprintf(fp10,"\n");

} 

//void ina_data(FILE *fp3, FILE *fp17, FILE *fp18, double time)
void ina_data(FILE *fp3, FILE *fp24, double time)
//void ina_data(FILE *fp3, double time)
{

	int i,j;
	
	fprintf(fp3,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		for (j=0;j<MEDIA_PATCH;j++){
			fprintf(fp3,"%lf ",var.ina[i][j]);
		}
	}
	fprintf(fp3,"\n");
	
	fprintf(fp24,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		fprintf(fp24,"%lf ",var.total_ina[i]);
	}
	fprintf(fp24,"\n");
/*	 	
	fprintf(fp17,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		for (j=0;j<MEDIA_PATCH;j++){
			fprintf(fp17,"%lf ",var.inak[i][j]);
		}
	}
	fprintf(fp17,"\n");
	 	
	fprintf(fp18,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		for (j=0;j<MEDIA_PATCH;j++){
			fprintf(fp18,"%lf ",var.inaca[i][j]);
		}
	}
	fprintf(fp18,"\n");
*/
}

//void na_active_data(FILE *fp20, FILE *fp21, FILE *fp22, double time)
void na_active_data(FILE *fp21, FILE *fp22, double time)
{

	int i,j;

/*	fprintf(fp20,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		for (j=0;j<MEDIA_PATCH;j++){
			fprintf(fp20,"%lf ",var.u[1][i][j]);
		}
	}
	fprintf(fp20,"\n");
*/
	fprintf(fp21,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		for (j=0;j<MEDIA_PATCH;j++){
			fprintf(fp21,"%lf ",var.u[2][i][j]);
		}
	}
	fprintf(fp21,"\n");
	
	fprintf(fp22,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		for (j=0;j<MEDIA_PATCH;j++){
			fprintf(fp22,"%lf ",var.u[3][i][j]);
		}
	}
	fprintf(fp22,"\n");

} 

void cai_data(FILE *fp20, double time)
{

	int i,j;

	fprintf(fp20,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		for (j=0;j<MEDIA_PATCH;j++){
			fprintf(fp20,"%lf ",var.u[18][i][j]);
		}
	}
	fprintf(fp20,"\n");

} 

//void ik_data(FILE *fp13, FILE *fp14, FILE *fp15, FILE *fp16, double time)
void ik_data(FILE *fp13, double time)
{

	int i,j;
	
	fprintf(fp13,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		for (j=0;j<MEDIA_PATCH;j++){
			fprintf(fp13,"%lf ",var.ito[i][j]);
		}
	}
	fprintf(fp13,"\n");
/*
	fprintf(fp14,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		for (j=0;j<MEDIA_PATCH;j++){
			fprintf(fp14,"%lf ",var.iki[i][j]);
		}
	}
	fprintf(fp14,"\n");
	 	
	fprintf(fp15,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		for (j=0;j<MEDIA_PATCH;j++){
			fprintf(fp15,"%lf ",var.iks[i][j]);
		}
	}
	fprintf(fp15,"\n");
	 	
	fprintf(fp16,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		for (j=0;j<MEDIA_PATCH;j++){
			fprintf(fp16,"%lf ",var.ikr[i][j]);
		}
	}
	fprintf(fp16,"\n");
*/	 	
}

void intra_i_data(FILE *fp6, double time)
{

	int i,j;

	fprintf(fp6,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		for (j=0;j<MEDIA_PATCH;j++){
			fprintf(fp6,"%lf ",var.inner_current[i][j]);
		}
	}
	fprintf(fp6,"\n");
 	
}

void cleft_v_data(FILE *fp7, FILE *fp8, double time)
{

	int i;

	fprintf(fp7,"%lf ",time);
	for(i=0;i<MEDIA_SITE-1;i++){
		fprintf(fp7,"%lf ",var.cleft_potential[i]);
	}
	fprintf(fp7,"\n");
	 	
	fprintf(fp8,"%lf ",time);
	for(i=0;i<MEDIA_SITE-1;i++){
		fprintf(fp8,"%lf %lf ",var.cleft_axial[i][0],var.cleft_axial[i][1]);
	}
	fprintf(fp8,"\n");
 	
}


void ilca_data(FILE *fp11, double time)
{

	int i,j;
	
	fprintf(fp11,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		for (j=0;j<MEDIA_PATCH;j++){
			fprintf(fp11,"%lf ",var.ilca[i][j]);
		}
	}
	fprintf(fp11,"\n");
	 	
}

void itotal_data(FILE *fp12, double time)
{

	int i,j;
	
	fprintf(fp12,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		for (j=0;j<MEDIA_PATCH;j++){
			fprintf(fp12,"%lf ",var.INa_total[i][j]+var.IK_total[i][j]+var.ICa_total[i][j]);
		}
	}
	fprintf(fp12,"\n");
	 	
}

void dvdt_data(FILE *fp19, double time)
{

	int i,j;
   
	fprintf(fp19,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		for (j=0;j<MEDIA_PATCH;j++){
			fprintf(fp19,"%lf ",var.dvdt[i][j]);
		}
	}
	fprintf(fp19,"\n");
      
}

