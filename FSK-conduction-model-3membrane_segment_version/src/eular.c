
#include "syspara.h"

void eular(int n, double h, double t)
{
    
	int i,j,k;  
	double ***xtemp;
        
	xtemp=(double***)malloc(sizeof(double**)*NN);
	if(xtemp==NULL) exit(1);
	for(i=0;i<NN;i++){
		xtemp[i]=(double**)malloc(sizeof(double*)*MEDIA_SITE);
		if(xtemp[i]==NULL) exit(1);
		for(j=0;j<MEDIA_SITE;j++){
			xtemp[i][j]=(double*)malloc(sizeof(double)*MEDIA_PATCH);
			if(xtemp[i][j]==NULL) exit(1);
		}
	}

	 // Calculation of input current into each unit.

	linear_vec();
	linear_solve();

	// Store of cleft potential
		for(i=0;i<MEDIA_SITE-1;i++){
			var.cleft_potential[i] = var.Rj*(var.solv[(MEDIA_PATCH-1)+MEDIA_PATCH*i] + var.solv[MEDIA_PATCH+MEDIA_PATCH*i]);
			var.cleft_axial[i][0] = 0.5*var.Rd*var.solv[(MEDIA_PATCH-1)+MEDIA_PATCH*i];
			var.cleft_axial[i][1] = 0.5*var.Rd*var.solv[MEDIA_PATCH+MEDIA_PATCH*i];
		}
	// Determination of membrane potential in JM from intracellular - extracellular potential
		for(i=0;i<MEDIA_SITE;i++){
			for(j=0;j<MEDIA_PATCH;j++){
				if(j==0 && i>0){
					var.inner_v[i][j] = var.u[0][i][j] + 0.5*var.Rd*var.solv[MEDIA_PATCH*i] + var.Rj*(var.solv[-1+MEDIA_PATCH*i]+var.solv[MEDIA_PATCH*i]);
				} else if(j==MEDIA_PATCH-1 && i<MEDIA_SITE-1){
					var.inner_v[i][j] = var.u[0][i][j] + 0.5*var.Rd*var.solv[j+MEDIA_PATCH*i] + var.Rj*(var.solv[j+MEDIA_PATCH*i]+var.solv[MEDIA_PATCH*(i+1)]);
				} else{
					var.inner_v[i][j] = var.u[0][i][j];
				}
			}
		}

	// Calculation of Ion channel currents
	for(i=0;i<MEDIA_SITE;i++){
		for(j=0;j<MEDIA_PATCH;j++){
				function(i,j,t);
			for(k=0;k<n;k++){
				xtemp[k][i][j] = var.k1[k][i][j];
			}	
		}
	}

	for(i=0;i<MEDIA_SITE;i++){
		for(j=0;j<MEDIA_PATCH;j++){
			for(k=0;k<n;k++){
				if(k == 0){
					if(j==0 || j==MEDIA_PATCH-1){
						var.inner_current[i][j] = 1000*h/(var.RGC*var.s2[j])*var.solv[j+MEDIA_PATCH*i];
						var.u[k][i][j] = var.u[k][i][j] + h*xtemp[k][i][j] + var.inner_current[i][j];
						if(j==0){
							var.disk_membrane[i][0] = var.u[0][i][0];
						} else if(j==MEDIA_PATCH-1){
							var.disk_membrane[i][1] = var.u[0][i][MEDIA_PATCH-1];
						}
						//printf("ion_current[%d][%d]=%e, in_current[%d][%d]=%e\n",i,j,xtemp[0][i][j],i,j,var.inner_current[i][j]);
					} else {
						var.inner_current[i][j] = 1000*h/(var.RGC*var.s2[j])*var.solv[j+MEDIA_PATCH*i];
						var.u[k][i][j] = var.u[k][i][j] + h*xtemp[k][i][j] + var.inner_current[i][j];
						//printf("ion_current[%d][%d]=%e, in_current[%d][%d]=%e\n",i,j,xtemp[0][i][j],i,j,var.inner_current[i][j]);
					}
				} else {
					var.u[k][i][j] = var.u[k][i][j] + h*xtemp[k][i][j];
				}
			}
		}
	}

        for(i=0;i<NN;i++){
                for(j=0;j<MEDIA_SITE;j++){
                        free(xtemp[i][j]);
                }
                free(xtemp[i]);
        }
        free(xtemp);

	for(i=0;i<MEDIA_SITE;i++){
		for(j=0;j<MEDIA_PATCH;j++){
			comp_iconcent2(i,j);
		}
	}
	
	// Determination of membrane potential in JM from intracellular - extracellular potential
/*		for(i=0;i<MEDIA_SITE-1;i++){
			for(j=0;j<MEDIA_PATCH;j++){
				if(j==0 && i>0){
					var.inner_v[i][j] = var.u[0][i][j] + 0.5*var.Rd[i]*var.solv[MEDIA_PATCH*i] + var.Rj[i]*(var.solv[-1+MEDIA_PATCH*i]+var.solv[MEDIA_PATCH*i]);
				} else if(j==MEDIA_PATCH-1 && i<MEDIA_SITE-1){
					var.inner_v[i][j] = var.u[0][i][j] + 0.5*var.Rd[i]*var.solv[j+MEDIA_PATCH*i] + var.Rj[i]*(var.solv[j+MEDIA_PATCH*i]+var.solv[MEDIA_PATCH*(i+1)]);
				} else{
					var.inner_v[i][j] = var.u[0][i][j];
				}
			}
		}
*/
}
