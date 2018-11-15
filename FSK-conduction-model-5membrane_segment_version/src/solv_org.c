/* solving the matrix equation A*x=b using LAPACK */
 
#include "syspara.h"

void linear_coeff()
{
	unsigned int i,j ,k; 
	int count,m,p,cflag;

/*	for (i=0; i<MAT_SIZE; i++){
  		for(j=0; j<MAT_SIZE; j++){
			var.AT[i+j*MAT_SIZE]=0.0;
		}
	}*/

// N = 1
			var.AT[0+0*MAT_SIZE]=0.5*var.Rix;
			var.AT[0+1*MAT_SIZE]=-var.p1;
			var.AT[0+2*MAT_SIZE]=0.0;
			var.AT[0+3*MAT_SIZE]=-var.p2;
			var.AT[0+4*MAT_SIZE]=-var.p3;
			var.AT[0+(MX*MEDIA_PATCH+2)*MAT_SIZE]=-var.p4;

			var.AT[1+0*MAT_SIZE]=0.0;
			var.AT[1+1*MAT_SIZE]=var.pa-var.p1;
			var.AT[1+2*MAT_SIZE]=0.0;
			var.AT[1+3*MAT_SIZE]=-var.p2;
			var.AT[1+4*MAT_SIZE]=var.pb-var.p3;
			var.AT[1+(MX*MEDIA_PATCH+2)*MAT_SIZE]=-var.p4;

			var.AT[2+0*MAT_SIZE]=0.0;
			var.AT[2+1*MAT_SIZE]=-var.p1;
			var.AT[2+2*MAT_SIZE]=0.5*var.Riy;
			var.AT[2+3*MAT_SIZE]=-var.p2;
			var.AT[2+4*MAT_SIZE]=-var.p3;
			var.AT[2+(MX*MEDIA_PATCH+2)*MAT_SIZE]=-var.p4;


			var.AT[3+0*MAT_SIZE]=0.0;
			var.AT[3+1*MAT_SIZE]=-var.p1;
			var.AT[3+2*MAT_SIZE]=0.0;
			var.AT[3+3*MAT_SIZE]=var.pe-var.p2;
			var.AT[3+4*MAT_SIZE]=-var.p3;
			var.AT[3+(MX*MEDIA_PATCH+2)*MAT_SIZE]=var.pf-var.p4;

// N = MX (N = 3)

			var.AT[(MX-1)*MEDIA_PATCH+((MX-2)*MEDIA_PATCH+1)*MAT_SIZE]=var.pb-var.p3;
			var.AT[(MX-1)*MEDIA_PATCH+((MX-1)*MEDIA_PATCH+0)*MAT_SIZE]=var.pa-var.p1;
			var.AT[(MX-1)*MEDIA_PATCH+((MX-1)*MEDIA_PATCH+1)*MAT_SIZE]=0.0;
			var.AT[(MX-1)*MEDIA_PATCH+((MX-1)*MEDIA_PATCH+2)*MAT_SIZE]=0.0;
			var.AT[(MX-1)*MEDIA_PATCH+((MX-1)*MEDIA_PATCH+3)*MAT_SIZE]=-var.p2;
			var.AT[(MX-1)*MEDIA_PATCH+((2*MX-1)*MEDIA_PATCH+2)*MAT_SIZE]=-var.p4;

			var.AT[(MX-1)*MEDIA_PATCH+1+((MX-2)*MEDIA_PATCH+1)*MAT_SIZE]=-var.p3;
			var.AT[(MX-1)*MEDIA_PATCH+1+((MX-1)*MEDIA_PATCH+0)*MAT_SIZE]=-var.p1;
			var.AT[(MX-1)*MEDIA_PATCH+1+((MX-1)*MEDIA_PATCH+1)*MAT_SIZE]=0.5*var.Rix;
			var.AT[(MX-1)*MEDIA_PATCH+1+((MX-1)*MEDIA_PATCH+2)*MAT_SIZE]=0.0;
			var.AT[(MX-1)*MEDIA_PATCH+1+((MX-1)*MEDIA_PATCH+3)*MAT_SIZE]=-var.p2;
			var.AT[(MX-1)*MEDIA_PATCH+1+((2*MX-1)*MEDIA_PATCH+2)*MAT_SIZE]=-var.p4;

			var.AT[(MX-1)*MEDIA_PATCH+2+((MX-2)*MEDIA_PATCH+1)*MAT_SIZE]=-var.p3;
			var.AT[(MX-1)*MEDIA_PATCH+2+((MX-1)*MEDIA_PATCH+0)*MAT_SIZE]=-var.p1;
			var.AT[(MX-1)*MEDIA_PATCH+2+((MX-1)*MEDIA_PATCH+1)*MAT_SIZE]=0.0;
			var.AT[(MX-1)*MEDIA_PATCH+2+((MX-1)*MEDIA_PATCH+2)*MAT_SIZE]=0.5*var.Riy;
			var.AT[(MX-1)*MEDIA_PATCH+2+((MX-1)*MEDIA_PATCH+3)*MAT_SIZE]=-var.p2;
			var.AT[(MX-1)*MEDIA_PATCH+2+((2*MX-1)*MEDIA_PATCH+2)*MAT_SIZE]=-var.p4;

			var.AT[(MX-1)*MEDIA_PATCH+3+((MX-2)*MEDIA_PATCH+1)*MAT_SIZE]=-var.p3;
			var.AT[(MX-1)*MEDIA_PATCH+3+((MX-1)*MEDIA_PATCH+0)*MAT_SIZE]=-var.p1;
			var.AT[(MX-1)*MEDIA_PATCH+3+((MX-1)*MEDIA_PATCH+1)*MAT_SIZE]=0.0;
			var.AT[(MX-1)*MEDIA_PATCH+3+((MX-1)*MEDIA_PATCH+2)*MAT_SIZE]=0.0;
			var.AT[(MX-1)*MEDIA_PATCH+3+((MX-1)*MEDIA_PATCH+3)*MAT_SIZE]=var.pe-var.p2;
			var.AT[(MX-1)*MEDIA_PATCH+3+((2*MX-1)*MEDIA_PATCH+2)*MAT_SIZE]=var.pf-var.p4;

// N=(MY-1)*MX+1 (N = 7)

			var.AT[MX*(MY-1)*MEDIA_PATCH+0+((MX*(MY-2))*MEDIA_PATCH+3)*MAT_SIZE]=-var.p4;
			var.AT[MX*(MY-1)*MEDIA_PATCH+0+((MX*(MY-1))*MEDIA_PATCH+0)*MAT_SIZE]=0.5*var.Rix;
			var.AT[MX*(MY-1)*MEDIA_PATCH+0+((MX*(MY-1))*MEDIA_PATCH+1)*MAT_SIZE]=-var.p1;
			var.AT[MX*(MY-1)*MEDIA_PATCH+0+((MX*(MY-1))*MEDIA_PATCH+2)*MAT_SIZE]=-var.p2;
			var.AT[MX*(MY-1)*MEDIA_PATCH+0+((MX*(MY-1))*MEDIA_PATCH+3)*MAT_SIZE]=0.0;
			var.AT[MX*(MY-1)*MEDIA_PATCH+0+((MX*(MY-1))*MEDIA_PATCH+4)*MAT_SIZE]=-var.p3;

			var.AT[MX*(MY-1)*MEDIA_PATCH+1+((MX*(MY-2))*MEDIA_PATCH+3)*MAT_SIZE]=-var.p4;
			var.AT[MX*(MY-1)*MEDIA_PATCH+1+((MX*(MY-1))*MEDIA_PATCH+0)*MAT_SIZE]=0.0;
			var.AT[MX*(MY-1)*MEDIA_PATCH+1+((MX*(MY-1))*MEDIA_PATCH+1)*MAT_SIZE]=var.pa-var.p1;
			var.AT[MX*(MY-1)*MEDIA_PATCH+1+((MX*(MY-1))*MEDIA_PATCH+2)*MAT_SIZE]=-var.p2;
			var.AT[MX*(MY-1)*MEDIA_PATCH+1+((MX*(MY-1))*MEDIA_PATCH+3)*MAT_SIZE]=0.0;
			var.AT[MX*(MY-1)*MEDIA_PATCH+1+((MX*(MY-1))*MEDIA_PATCH+4)*MAT_SIZE]=var.pb-var.p3;

			var.AT[MX*(MY-1)*MEDIA_PATCH+2+((MX*(MY-2))*MEDIA_PATCH+3)*MAT_SIZE]=var.pf-var.p4;
			var.AT[MX*(MY-1)*MEDIA_PATCH+2+((MX*(MY-1))*MEDIA_PATCH+0)*MAT_SIZE]=0.0;
			var.AT[MX*(MY-1)*MEDIA_PATCH+2+((MX*(MY-1))*MEDIA_PATCH+1)*MAT_SIZE]=-var.p1;
			var.AT[MX*(MY-1)*MEDIA_PATCH+2+((MX*(MY-1))*MEDIA_PATCH+2)*MAT_SIZE]=var.pe-var.p2;
			var.AT[MX*(MY-1)*MEDIA_PATCH+2+((MX*(MY-1))*MEDIA_PATCH+3)*MAT_SIZE]=0.0;
			var.AT[MX*(MY-1)*MEDIA_PATCH+2+((MX*(MY-1))*MEDIA_PATCH+4)*MAT_SIZE]=-var.p3;

			var.AT[MX*(MY-1)*MEDIA_PATCH+3+((MX*(MY-2))*MEDIA_PATCH+3)*MAT_SIZE]=-var.p4;
			var.AT[MX*(MY-1)*MEDIA_PATCH+3+((MX*(MY-1))*MEDIA_PATCH+0)*MAT_SIZE]=0.0;
			var.AT[MX*(MY-1)*MEDIA_PATCH+3+((MX*(MY-1))*MEDIA_PATCH+1)*MAT_SIZE]=-var.p1;
			var.AT[MX*(MY-1)*MEDIA_PATCH+3+((MX*(MY-1))*MEDIA_PATCH+2)*MAT_SIZE]=-var.p2;
			var.AT[MX*(MY-1)*MEDIA_PATCH+3+((MX*(MY-1))*MEDIA_PATCH+3)*MAT_SIZE]=0.5*var.Riy;
			var.AT[MX*(MY-1)*MEDIA_PATCH+3+((MX*(MY-1))*MEDIA_PATCH+4)*MAT_SIZE]=-var.p3;

// N = 9

			var.AT[(MX*MY-1)*MEDIA_PATCH+0+((MX*(MY-1)-1)*MEDIA_PATCH+3)*MAT_SIZE]=-var.p4;
			var.AT[(MX*MY-1)*MEDIA_PATCH+0+((MX*MY-2)*MEDIA_PATCH+1)*MAT_SIZE]=var.pb-var.p3;
			var.AT[(MX*MY-1)*MEDIA_PATCH+0+((MX*MY-1)*MEDIA_PATCH+0)*MAT_SIZE]=var.pa-var.p1;
			var.AT[(MX*MY-1)*MEDIA_PATCH+0+((MX*MY-1)*MEDIA_PATCH+1)*MAT_SIZE]=0.0;
			var.AT[(MX*MY-1)*MEDIA_PATCH+0+((MX*MY-1)*MEDIA_PATCH+2)*MAT_SIZE]=-var.p2;
			var.AT[(MX*MY-1)*MEDIA_PATCH+0+((MX*MY-1)*MEDIA_PATCH+3)*MAT_SIZE]=0.0;

			var.AT[(MX*MY-1)*MEDIA_PATCH+1+((MX*(MY-1)-1)*MEDIA_PATCH+3)*MAT_SIZE]=-var.p4;
			var.AT[(MX*MY-1)*MEDIA_PATCH+1+((MX*MY-2)*MEDIA_PATCH+1)*MAT_SIZE]=-var.p3;
			var.AT[(MX*MY-1)*MEDIA_PATCH+1+((MX*MY-1)*MEDIA_PATCH+0)*MAT_SIZE]=-var.p1;
			var.AT[(MX*MY-1)*MEDIA_PATCH+1+((MX*MY-1)*MEDIA_PATCH+1)*MAT_SIZE]=0.5*var.Rix;
			var.AT[(MX*MY-1)*MEDIA_PATCH+1+((MX*MY-1)*MEDIA_PATCH+2)*MAT_SIZE]=-var.p2;
			var.AT[(MX*MY-1)*MEDIA_PATCH+1+((MX*MY-1)*MEDIA_PATCH+3)*MAT_SIZE]=0.0;

			var.AT[(MX*MY-1)*MEDIA_PATCH+2+((MX*(MY-1)-1)*MEDIA_PATCH+3)*MAT_SIZE]=var.pf-var.p4;
			var.AT[(MX*MY-1)*MEDIA_PATCH+2+((MX*MY-2)*MEDIA_PATCH+1)*MAT_SIZE]=-var.p3;
			var.AT[(MX*MY-1)*MEDIA_PATCH+2+((MX*MY-1)*MEDIA_PATCH+0)*MAT_SIZE]=-var.p1;
			var.AT[(MX*MY-1)*MEDIA_PATCH+2+((MX*MY-1)*MEDIA_PATCH+1)*MAT_SIZE]=0.0;
			var.AT[(MX*MY-1)*MEDIA_PATCH+2+((MX*MY-1)*MEDIA_PATCH+2)*MAT_SIZE]=var.pe-var.p2;
			var.AT[(MX*MY-1)*MEDIA_PATCH+2+((MX*MY-1)*MEDIA_PATCH+3)*MAT_SIZE]=0.0;
			
			var.AT[(MX*MY-1)*MEDIA_PATCH+3+((MX*(MY-1)-1)*MEDIA_PATCH+3)*MAT_SIZE]=-var.p4;
			var.AT[(MX*MY-1)*MEDIA_PATCH+3+((MX*MY-2)*MEDIA_PATCH+1)*MAT_SIZE]=-var.p3;
			var.AT[(MX*MY-1)*MEDIA_PATCH+3+((MX*MY-1)*MEDIA_PATCH+0)*MAT_SIZE]=-var.p1;
			var.AT[(MX*MY-1)*MEDIA_PATCH+3+((MX*MY-1)*MEDIA_PATCH+1)*MAT_SIZE]=0.0;
			var.AT[(MX*MY-1)*MEDIA_PATCH+3+((MX*MY-1)*MEDIA_PATCH+2)*MAT_SIZE]=-var.p2;
			var.AT[(MX*MY-1)*MEDIA_PATCH+3+((MX*MY-1)*MEDIA_PATCH+3)*MAT_SIZE]=0.5*var.Riy;

  		for(k=1;k<(MX-1);k++){ // k = 2,...,MX-1

			var.AT[k*MEDIA_PATCH+((k-1)*MEDIA_PATCH+1)*MAT_SIZE]=var.pb-var.p3;
			var.AT[k*MEDIA_PATCH+(k*MEDIA_PATCH+0)*MAT_SIZE]=var.pa-var.p1;
			var.AT[k*MEDIA_PATCH+(k*MEDIA_PATCH+1)*MAT_SIZE]=-var.p1;
			var.AT[k*MEDIA_PATCH+(k*MEDIA_PATCH+2)*MAT_SIZE]=0.0;
			var.AT[k*MEDIA_PATCH+(k*MEDIA_PATCH+3)*MAT_SIZE]=-var.p2;
			var.AT[k*MEDIA_PATCH+(k*MEDIA_PATCH+4)*MAT_SIZE]=-var.p3;
			var.AT[k*MEDIA_PATCH+((k+MX)*MEDIA_PATCH+2)*MAT_SIZE]=-var.p4;

			var.AT[k*MEDIA_PATCH+1+((k-1)*MEDIA_PATCH+1)*MAT_SIZE]=-var.p3;
			var.AT[k*MEDIA_PATCH+1+(k*MEDIA_PATCH+0)*MAT_SIZE]=-var.p1;
			var.AT[k*MEDIA_PATCH+1+(k*MEDIA_PATCH+1)*MAT_SIZE]=var.pa-var.p1;
			var.AT[k*MEDIA_PATCH+1+(k*MEDIA_PATCH+2)*MAT_SIZE]=0.0;
			var.AT[k*MEDIA_PATCH+1+(k*MEDIA_PATCH+3)*MAT_SIZE]=-var.p2;
			var.AT[k*MEDIA_PATCH+1+(k*MEDIA_PATCH+4)*MAT_SIZE]=var.pb-var.p3;
			var.AT[k*MEDIA_PATCH+1+((k+MX)*MEDIA_PATCH+2)*MAT_SIZE]=-var.p4;

			var.AT[k*MEDIA_PATCH+2+((k-1)*MEDIA_PATCH+1)*MAT_SIZE]=-var.p3;
			var.AT[k*MEDIA_PATCH+2+(k*MEDIA_PATCH+0)*MAT_SIZE]=-var.p1;
			var.AT[k*MEDIA_PATCH+2+(k*MEDIA_PATCH+1)*MAT_SIZE]=-var.p1;
			var.AT[k*MEDIA_PATCH+2+(k*MEDIA_PATCH+2)*MAT_SIZE]=0.5*var.Riy;
			var.AT[k*MEDIA_PATCH+2+(k*MEDIA_PATCH+3)*MAT_SIZE]=-var.p2;
			var.AT[k*MEDIA_PATCH+2+(k*MEDIA_PATCH+4)*MAT_SIZE]=-var.p3;
			var.AT[k*MEDIA_PATCH+2+((k+MX)*MEDIA_PATCH+2)*MAT_SIZE]=-var.p4;

			var.AT[k*MEDIA_PATCH+3+((k-1)*MEDIA_PATCH+1)*MAT_SIZE]=-var.p3;
			var.AT[k*MEDIA_PATCH+3+(k*MEDIA_PATCH+0)*MAT_SIZE]=-var.p1;
			var.AT[k*MEDIA_PATCH+3+(k*MEDIA_PATCH+1)*MAT_SIZE]=-var.p1;
			var.AT[k*MEDIA_PATCH+3+(k*MEDIA_PATCH+2)*MAT_SIZE]=0.0;
			var.AT[k*MEDIA_PATCH+3+(k*MEDIA_PATCH+3)*MAT_SIZE]=var.pe-var.p2;
			var.AT[k*MEDIA_PATCH+3+(k*MEDIA_PATCH+4)*MAT_SIZE]=-var.p3;
			var.AT[k*MEDIA_PATCH+3+((k+MX)*MEDIA_PATCH+2)*MAT_SIZE]=var.pf-var.p4;
		}

		for(k=MX*(MY-1)+1;k<(MX*MY-1);k++){ // k = MX(MY-1)+1 ... MX*MY-1

			var.AT[k*MEDIA_PATCH+((k-MX)*MEDIA_PATCH+3)*MAT_SIZE]=-var.p4;
			var.AT[k*MEDIA_PATCH+((k-1)*MEDIA_PATCH+1)*MAT_SIZE]=var.pb-var.p3;
			var.AT[k*MEDIA_PATCH+(k*MEDIA_PATCH+0)*MAT_SIZE]=var.pa-var.p1;
			var.AT[k*MEDIA_PATCH+(k*MEDIA_PATCH+1)*MAT_SIZE]=-var.p1;
			var.AT[k*MEDIA_PATCH+(k*MEDIA_PATCH+2)*MAT_SIZE]=-var.p2;
			var.AT[k*MEDIA_PATCH+(k*MEDIA_PATCH+3)*MAT_SIZE]=0.0;
			var.AT[k*MEDIA_PATCH+(k*MEDIA_PATCH+4)*MAT_SIZE]=-var.p3;
			
			var.AT[k*MEDIA_PATCH+1+((k-MX)*MEDIA_PATCH+3)*MAT_SIZE]=-var.p4;
			var.AT[k*MEDIA_PATCH+1+((k-1)*MEDIA_PATCH+1)*MAT_SIZE]=-var.p3;
			var.AT[k*MEDIA_PATCH+1+(k*MEDIA_PATCH+0)*MAT_SIZE]=-var.p1;
			var.AT[k*MEDIA_PATCH+1+(k*MEDIA_PATCH+1)*MAT_SIZE]=var.pa-var.p1;
			var.AT[k*MEDIA_PATCH+1+(k*MEDIA_PATCH+2)*MAT_SIZE]=-var.p2;
			var.AT[k*MEDIA_PATCH+1+(k*MEDIA_PATCH+3)*MAT_SIZE]=0.0;
			var.AT[k*MEDIA_PATCH+1+(k*MEDIA_PATCH+4)*MAT_SIZE]=var.pb-var.p3;
			
			var.AT[k*MEDIA_PATCH+2+((k-MX)*MEDIA_PATCH+3)*MAT_SIZE]=var.pf-var.p4;
			var.AT[k*MEDIA_PATCH+2+((k-1)*MEDIA_PATCH+1)*MAT_SIZE]=-var.p3;
			var.AT[k*MEDIA_PATCH+2+(k*MEDIA_PATCH+0)*MAT_SIZE]=-var.p1;
			var.AT[k*MEDIA_PATCH+2+(k*MEDIA_PATCH+1)*MAT_SIZE]=-var.p1;
			var.AT[k*MEDIA_PATCH+2+(k*MEDIA_PATCH+2)*MAT_SIZE]=var.pe-var.p2;
			var.AT[k*MEDIA_PATCH+2+(k*MEDIA_PATCH+3)*MAT_SIZE]=0.0;
			var.AT[k*MEDIA_PATCH+2+(k*MEDIA_PATCH+4)*MAT_SIZE]=-var.p3;
			
			var.AT[k*MEDIA_PATCH+3+((k-MX)*MEDIA_PATCH+3)*MAT_SIZE]=-var.p4;
			var.AT[k*MEDIA_PATCH+3+((k-1)*MEDIA_PATCH+1)*MAT_SIZE]=-var.p3;
			var.AT[k*MEDIA_PATCH+3+(k*MEDIA_PATCH+0)*MAT_SIZE]=-var.p1;
			var.AT[k*MEDIA_PATCH+3+(k*MEDIA_PATCH+1)*MAT_SIZE]=-var.p1;
			var.AT[k*MEDIA_PATCH+3+(k*MEDIA_PATCH+2)*MAT_SIZE]=-var.p2;
			var.AT[k*MEDIA_PATCH+3+(k*MEDIA_PATCH+3)*MAT_SIZE]=0.5*var.Riy;
			var.AT[k*MEDIA_PATCH+3+(k*MEDIA_PATCH+4)*MAT_SIZE]=-var.p3;
		}
	

  		for(k=1;k<(MY-1);k++){ // MY = 4,...,MX*(MY-2)

			var.AT[k*MX*MEDIA_PATCH+0+((k-1)*MX*MEDIA_PATCH+3)*MAT_SIZE]=-var.p4;
			var.AT[k*MX*MEDIA_PATCH+0+(k*MX*MEDIA_PATCH+0)*MAT_SIZE]=0.5*var.Rix;
			var.AT[k*MX*MEDIA_PATCH+0+(k*MX*MEDIA_PATCH+1)*MAT_SIZE]=-var.p1;
			var.AT[k*MX*MEDIA_PATCH+0+(k*MX*MEDIA_PATCH+2)*MAT_SIZE]=-var.p2;
			var.AT[k*MX*MEDIA_PATCH+0+(k*MX*MEDIA_PATCH+3)*MAT_SIZE]=-var.p2;
			var.AT[k*MX*MEDIA_PATCH+0+(k*MX*MEDIA_PATCH+4)*MAT_SIZE]=-var.p3;
			var.AT[k*MX*MEDIA_PATCH+0+((k+1)*MX*MEDIA_PATCH+2)*MAT_SIZE]=-var.p4;

			var.AT[k*MX*MEDIA_PATCH+1+((k-1)*MX*MEDIA_PATCH+3)*MAT_SIZE]=-var.p4;
			var.AT[k*MX*MEDIA_PATCH+1+(k*MX*MEDIA_PATCH+0)*MAT_SIZE]=0.0;
			var.AT[k*MX*MEDIA_PATCH+1+(k*MX*MEDIA_PATCH+1)*MAT_SIZE]=var.pa-var.p1;
			var.AT[k*MX*MEDIA_PATCH+1+(k*MX*MEDIA_PATCH+2)*MAT_SIZE]=-var.p2;
			var.AT[k*MX*MEDIA_PATCH+1+(k*MX*MEDIA_PATCH+3)*MAT_SIZE]=-var.p2;
			var.AT[k*MX*MEDIA_PATCH+1+(k*MX*MEDIA_PATCH+4)*MAT_SIZE]=var.pb-var.p3;
			var.AT[k*MX*MEDIA_PATCH+1+((k+1)*MX*MEDIA_PATCH+2)*MAT_SIZE]=-var.p4;

			var.AT[k*MX*MEDIA_PATCH+2+((k-1)*MX*MEDIA_PATCH+3)*MAT_SIZE]=var.pf-var.p4;
			var.AT[k*MX*MEDIA_PATCH+2+(k*MX*MEDIA_PATCH+0)*MAT_SIZE]=0.0;
			var.AT[k*MX*MEDIA_PATCH+2+(k*MX*MEDIA_PATCH+1)*MAT_SIZE]=-var.p1;
			var.AT[k*MX*MEDIA_PATCH+2+(k*MX*MEDIA_PATCH+2)*MAT_SIZE]=var.pe-var.p2;
			var.AT[k*MX*MEDIA_PATCH+2+(k*MX*MEDIA_PATCH+3)*MAT_SIZE]=-var.p2;
			var.AT[k*MX*MEDIA_PATCH+2+(k*MX*MEDIA_PATCH+4)*MAT_SIZE]=-var.p3;
			var.AT[k*MX*MEDIA_PATCH+2+((k+1)*MX*MEDIA_PATCH+2)*MAT_SIZE]=-var.p4;

			var.AT[k*MX*MEDIA_PATCH+3+((k-1)*MX*MEDIA_PATCH+3)*MAT_SIZE]=-var.p4;
			var.AT[k*MX*MEDIA_PATCH+3+(k*MX*MEDIA_PATCH+0)*MAT_SIZE]=0.0;
			var.AT[k*MX*MEDIA_PATCH+3+(k*MX*MEDIA_PATCH+1)*MAT_SIZE]=-var.p1;
			var.AT[k*MX*MEDIA_PATCH+3+(k*MX*MEDIA_PATCH+2)*MAT_SIZE]=-var.p2;
			var.AT[k*MX*MEDIA_PATCH+3+(k*MX*MEDIA_PATCH+3)*MAT_SIZE]=var.pe-var.p2;
			var.AT[k*MX*MEDIA_PATCH+3+(k*MX*MEDIA_PATCH+4)*MAT_SIZE]=-var.p3;
			var.AT[k*MX*MEDIA_PATCH+3+((k+1)*MX*MEDIA_PATCH+2)*MAT_SIZE]=var.pf-var.p4;
		}

		for (k=1; k < (MY-1); k++) { // k = (k+1)MX, (k+2)MX, ...., (MY-1)MX

			var.AT[((k+1)*MX-1)*MEDIA_PATCH+0+((k*MX-1)*MEDIA_PATCH+3)*MAT_SIZE]=-var.p4;
			var.AT[((k+1)*MX-1)*MEDIA_PATCH+0+(((k+1)*MX-2)*MEDIA_PATCH+1)*MAT_SIZE]=var.pb-var.p3;
			var.AT[((k+1)*MX-1)*MEDIA_PATCH+0+(((k+1)*MX-1)*MEDIA_PATCH+0)*MAT_SIZE]=var.pa-var.p1;
			var.AT[((k+1)*MX-1)*MEDIA_PATCH+0+(((k+1)*MX-1)*MEDIA_PATCH+1)*MAT_SIZE]=0.0;
			var.AT[((k+1)*MX-1)*MEDIA_PATCH+0+(((k+1)*MX-1)*MEDIA_PATCH+2)*MAT_SIZE]=-var.p2;
			var.AT[((k+1)*MX-1)*MEDIA_PATCH+0+(((k+1)*MX-1)*MEDIA_PATCH+3)*MAT_SIZE]=-var.p2;
			var.AT[((k+1)*MX-1)*MEDIA_PATCH+0+(((k+2)*MX-1)*MEDIA_PATCH+2)*MAT_SIZE]=-var.p4;

			var.AT[((k+1)*MX-1)*MEDIA_PATCH+1+((k*MX-1)*MEDIA_PATCH+3)*MAT_SIZE]=-var.p4;
			var.AT[((k+1)*MX-1)*MEDIA_PATCH+1+(((k+1)*MX-2)*MEDIA_PATCH+1)*MAT_SIZE]=-var.p3;
			var.AT[((k+1)*MX-1)*MEDIA_PATCH+1+(((k+1)*MX-1)*MEDIA_PATCH+0)*MAT_SIZE]=-var.p1;
			var.AT[((k+1)*MX-1)*MEDIA_PATCH+1+(((k+1)*MX-1)*MEDIA_PATCH+1)*MAT_SIZE]=0.5*var.Rix;
			var.AT[((k+1)*MX-1)*MEDIA_PATCH+1+(((k+1)*MX-1)*MEDIA_PATCH+2)*MAT_SIZE]=-var.p2;
			var.AT[((k+1)*MX-1)*MEDIA_PATCH+1+(((k+1)*MX-1)*MEDIA_PATCH+3)*MAT_SIZE]=-var.p2;
			var.AT[((k+1)*MX-1)*MEDIA_PATCH+1+(((k+2)*MX-1)*MEDIA_PATCH+2)*MAT_SIZE]=-var.p4;

			var.AT[((k+1)*MX-1)*MEDIA_PATCH+2+((k*MX-1)*MEDIA_PATCH+3)*MAT_SIZE]=var.pf-var.p4;
			var.AT[((k+1)*MX-1)*MEDIA_PATCH+2+(((k+1)*MX-2)*MEDIA_PATCH+1)*MAT_SIZE]=-var.p3;
			var.AT[((k+1)*MX-1)*MEDIA_PATCH+2+(((k+1)*MX-1)*MEDIA_PATCH+0)*MAT_SIZE]=-var.p1;
			var.AT[((k+1)*MX-1)*MEDIA_PATCH+2+(((k+1)*MX-1)*MEDIA_PATCH+1)*MAT_SIZE]=0.0;
			var.AT[((k+1)*MX-1)*MEDIA_PATCH+2+(((k+1)*MX-1)*MEDIA_PATCH+2)*MAT_SIZE]=var.pe-var.p2;
			var.AT[((k+1)*MX-1)*MEDIA_PATCH+2+(((k+1)*MX-1)*MEDIA_PATCH+3)*MAT_SIZE]=-var.p2;
			var.AT[((k+1)*MX-1)*MEDIA_PATCH+2+(((k+2)*MX-1)*MEDIA_PATCH+2)*MAT_SIZE]=-var.p4;

			var.AT[((k+1)*MX-1)*MEDIA_PATCH+3+((k*MX-1)*MEDIA_PATCH+3)*MAT_SIZE]=-var.p4;
			var.AT[((k+1)*MX-1)*MEDIA_PATCH+3+(((k+1)*MX-2)*MEDIA_PATCH+1)*MAT_SIZE]=-var.p3;
			var.AT[((k+1)*MX-1)*MEDIA_PATCH+3+(((k+1)*MX-1)*MEDIA_PATCH+0)*MAT_SIZE]=-var.p1;
			var.AT[((k+1)*MX-1)*MEDIA_PATCH+3+(((k+1)*MX-1)*MEDIA_PATCH+1)*MAT_SIZE]=0.0;
			var.AT[((k+1)*MX-1)*MEDIA_PATCH+3+(((k+1)*MX-1)*MEDIA_PATCH+2)*MAT_SIZE]=-var.p2;
			var.AT[((k+1)*MX-1)*MEDIA_PATCH+3+(((k+1)*MX-1)*MEDIA_PATCH+3)*MAT_SIZE]=var.pe-var.p2;
			var.AT[((k+1)*MX-1)*MEDIA_PATCH+3+(((k+2)*MX-1)*MEDIA_PATCH+2)*MAT_SIZE]=var.pf-var.p4;

		}


	for (i=1;i<(MY-1);i++){	 // for N = ....
  		for(j=1;j<MX-1; j++){

			var.AT[(i*MX+j)*MEDIA_PATCH+0+(((i-1)*MX+j)*MEDIA_PATCH+3)*MAT_SIZE]=-var.p4;
			var.AT[(i*MX+j)*MEDIA_PATCH+0+((i*MX+j-1)*MEDIA_PATCH+1)*MAT_SIZE]=var.pb-var.p3;
			var.AT[(i*MX+j)*MEDIA_PATCH+0+((i*MX+j)*MEDIA_PATCH+0)*MAT_SIZE]=var.pa-var.p1;
			var.AT[(i*MX+j)*MEDIA_PATCH+0+((i*MX+j)*MEDIA_PATCH+1)*MAT_SIZE]=-var.p1;
			var.AT[(i*MX+j)*MEDIA_PATCH+0+((i*MX+j)*MEDIA_PATCH+2)*MAT_SIZE]=-var.p2;
			var.AT[(i*MX+j)*MEDIA_PATCH+0+((i*MX+j)*MEDIA_PATCH+3)*MAT_SIZE]=-var.p2;
			var.AT[(i*MX+j)*MEDIA_PATCH+0+((i*MX+j)*MEDIA_PATCH+4)*MAT_SIZE]=-var.p3;
			var.AT[(i*MX+j)*MEDIA_PATCH+0+(((i+1)*MX+j)*MEDIA_PATCH+2)*MAT_SIZE]=-var.p4;
			
			var.AT[(i*MX+j)*MEDIA_PATCH+1+(((i-1)*MX+j)*MEDIA_PATCH+3)*MAT_SIZE]=-var.p4;
			var.AT[(i*MX+j)*MEDIA_PATCH+1+((i*MX+j-1)*MEDIA_PATCH+1)*MAT_SIZE]=-var.p3;
			var.AT[(i*MX+j)*MEDIA_PATCH+1+((i*MX+j)*MEDIA_PATCH+0)*MAT_SIZE]=-var.p1;
			var.AT[(i*MX+j)*MEDIA_PATCH+1+((i*MX+j)*MEDIA_PATCH+1)*MAT_SIZE]=var.pa-var.p1;
			var.AT[(i*MX+j)*MEDIA_PATCH+1+((i*MX+j)*MEDIA_PATCH+2)*MAT_SIZE]=-var.p2;
			var.AT[(i*MX+j)*MEDIA_PATCH+1+((i*MX+j)*MEDIA_PATCH+3)*MAT_SIZE]=-var.p2;
			var.AT[(i*MX+j)*MEDIA_PATCH+1+((i*MX+j)*MEDIA_PATCH+4)*MAT_SIZE]=var.pb-var.p3;
			var.AT[(i*MX+j)*MEDIA_PATCH+1+(((i+1)*MX+j)*MEDIA_PATCH+2)*MAT_SIZE]=-var.p4;
			
			var.AT[(i*MX+j)*MEDIA_PATCH+2+(((i-1)*MX+j)*MEDIA_PATCH+3)*MAT_SIZE]=var.pf-var.p4;
			var.AT[(i*MX+j)*MEDIA_PATCH+2+((i*MX+j-1)*MEDIA_PATCH+1)*MAT_SIZE]=-var.p3;
			var.AT[(i*MX+j)*MEDIA_PATCH+2+((i*MX+j)*MEDIA_PATCH+0)*MAT_SIZE]=-var.p1;
			var.AT[(i*MX+j)*MEDIA_PATCH+2+((i*MX+j)*MEDIA_PATCH+1)*MAT_SIZE]=-var.p1;
			var.AT[(i*MX+j)*MEDIA_PATCH+2+((i*MX+j)*MEDIA_PATCH+2)*MAT_SIZE]=var.pe-var.p2;
			var.AT[(i*MX+j)*MEDIA_PATCH+2+((i*MX+j)*MEDIA_PATCH+3)*MAT_SIZE]=-var.p2;
			var.AT[(i*MX+j)*MEDIA_PATCH+2+((i*MX+j)*MEDIA_PATCH+4)*MAT_SIZE]=-var.p3;
			var.AT[(i*MX+j)*MEDIA_PATCH+2+(((i+1)*MX+j)*MEDIA_PATCH+2)*MAT_SIZE]=-var.p4;
			
			var.AT[(i*MX+j)*MEDIA_PATCH+3+(((i-1)*MX+j)*MEDIA_PATCH+3)*MAT_SIZE]=-var.p4;
			var.AT[(i*MX+j)*MEDIA_PATCH+3+((i*MX+j-1)*MEDIA_PATCH+1)*MAT_SIZE]=-var.p3;
			var.AT[(i*MX+j)*MEDIA_PATCH+3+((i*MX+j)*MEDIA_PATCH+0)*MAT_SIZE]=-var.p1;
			var.AT[(i*MX+j)*MEDIA_PATCH+3+((i*MX+j)*MEDIA_PATCH+1)*MAT_SIZE]=-var.p1;
			var.AT[(i*MX+j)*MEDIA_PATCH+3+((i*MX+j)*MEDIA_PATCH+2)*MAT_SIZE]=-var.p2;
			var.AT[(i*MX+j)*MEDIA_PATCH+3+((i*MX+j)*MEDIA_PATCH+3)*MAT_SIZE]=var.pe-var.p2;
			var.AT[(i*MX+j)*MEDIA_PATCH+3+((i*MX+j)*MEDIA_PATCH+4)*MAT_SIZE]=-var.p3;
			var.AT[(i*MX+j)*MEDIA_PATCH+3+(((i+1)*MX+j)*MEDIA_PATCH+2)*MAT_SIZE]=var.pf-var.p4;

		}
	}

/*
	for (i=0; i<MAT_SIZE; i++){
		for(j=0; j<MAT_SIZE; j++){
			printf("%2.0lf,",var.AT[i+MAT_SIZE*j]);
		}
		printf("\n");
	}
*/

	// count はATの非ゼロ要素数のカウント
	count=0;
	for (i=0; i<MAT_SIZE; i++){
		for(j=0; j<MAT_SIZE; j++){
			if(var.AT[i+MAT_SIZE*j] != 0.0 ){
				count +=1;
			}
		}
	}
	printf("#all elements=%d\n",count);
	var.par_k = count;	
	//exit(0);


}


/* solving the matrix equation A*x=b using LAPACK */

void linear_vec()
{

	long int i,j,k,N;

	for(k=0;k<MY;k++){
		N=k*MX;
		for(i=0;i<MX;i++){
			for(j=0;j<MEDIA_PATCH;j++){
			
			if (j == 0){
				if (i == 0) {
					var.b[j+(i+N)*MEDIA_PATCH] = (var.pgxy*var.Riy-1.0)*var.u[0][i+N][j]+var.pgxy*var.Riy*var.u[0][i+N][j+1]+var.pgxy*var.Rix*(var.u[0][i+N][j+2]+var.u[0][i+N][j+3]);
					//printf("m[%d][%d]=%d-->b[%d]\n",N+i,j,10,j+(N+i)*MEDIA_PATCH);
				} else {
					var.b[j+(i+N)*MEDIA_PATCH] = var.pd*var.u[0][i-1+N][j+1]+(var.pgxy*var.Riy-var.pc)*var.u[0][i+N][j] + var.pgxy*var.Riy*var.u[0][i+N][j+1]+var.pgxy*var.Rix*(var.u[0][i+N][j+2]+var.u[0][i+N][j+3]);
					//printf("m[%d][%d]=%d-->b[%d]\n",N+i,j,11,j+(N+i)*MEDIA_PATCH);
				}
			}

			if (j == 1){
				if(i==MX-1) {
					var.b[j+(i+N)*MEDIA_PATCH] = var.pgxy*var.Riy*var.u[0][i+N][j-1] + (var.pgxy*var.Riy-1.0)*var.u[0][i+N][j]+var.pgxy*var.Rix*(var.u[0][N+i][j+1]+var.u[0][N+i][j+2]); 
					//printf("m[%d][%d]=%d-->b[%d]\n",N+i,j,20,j+(N+i)*MEDIA_PATCH);
				} else {
					var.b[j+(i+N)*MEDIA_PATCH] = var.pgxy*var.Riy*var.u[0][N+i][j-1] + (var.pgxy*var.Riy-var.pc)*var.u[0][N+i][j]+var.pgxy*var.Rix*(var.u[0][i+N][j+1]+var.u[0][i+N][j+2])+var.pd*var.u[0][N+i+1][j-1]; 
					//printf("m[%d][%d]=%d-->b[%d]\n",N+i,j,22,j+(N+i)*MEDIA_PATCH);
				}
			}

			if (j == 2){
				if(k==0){
					var.b[j+(i+N)*MEDIA_PATCH] = var.pgxy*var.Riy*(var.u[0][N+i][j-2]+var.u[0][N+i][j-1])+(var.pgxy*var.Rix-1.0)*var.u[0][N+i][j]+var.pgxy*var.Rix*var.u[0][N+i][j+1]; 
					//printf("m[%d][%d]=%d-->b[%d]\n",N+i,j,30,j+(N+i)*MEDIA_PATCH);
				} else {
					var.b[j+(i+N)*MEDIA_PATCH] = var.pgxy*var.Riy*(var.u[0][N+i][j-2]+var.u[0][N+i][j-1])+(var.pgxy*var.Rix-var.pg)*var.u[0][N+i][j]+var.pgxy*var.Rix*var.u[0][N+i][j+1]+var.ph*var.u[0][N-MX+i][j+1]; 
					//printf("m[%d][%d]=%d-->b[%d]\n",N+i,j,33,j+(N+i)*MEDIA_PATCH);
				}
			}

			if (j == 3){
				if(k==MY-1){
					var.b[j+(N+i)*MEDIA_PATCH] = var.pgxy*var.Riy*(var.u[0][N+i][j-3]+var.u[0][N+i][j-2])+var.pgxy*var.Rix*var.u[0][N+i][j-1]+(var.pgxy*var.Rix-1.0)*var.u[0][N+i][j]; 
					//printf("m[%d][%d]=%d-->b[%d]\n",N+i,j,40,j+(N+i)*MEDIA_PATCH);
				} else {
					var.b[j+(N+i)*MEDIA_PATCH] = var.pgxy*var.Riy*(var.u[0][N+i][j-3]+var.u[0][N+i][j-2])+var.pgxy*var.Rix*var.u[0][N+i][j-1]+(var.pgxy*var.Rix-var.pg)*var.u[0][N+i][j]+var.ph*var.u[0][N+MX+i][j-1]; 
					//printf("m[%d][%d]=%d-->b[%d]\n",N+i,j,44,j+(N+i)*MEDIA_PATCH);
				}
			}

			} 

		}
	}
	//exit(0);

}


void linear_solve()
{

	long long n;
	MKL_INT i,j,k;

	n=MAT_SIZE;
	var.phase = 33;
	var.iparm[7] = 2; /* Max numbers of iterative refinement steps. */

	//PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, x, &error);
	PARDISO_64 (var.pt, &var.maxfct, &var.mnum, &var.mtype, &var.phase, &n, var.AT2, var.ia, var.ja, &var.idum, &var.nrhs, var.iparm, &var.msglvl, var.b, var.solv, &var.error);
		if (var.error != 0) {
			printf("\nERROR during solution: %ld", var.error);
			exit(3);
		}

/*	for (i = 0; i < n; i++) {
		printf("\n k=%d x [%d] = %16.15e",i,i,var.solv[i]);
	}
	printf ("\n");
	*/
//	exit(0);

/* -------------------------------------------------------------------- */
/* .. Termination and release of memory. */
/* -------------------------------------------------------------------- */
//	var.phase = -1; /* Release internal memory. */
//	PARDISO (var.pt, &var.maxfct, &var.mnum, &var.mtype, &var.phase, &n, &var.ddum, var.ia, var.ja, &var.idum, &var.nrhs, var.iparm, &var.msglvl, &var.ddum, &var.ddum, &var.error);
//	return 0;

}

