#include "syspara.h"

void static_paras(FILE *fp1)
{
	int i,j,k;

	// general
		var.RTonF = R*T/F;

	// Cell Geometry */

		var.RGC = 2;  // The ratio of the actual surface area to the geometrical surface area (cm)
		var.length = var.length;  // Length of the cell (cm)
		var.a = var.a;     // Radius of the cell (cm)

		//var.vcell = 1000*M_PI*var.a*var.a*var.length; // Cell Volume:3.801e-5 (uL)
		if(MEDIA_PATCH==5){
			for(i=0;i<MEDIA_PATCH;i++){
				if(i==0 || i==MEDIA_PATCH-1) {
					var.vcell[i] = 1000*M_PI*var.a*var.a*((var.length/(MEDIA_PATCH-2))*(1.0/4.0)); // (length*(1/(unit-2)*(1/4)))
				} else if (i==1 || i== MEDIA_PATCH -2 ){
					var.vcell[i] = 1000*M_PI*var.a*var.a*((var.length/(MEDIA_PATCH-2))*(3.0/4.0)); // (lenght*(1/(unit-2)*(3/4)))
				} else {
					var.vcell[i] = 1000*M_PI*var.a*var.a*var.length/(MEDIA_PATCH-2); // (lenght*(1/(unit-2)))
				}
			}
		} else if(MEDIA_PATCH==3){
			for(i=0;i<MEDIA_PATCH;i++){
				if(i==0 || i==MEDIA_PATCH-1) {
					var.vcell[i] = 1000*M_PI*var.a*var.a*var.length*(1.0/12.0); // (length*(1/(unit-2)*(1/4)))
				} else {
					var.vcell[i] = 1000*M_PI*var.a*var.a*var.length*(5.0/6.0); // (lenght*(1/(unit-2)))
				}
			}
		}
		
	// geometric membrane area: 7.671e-5 (cm^2)
		//var.ageo = M_PI*var.a*var.a+2.0*M_PI*var.a*var.length;

		for(i=0;i<MEDIA_PATCH;i++){
			if(i==0 || i == MEDIA_PATCH-1){
				var.ageo[i] = var.Srate*M_PI*var.a*var.a;
			} else {
				var.ageo[i] = 2.0*M_PI*var.a*var.length/(MEDIA_PATCH-2);
			}
		}
		
		for(i=0;i<MEDIA_PATCH;i++){
			var.acap[i] = var.ageo[i]*var.RGC;          // Capacitive membrane area: 1.534e-4 cm^2 (cm^2)
			var.vmyo[i] = var.vcell[i]*0.68;      // Myoplasm volume (uL) = 68% for Cell volume
			var.vmito[i] = var.vcell[i]*0.26;     // Mitochondria volume (uL) = 26% for cell volume
			var.vsr[i] = var.vcell[i]*0.06;       // SR volume (uL)
			var.vnsr[i] = var.vcell[i]*0.0552;    // NSR volume (uL)
			var.vjsr[i] = var.vcell[i]*0.0048;    // JSR volume (uL)
			var.vcleft[i] = var.vcell[i]*0.12/0.88;  // Cleft volume (uL)
		}

	// Ion Valences
		var.zna = 1;  // Na valence
		var.zk = 1;   // K valence
		var.zca = 2;  // Ca valence
	
	// Extracellular ion concentrations
		//var.nao = 132.0;     // Na (mM) Original value
		//var.nao = 140.0;     // Na (mM) Correct value
		//var.ko = 4.5;      // K (mM)
		//var.cao = 1.8;     // Ca (mM)
		var.Na_bulk = 140.0;     // Na (mM) Correct value
		var.K_bulk = 5.4;      // K (mM)
		var.Ca_bulk = 1.8;     // Ca (mM)

	// Max conductanve (constant)
		var.gna = 16;	// (mS/uF);
		var.gcat = 0.05;	// (mS/uF);
		var.gkna = 0.12848;
		var.kdkna = 66.0; 
	
		switch(var.na_switch){
			case 1: // heterogeneous %
				var.gna_all = var.narate*var.gna*(2.0*M_PI*var.a*var.a + 2.0*M_PI*var.a*var.length);
				for(i=0;i<MEDIA_SITE;i++){
					for(j=0;j<MEDIA_PATCH;j++){
						if (j==0||j==MEDIA_PATCH-1){
							var.gna_local[i][j] = (0.5*var.jmrate1*var.gna_all)/(M_PI*var.a*var.a);
						} else {
							var.gna_local[i][j] = (var.lmrate1*var.gna_all/(MEDIA_PATCH-2))/(2.0*M_PI*var.a*var.length/(MEDIA_PATCH-2));
						}
					}
				}
				break;
			case 2: // uniform
				var.gna_all = var.narate*var.gna*(2.0*M_PI*var.a*var.a + 2.0*M_PI*var.a*var.length);
				for(i=0;i<MEDIA_SITE;i++){
					for(j=0;j<MEDIA_PATCH;j++){
						var.gna_local[i][j] = var.narate*var.gna;
					}
				}
				break;
			case 3: // heterogeneous myofiber and local differences of NaCh in a myofiber case  %
				var.gna_all = var.narate*var.gna*(2.0*M_PI*var.a*var.a + 2.0*M_PI*var.a*var.length);
				for(i=0;i<30;i++){
					for(j=0;j<MEDIA_PATCH;j++){
						if (j==0||j==MEDIA_PATCH-1){
							var.gna_local[i][j] = (0.5*0.3*var.gna_all)/(M_PI*var.a*var.a); 
						} else {
							var.gna_local[i][j] = (0.7*var.gna_all/(MEDIA_PATCH-2))/(2.0*M_PI*var.a*var.length/(MEDIA_PATCH-2));
						}
					}
				}
				for(i=30;i<150;i++){
					for(j=0;j<MEDIA_PATCH;j++){
						if (j==0||j==MEDIA_PATCH-1){
							var.gna_local[i][j] = (0.5*var.jmrate1*var.gna_all)/(M_PI*var.a*var.a); 
						} else {
							var.gna_local[i][j] = (var.lmrate1*var.gna_all/(MEDIA_PATCH-2))/(2.0*M_PI*var.a*var.length/(MEDIA_PATCH-2));
						}
					}
				}
				for(i=150;i<MEDIA_SITE;i++){
					for(j=0;j<MEDIA_PATCH;j++){
						if (j==0||j==MEDIA_PATCH-1){
							var.gna_local[i][j] = (0.5*var.jmrate2*var.gna_all)/(M_PI*var.a*var.a); 
						} else {
							var.gna_local[i][j] = (var.lmrate2*var.gna_all/(MEDIA_PATCH-2))/(2.0*M_PI*var.a*var.length/(MEDIA_PATCH-2));
						}
					}
				}
				break;
			case 4: // homogeneous myofiber and local differences of NaCh in a myofiber case  %
            	break;
		}

		for(i=0;i<MEDIA_SITE;i++){
			for(k=0;k<MEDIA_PATCH;k++){
				printf("gna[%d][%d]=%lf\n",i,k,var.gna_local[i][k]);
				fprintf(fp1,"gna[%d][%d]=%e\n",i,k,var.gna_local[i][k]);
			}
		}
		
	// classification by cell type
	//	var.gkr = 0.02614;                      // real value: gkr*ikrf
	//	var.gkr_max = var.gkr*var.ikrf;
	//	var.gks = 0.433;                        // real value: gks_max*iksf
	//	var.gks_max = var.gks*var.iksf; 
	//	var.gkp = 0.02;                         // real value: gkp_max*ikxf
	//	var.gkp_max = var.gkp*var.ikxf; 
	//	var.gitodv = 0.5;                       // real value: gitodv*itof
	//	var.gitodv_max = var.gitodv*var.itof; 
		
	switch(var.apd_switch){
		case 1: // Only Endo cell fiber
			for(i=0;i<MEDIA_SITE;i++){
				for(j=0;j<MEDIA_PATCH;j++){
					var.ikrf[i][j] = 0.6400;
					var.iksf[i][j] = 0.1903;
					var.ikxf[i][j] = 1.2699;
					var.itof[i][j] = 0.1131;
				}
			}
			printf("ikrf[0]=%lf,iksf[0]=%lf,ikxf[0]=%lf,itof[0]=%lf\n",
				var.ikrf[0][2],var.iksf[0][2],var.ikxf[0][2],var.itof[0][2]);
			fprintf(fp1,"ikrf[0]=%lf,iksf[0]=%lf,ikxf[0]=%lf,itof[0]=%lf\n",
				var.ikrf[0][2],var.iksf[0][2],var.ikxf[0][2],var.itof[0][2]);
			break;
		case 2: // Only M cell fiber
			for(i=0;i<MEDIA_SITE;i++){
				for(j=0;j<MEDIA_PATCH;j++){
					var.ikrf[i][j] = 0.5382;
					var.iksf[i][j] = 0.0951;
					var.ikxf[i][j] = 1.2699;
					var.itof[i][j] = 0.1600;
				}
			}
			printf("ikrf[0]=%lf,iksf[0]=%lf,ikxf[0]=%lf,itof[0]=%lf\n",
				var.ikrf[0][2],var.iksf[0][2],var.ikxf[0][2],var.itof[0][2]);
			fprintf(fp1,"ikrf[0]=%lf,iksf[0]=%lf,ikxf[0]=%lf,itof[0]=%lf\n",
				var.ikrf[0][2],var.iksf[0][2],var.ikxf[0][2],var.itof[0][2]);
			break;
		case 3: // Only Epi cell fiber
			for(i=0;i<MEDIA_SITE;i++){
				for(j=0;j<MEDIA_PATCH;j++){
					var.ikrf[i][j] = 1.2800;
					var.iksf[i][j] = 0.1903;
					var.ikxf[i][j] = 0.5040;
					var.itof[i][j] = 0.2263;
				}
			}
			printf("ikrf[0]=%lf,iksf[0]=%lf,ikxf[0]=%lf,itof[0]=%lf\n",
				var.ikrf[0][2],var.iksf[0][2],var.ikxf[0][2],var.itof[0][2]);
			fprintf(fp1,"ikrf[0]=%lf,iksf[0]=%lf,ikxf[0]=%lf,itof[0]=%lf\n",
				var.ikrf[0][2],var.iksf[0][2],var.ikxf[0][2],var.itof[0][2]);
			break;
	}	


	// L-type calcium current
		var.kmca = 0.0006;     // Half-saturation concentration of Ca channel (mM)
		var.pca = 0.00054;     // Permiability of membrane to Ca (cm/s)
		var.gacai = 1;         // Activity coefficient of Ca
		var.gacao = 0.341;     // Activity coefficient of Ca
		var.pna = 0.000000675; // Permiability of membrane to Na (cm/s)
		var.pk = 0.000000193;  // Permiability of membrane to K (cm/s)
		//var.ratgca = 0.5;      // rate ?
		var.fcarat = 0.25;	   // rate of the slow gate during CDI process
		var.aCDI = 0.00128;	   // rate constant of inactivation
		var.bCDI = 0.003;	   // rate constant of recovery	

	// Rapid Activated Potassium Current: IKr
		//var.gkr = 0.02614*sqrt(var.ko/5.4);

	// Slow Activated Potassium Current: IKs
		var.prnak = 0.01833;     // Na/K Permiability Ratio

	// Inward rectifier potassium current: IK1
		//var.gki = 0.75*sqrt(var.ko/5.4);
		var.c1_ki = exp(-0.2385*59.215);
		var.c2_ki = exp(0.08032*5.476);
		var.c3_ki = exp(-0.06175*594.31);
		var.c4_ki = exp(-0.5143*4.753);

	// Sodium-Calcium Exchanger V-S
		var.c1 = 0.00025;   // Scaling factor for inaca (uA/uF)
		var.c2 = 0.0001;    // Half-saturation concentration of NaCa exhanger (mM)
		var.gammas = 0.15;  // Position of energy barrier controlling voltage dependance of inaca

	// Sodium-Potassium Pump
		var.ibarnak = 2.25;   // Max. current through Na-K pump (uA/uF)
		var.kmnai = 10;    // Half-saturation concentration of NaK pump (mM)
		var.kmko = 1.5;    // Half-saturation concentration of NaK pump (mM)
		//var.sigma = (exp(var.nao/67.3)-1.0)/7.0;

	// Nonspecific Ca-activated Current
		var.pnsca = 0.000000175;   // Permiability of channel to Na and K (cm/s)
		var.kmnsca = 0.0012;       // Half-saturation concentration of NSCa channel (mM)
		var.ganai = 0.75;      // Activity coefficient of Na
		var.ganao = 0.75;      // Activity coefficient of Na
		var.gaki = 0.75;       // Activity coefficient of K
		var.gako = 0.75;       // Activity coefficient of K

	
	// Sarcolemmal Ca Pump
		var.ibarpca = 1.15;  // Max. Ca current through sarcolemmal Ca pump (uA/uF)
		var.kmpca = 0.0005;  // Half-saturation concentration of sarcolemmal Ca pump (mM)

	// Ca Background Current 
		var.gcab = 0.003016; // Max. conductance of Ca background (mS/uF)

	// Na Background Current 
		var.gnab = 0.004;    // Max. conductance of Na background (mS/uF)

	// NSR Ca Ion Concentration Changes 
		var.kmup = 0.00092;   // Half-saturation concentration of iup (mM)
		var.iupbar = 0.00875; // Max. current through iup channel (mM/ms)
		var.nsrbar = 15;      // Max. [Ca] in NSR (mM)

	// JSR Ca Ion Concentration Changes 
		var.tauon = 2.0;        // Time constant of activation of Ca release from JSR (ms)
		var.tauoff = 2.0;       // Time constant of deactivation of Ca release from JSR (ms)
		var.csqnth = 8.75;    // Threshold for release of Ca from CSQN due to JSR ovreload (mM)
		var.gmaxrel = 150;    // Max. rate constant of Ca release from JSR due to overload (ms^-1)
		var.grelbarjsrol = 0; // Rate constant of Ca release from JSR due to overload (ms^-1)
		var.swspontan = 0;    // switch of spontaneous release
		var.csqnbar = 10;     // Max. [Ca] buffered in CSQN (mM)
		var.kmcsqn = 0.8;     // Equalibrium constant of buffering for CSQN (mM)

		for (i=0;i<MEDIA_SITE;i++){
			for (k=0;k<MEDIA_PATCH;k++){
			var.csqn[i][k] = 6.97978;
			var.trpn[i][k] = 0.0143923;
			var.cmdn[i][k] = 0.00257849;
	
		// Another parameter initial setting
			var.t_cicr[i][k] = 25;
			var.t_overload[i][k] = 25;
			var.boolien[i][k] = 1;
			var.dICa_total[i][k] = 0;
			}
		}

	// Translocation of Ca Ions from NSR to JSR
		var.tautr = 180;      // Time constant of Ca transfer from NSR to JSR (ms)

	// Myoplasmic Ca Ion Concentration Changes 
		var.cmdnbar = 0.050;   // Max. [Ca] buffered in CMDN (mM)
		var.trpnbar = 0.070;   // Max. [Ca] buffered in TRPN (mM)
		var.kmcmdn = 0.00238;  // Equalibrium constant of buffering for CMDN (mM)
		var.kmtrpn = 0.0005;   // Equalibrium constant of buffering for TRPN (mM)

	// Extracellular Ion Concentration Changes 
		var.tau_diff = 1000; // Diffusion Constant for Ion Movement from Bulk Medium to Cleft Space

}

