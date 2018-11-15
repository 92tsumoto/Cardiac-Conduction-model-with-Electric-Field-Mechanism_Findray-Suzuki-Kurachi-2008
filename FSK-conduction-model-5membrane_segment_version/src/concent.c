#include "syspara.h"

// The units of dnai is in mM.  Note that naiont should be multiplied by the
// cell capacitance to get the correct units.  Since cell capacitance = 1 uF/cm^2,
// it doesn't explicitly appear in the equation below.
// This holds true for the calculation of dki and dcai.

// Na concentration dose not use a function.
// K concentration dose not use a function.
inline void comp_iconcent (int site, int patch)
{

	conc_nsr(site,patch);
	conc_jsr(site,patch);
	conc_itr(site,patch);
	conc_cai(site,patch);

}

inline void conc_nsr (int site, int patch)

// NSR Ca Ion Concentration Changes 
// iupbar = 0.00875;  // Max. current through iup channel (mM/ms)
// iup;               // Ca uptake from myo. to NSR (mM/ms)
// kleak;             // Rate constant of Ca leakage from NSR (ms^-1)
// ileak;             // Ca leakage from NSR to myo. (mM/ms)
// dnsr;              // Change in [Ca] in the NSR (mM)
// kmup = 0.00092;    // Half-saturation concentration of iup (mM)
// nsrbar = 15;       // Max. [Ca] in NSR (mM) 
// x[18] = [Ca]i (cai)
// x[19] = [Ca]_NSR (nsr)
// x[20] = [Ca]_JSR (jsr)

{

	double kleak;

	kleak = var.iupbar/var.nsrbar;
	var.ileak[site][patch] = kleak*var.u[19][site][patch];

	var.iup[site][patch] = var.iupbar*var.u[18][site][patch]/(var.u[18][site][patch]+var.kmup);

	//dnsr = dt*(iup-ileak-itr*vjsr/vnsr);
	//nsr = nsr+dnsr;

}

inline void conc_jsr (int site, int patch)

// JSR Ca Ion Concentration Changes 
// djsr;     // Change in [Ca] in the JSR (mM)
// tauon = 2;        // Time constant of activation of Ca release from JSR (ms)
// tauoff = 2;       // Time constant of deactivation of Ca release from JSR (ms)
// tcicr;            // t=0 at time of CICR (ms)
// irelcicr;         // Ca release from JSR to myo. due to CICR (mM/ms)
// csqnth = 8.75;    // Threshold for release of Ca from CSQN due to JSR ovreload (mM)
// gmaxrel = 150;    // Max. rate constant of Ca release from JSR due to overload (ms^-1)
// grelbarjsrol;     // Rate constant of Ca release from JSR due to overload (ms^-1)
// greljsrol;        // Rate constant of Ca release from JSR due to CICR (ms^-1)
// tjsrol;           // t=0 at time of JSR overload (ms)
// ireljsrol;        // Ca release from JSR to myo. due to JSR overload (mM/ms)
// swspontan = 0;    // switch of spontaneous release
// csqnbar = 10;     // Max. [Ca] buffered in CSQN (mM)
// kmcsqn = 0.8;     // Equalibrium constant of buffering for CSQN (mM)
// bjsr;             // b Variable for analytical computation of [Ca] in JSR (mM)
// cjsr;             // c Variable for analytical computation of [Ca] in JSR (mM)
// on;               // Time constant of activation of Ca release from JSR (ms)
// off;              // Time constant of deactivation of Ca release from JSR (ms)
// magrel;           // Magnitude of Ca release
// dICa;          // Rate of change of Ca entry
// dICa_new;       // New rate of change of Ca entry
// ICa_total_old;        // Old rate of change of Ca entry
// x[18] = [Ca]i (cai)
// x[19] = [Ca]_NSR (nsr)
// x[20] = [Ca]_JSR (jsr)
// t_cicr = time of CICR (tcicr)
// t_overload = time of JSR overload (tjsrol)

{

	double magrel;
	double on,off;
	double greljsrol;
	double test;

	var.dICa_total_new[site][patch] = (var.ICa_total[site][patch] - var.ICa_total_old[site][patch])/var.dt;
	//printf("dICa_total_new=%lf,ICa_total=%lf\n",var.dICa_total_new,var.ICa_total);

	if( var.u[0][site][patch] > -35.0 && var.dICa_total_new[site][patch] < var.dICa_total[site][patch] && var.boolien[site][patch] == 0){
		test=var.t_cicr[site][patch];
		var.boolien[site][patch] = 1;
		var.t_cicr[site][patch] = 0;
		//if(patch==0) printf("reset t_cicr[%d][%d]\n",site,patch);
		//printf("reset t_cicr[%d][%d]=%lf\n",site,patch,test);
	}

	on = 1.0/(1.0+exp((-var.t_cicr[site][patch]+4.0)/0.5));
	off = (1.0-1.0/(1.0+exp((-var.t_cicr[site][patch]+4.0)/0.5)));

	magrel = 1.0/(1.0+exp(((var.ilca[site][patch]+var.icab[site][patch]+var.ipca[site][patch]-2.0*var.inaca[site][patch]+var.icat[site][patch])+5.0)/0.9));

	var.Irel_cicr[site][patch] = var.gmaxrel*on*off*magrel*(var.u[20][site][patch]-var.u[18][site][patch]);
	
	var.t_cicr[site][patch] += var.dt;

	greljsrol = var.grelbarjsrol*(1.0-exp(-var.t_overload[site][patch]/var.tauon))*exp(-var.t_overload[site][patch]/var.tauoff);
	var.Irel_jsr_overload[site][patch] = greljsrol*(var.u[20][site][patch]-var.u[18][site][patch]);
	
	var.t_overload[site][patch] += var.dt;

	var.csqn[site][patch] = var.csqnbar*(var.u[20][site][patch]/(var.u[20][site][patch] + var.kmcsqn));

	//var.djsr = var.dt*(var.itr-var.Irel_cicr-var.Irel_jsr_overload);
	//var.bjsr = var.csqnbar - var.csqn - var.djsr - var.jsr + var.kmcsqn;
	//var.cjsr = var.kmcsqn*(var.csqn + var.djsr + var.jsr);
	//var.jsr = (sqrt(var.bjsr*var.bjsr + 4.0*var.cjsr) - var.bjsr)/2.0;
}

inline void conc_itr (int site, int patch)

// Translocation of Ca Ions from NSR to JSR
// itr;      // Translocation current of Ca ions from NSR to JSR (mM/ms)
// tautr = 180;      // Time constant of Ca transfer from NSR to JSR (ms)
// x[19] = [Ca]_NSR (nsr)
// x[20] = [Ca]_JSR (jsr)
{
	var.itr[site][patch] = (var.u[19][site][patch]-var.u[20][site][patch])/var.tautr; 
}

inline void conc_cai (int site, int patch)

// Myoplasmic Ca Ion Concentration Changes 
// dcai;    // Change in myoplasmic Ca concentration (mM)
// catotal; // Total myoplasmic Ca concentration (mM)
// bmyo;    // b Variable for analytical computation of [Ca] in myoplasm (mM)
// cmyo;    // c Variable for analytical computation of [Ca] in myoplasm (mM)
// dmyo;    // d Variable for analytical computation of [Ca] in myoplasm (mM)
// gpig;    // Tribute to all the guinea pigs killed for the advancement of knowledge
// cmdnbar = 0.050;   // Max. [Ca] buffered in CMDN (mM)
// trpnbar = 0.070;   // Max. [Ca] buffered in TRPN (mM)
// kmcmdn = 0.00238;  // Equalibrium constant of buffering for CMDN (mM)
// kmtrpn = 0.0005;   // Equalibrium constant of buffering for TRPN (mM)

{

	//double trpn,trpnbar,cmdn;
	//dcai/dt = -(((var.caiont*var.acap)/(var.vmyo*var.zca*F))
	//	+((var.iup-var.ileak)*var.vnsr/var.vmyo)
	//  -(var.irelcicr*var.vjsr/var.vmyo)-(var.ireljsrol*var.vjsr/var.vmyo));
	var.trpn[site][patch] = var.trpnbar*(var.u[18][site][patch]/(var.u[18][site][patch]+var.kmtrpn));
	var.cmdn[site][patch] = var.cmdnbar*(var.u[18][site][patch]/(var.u[18][site][patch]+var.kmcmdn));

	//catotal = trpn + cmdn + dcai + x[20];

	//bmyo = var.cmdnbar + var.trpnbar + var.kmtrpn + var.kmcmdn - catotal;

	//cmyo = (var.kmcmdn*var.kmtrpn)
	//		-(catotal*(var.kmtrpn + var.kmcmdn))
	//		+(var.trpnbar*var.kmcmdn)
	//		+(var.cmdnbar*var.kmtrpn);

	//dmyo = -var.kmtrpn*var.kmcmdn*catotal;

	//gpig = sqrt(bmyo*bmyo-3.0*cmyo);

	//var.cai = (2*gpig/3)*cos(acos((9*bmyo*cmyo-2*bmyo*bmyo*bmyo-27*dmyo)/
	//			(2*pow((bmyo*bmyo-3*cmyo),1.5)))/3)-(bmyo/3);

}

inline void comp_iconcent2 (int site, int patch)
{

// Update Ca_JSR concentration x[20] = [Ca]_JSR
	var.bjsr[site][patch] = var.csqnbar - var.csqn[site][patch] - var.u[20][site][patch] + var.kmcsqn;
	var.cjsr[site][patch] = var.kmcsqn*(var.csqn[site][patch] + var.u[20][site][patch]);
	var.u[20][site][patch] = (sqrt(var.bjsr[site][patch]*var.bjsr[site][patch]+4.0*var.cjsr[site][patch])-var.bjsr[site][patch])/2.0;

// Update intracellular Ca concentration x[18] = [Ca]_i

	var.Ca_total[site][patch] = var.trpn[site][patch]+var.cmdn[site][patch]+var.u[18][site][patch];
	var.bmyo[site][patch] = var.cmdnbar+var.trpnbar+var.kmtrpn+var.kmcmdn-var.Ca_total[site][patch];
	var.cmyo[site][patch] = (var.kmcmdn*var.kmtrpn)-(var.Ca_total[site][patch]*(var.kmtrpn+var.kmcmdn))+(var.trpnbar*var.kmcmdn)+(var.cmdnbar*var.kmtrpn);
	var.dmyo[site][patch] = -var.kmtrpn*var.kmcmdn*var.Ca_total[site][patch];
	var.gpig = sqrt(var.bmyo[site][patch]*var.bmyo[site][patch]-3.0*var.cmyo[site][patch]);

	var.u[18][site][patch] = (2.0*var.gpig/3.0)*cos(acos((9.0*var.bmyo[site][patch]*var.cmyo[site][patch]
			-2.0*var.bmyo[site][patch]*var.bmyo[site][patch]*var.bmyo[site][patch]-27.0*var.dmyo[site][patch])/
			(2.0*pow((var.bmyo[site][patch]*var.bmyo[site][patch]-3.0*var.cmyo[site][patch]),1.5)))/3.0)-(var.bmyo[site][patch]/3.0);

}

