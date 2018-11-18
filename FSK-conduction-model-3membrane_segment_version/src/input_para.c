#include <string.h>
#include "syspara.h"

void input_para(FILE *fpin)
{
	int i,ii;

	fscanf(fpin,"%d",&var.na_switch);
	fscanf(fpin,"%lf",&var.jmrate1);
	fscanf(fpin,"%lf",&var.lmrate1);
	fscanf(fpin,"%lf",&var.jmrate2);
	fscanf(fpin,"%lf",&var.lmrate2);
	//fscanf(fpin,"%lf",&var.vshift);
	fscanf(fpin,"%lf",&var.a);
	fscanf(fpin,"%lf",&var.length);
	fscanf(fpin,"%lf",&var.Rmyo);
	fscanf(fpin,"%lf",&var.Gj);
	fscanf(fpin,"%lf",&var.cleft_R);
	fscanf(fpin,"%lf",&var.cleft_width);
	fscanf(fpin,"%lf",&var.BCL);
	fscanf(fpin,"%lf",&var.Istim_base);
	fscanf(fpin,"%d",&var.apd_switch);
	fscanf(fpin,"%d",&var.datas);
	for (ii = 0; ii < var.datas; ii++){
		for (i=0;i<NN;i++){
			fscanf(fpin,"%lf",&var.x0[ii][i]);
		}
		fscanf(fpin, "%lf", &var.tsign[ii]);
		fscanf(fpin, "%lf", &var.tend[ii]);
	}
	fscanf(fpin,"%lf",&var.narate);
	fscanf(fpin,"%lf",&var.Srate);
	fscanf(fpin,"%lf",&var.ratgca);
	fscanf(fpin,"%d",&var.l);
	fscanf(fpin,"%d",&var.m);
	fscanf(fpin,"%d",&var.write);
	fscanf(fpin,"%d",&var.out_data);

}

