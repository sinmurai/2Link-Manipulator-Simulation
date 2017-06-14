// 2 Link Manipulator Control

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ManipulatorFunc.h"

void main(void)
{
	FILE *fo;
	struct LinkParam Link;
	struct PdParam Cont;
	double x[4];
	// x[0] : theta 1            //
	// x[1] : theta 2            //
	// x[2] : theta velocity 1   //
	// x[3] : theta velocity 2   //
	int i;
	double t,TSTEP;
	//
	//パラメータ初期化・設定
	Link.m1 = 1.0;  Link.m2 = 1.0;
	Link.L1 = 0.1;  Link.L2 = 0.1;
	Link.J1 = 0.01; Link.J2 = 0.01;
	Link.g  = 9.8; 
	//
	Cont.cp1 = 360.0;  Cont.cp2 = 360.0;
	Cont.cd1 = 15.0;  Cont.cd2 = 13.0;
	Cont.dv1 = 30*pi/180.0;
	Cont.dv2 = 30*pi/180.0;
	//
	for (i=0;i<4;i++){
		x[i] = 0.0;
	}
	TSTEP = 0.01;
	//出力ファイル設定
	fo = fopen("outdata.csv","w");
	//ファイル出力
	for (i=0;i<4;i++){
		fprintf(fo,"%f,",x[i]);
	}fprintf(fo,"\n");

	//シミュレーション
	for (t=0;t < LTIME; t += TSTEP){
		//数値積分
		Runge_Kutta(Link,Cont,x,t,TSTEP);
		//ファイル出力
		for (i=0;i<4;i++){
			fprintf(fo,"%f,",x[i]);
		}fprintf(fo,"\n");
	}
	fclose(fo);

}
