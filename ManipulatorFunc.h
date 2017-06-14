#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define pi 3.1416
#define POW2(x) x*x
#define LTIME 50 //シミュレーション最終時間 
 


//リンクパラメータ
struct LinkParam{
	double m1, m2;    //Mass
	double L1, L2;    //Link length
	double J1, J2;    //Inertia
	double g;         //Gravity
};
//PDコントローラパラメータ
struct PdParam{
	double cp1, cp2;
	double cd1, cd2;
	double dv1, dv2;  //desired value
};


double Traj_sin(double t, double dv)
{// tf : 目標値追従時間
 // t  : 現在の時間
 // dv : 目標値
	double tf;
	//パラメータ設定
	tf = LTIME * 0.8;//最終時間よりちょっと早めに設定
	if (t > tf)t = tf;
	//
	return (dv*(1+sin(3.0/2.0*pi+t/tf*pi))/2.0);
}

//ロボットアーム微分方程式
void Dif_eq(struct LinkParam model, struct PdParam Controller, double x[], double dotx[], double t)
{
	double m1,m2,L1,L2,Lc1,Lc2,J1,J2,g;
	double cd1,cd2,cp1,cp2;
	double dv1,dv2;
	//
	double M[2][2];
	double det;
	double C1,C2;
	double K1,K2;
	double T1,T2;
	//パラメータ設定
	m1 = model.m1;   m2 = model.m2;
	L1 = model.L1;   L2 = model.L2; 
	Lc1 = L1/2.0;    Lc2 = L2/2.0;
	J1 = model.J1;   J2 = model.J2;
	g  = model.g;
	cp1 = Controller.cp1; cp2 = Controller.cp2;
	cd1 = Controller.cd1; cd2 = Controller.cd2;
	dv1 = Traj_sin(t,Controller.dv1); dv2 = Traj_sin(t, Controller.dv2);
	//慣性項
	M[0][0] = J1+J2+m1*POW2(Lc1)+m2*POW2(L1)+m2*POW2(Lc2)+2*m2*L1*Lc2*cos(x[1]);
	M[0][1] = J2+m2*POW2(Lc2)+m2*L1*Lc2*cos(x[1]); 
	M[1][0] = J2+m2*POW2(Lc2)+m2*L1*Lc2*cos(x[1]); 
	M[1][1] = J2+m2*POW2(Lc2);
	det     = M[0][0]*M[1][1];
	det    -= M[0][1]*M[1][0];
	//コリオリ項
	C1      = -m2*L1*Lc2*sin(x[1])*(2*x[2]*x[3]+POW2(x[3]));
	C2      = m2*POW2(x[2])*L1*Lc2*sin(x[1]);
	//重力項
	K1      = -m1*g*Lc1*sin(x[0])-m2*g*(L1*sin(x[0])+Lc2*sin(x[0]+x[1]));
	K2      = -m2*g*Lc2*sin(x[0]+x[1]);
	//コントローラ
	T1      = cp1*(dv1 - x[0]) - cd1;
	T2      = cp2*(dv2 - x[1]) - cd2;

	//運動方程式
	dotx[0] = x[2];
	dotx[1] = x[3];
	dotx[2] = ((M[1][1])*(-C1-K1+T1) + (-M[0][1]*(-C2-K2+T2))) / det;
	dotx[3] = ((-M[1][0])*(-C1-K1+T1) + (M[0][0]*(-C2-K2+T2))) / det;
}

//ルンゲクッタ法
void Runge_Kutta(struct LinkParam model, struct PdParam Controller, double x[], double t, double TSTEP)
{
	//パラメータ設定
	double k1[4], k2[4], k3[4], k4[4];
	double dotx[4];
	double Rx[4], REt, REx[4];
	int i;
	// 初期化
	for (i=0;i<4;i++){
		k1[i] = 0.0; k2[i] = 0.0; k3[i] = 0.0; k4[i] = 0.0;
		dotx[i] = 0.0; Rx[i] = 0.0; REx[i] = 0.0;
	}
	REt = t; 
	
	// k1[] 
	Dif_eq(model, Controller, x, dotx, REt); 
	for(i=0;i<4;i++){ 
		k1[i] = TSTEP*dotx[i];
		REx[i] = x[i] + k1[i]/2.0;
	}
	REt = t + TSTEP/2.0;
	
	// k2[]
	Dif_eq(model, Controller, REx, dotx, REt); 
	for(i=0;i<4;i++){ 
		k2[i] = TSTEP*dotx[i];
		REx[i] = x[i] + k2[i]/2.0;
	}
	REt = t + TSTEP/2.0;

	// k3[]
	Dif_eq(model, Controller, REx, dotx, REt); 
	for(i=0;i<4;i++){ 
		k3[i] = TSTEP*dotx[i];
		REx[i] = x[i] + k3[i];
	}
	REt = t + TSTEP;

	// k4[] 
	Dif_eq(model, Controller, REx, dotx, REt); 
	for(i=0;i<4;i++){ 
		k4[i] = TSTEP*dotx[i];
	}

	// change x[] 
	for(i=0;i<4;i++){ 
		x[i] = x[i] + (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i])/6.0;
	}
}