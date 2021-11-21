/* Jean Philippe EIMER - 1999 *
 *  phil.eimer@9online.fr     */

#include <math.h>

#include <iostream>

#include "include.h"

using namespace std;


void	Orbite::Init(Astre *pAst)
{
	static const double Vorb[9][6][4]={ \
		{{MER_A,0,0,0},MER_VE,MER_VI,MER_VL,MER_VO,MER_VOB}, \
		{{VEN_A,0,0,0},VEN_VE,VEN_VI,VEN_VL,VEN_VO,VEN_VOB}, \
		{{MAR_A,0,0,0},MAR_VE,MAR_VI,MAR_VL,MAR_VO,MAR_VOB}, \
		{{JUP_A,0,0,0},JUP_VE,JUP_VI,JUP_VL,JUP_VO,JUP_VOB}, \
		{{SAT_A,0,0,0},SAT_VE,SAT_VI,SAT_VL,SAT_VO,SAT_VOB}, \
		{{URA_A,0,0,0},URA_VE,URA_VI,URA_VL,URA_VO,URA_VOB}, \
		{{NEP_A,0,0,0},NEP_VE,NEP_VI,NEP_VL,NEP_VO,NEP_VOB}, \
		{{PLU_A,0,0,0},PLU_VE,PLU_VI,PLU_VL,PLU_VO,PLU_VOB}, \
		{{SOL_A,0,0,0},SOL_VE,SOL_VI,SOL_VL,SOL_VOB,SOL_VOB}};

	
	pA=pAst;

	switch(pA->astreid)
	{
		case	SOL:
			*this=Vorb[PLU+1];
			break;

		case	TER:
			*this=pA->pE->Soleil.ElemOrbit;
			L+=180.0;
			L=K360(L);
			Omega+=180.0;
			Omega=K360(Omega);
			omegab+=180.0;
			omegab=K360(omegab);
			break;

		default:
			*this=Vorb[pA->astreid];
	}
}


// Initialisation des vecteurs
Orbite	Orbite::operator =(const double Vinit[6][4])
{
	a=Vinit[0][0];
	Ve=Vinit[1];
	Vi=Vinit[2];
	VL=Vinit[3];
	VO=Vinit[4];
	Vob=Vinit[5];

	return(*this);
}


// Copie d'une structure Orbite
Orbite	&Orbite::operator =(const Orbite &O2)
{
	pA=O2.pA;
	
	Ve=O2.Ve;
	Vi=O2.Vi;
	VL=O2.VL;
	VO=O2.VO;
	Vob=O2.Vob;

	a=O2.a;
	e=O2.e;
	i=O2.i;
	L=O2.L;
	Omega=O2.Omega;
	omegab=O2.omegab;
	deltaL=O2.deltaL;

	aper=O2.aper;
	eper=O2.eper;
	Lper=O2.Lper;
	vkper=O2.vkper;
	rper=O2.rper;
	lper=O2.lper;
	bper=O2.bper;

	M=O2.M;
	MRsp=O2.MRsp;
	omega=O2.omega;

	E=O2.E;
	ER=O2.ER;

	r=O2.r;
	v=O2.v;

	return(*this);
}


// Calcul des �l�ments orbitaux en fonction du temps
void	Orbite::Elements(double t)
{
	// calcul � partir des vecteurs
	e=Ve.DotT(t);
	i=Vi.DotT(t);
	L=VL.DotTprecis(t);
	Omega=VO.DotT(t);
	omegab=Vob.DotT(t);
	deltaL=VL.d[1]*360.0+VL.d[2]+VL.d[3];

	L=K360(L);
	Omega=K360(Omega);
	omegab=K360(omegab);

	MRsp=DEG2RAD(L-omegab);	// anomalie moyenne sans perturbation, en radians
}


void	Orbite::Anomalies(void)
{
	// application des perturbations
	a+=aper;
	omegab+=vkper/e;
	e+=eper;
	L+=Lper;

	// calcul des autres �l�ments avec perturbations
	AnomalieMoyenne();
	Kepler();
	AnomalieVraie();
}


// calcul de M et omega
double	Orbite::AnomalieMoyenne(void)
{
	M=L-omegab;
	omega=omegab-Omega;

	M=K360(M);
	omega=K360(omega);

	return(M);
}


// Equation de Kepler
double	Orbite::Kepler(void)
{
	double	MR=DEG2RAD(M);
	double	E2=ER;

	do
	{
		ER=E2;
		E2=ER-(ER-e*sin(ER)-MR)/(1.0-e*cos(ER));
	}
	while(fabs(E2-ER) > PRECISION_KEPLER);

	E=RAD2DEG(E2);
	E=K360(E);

	ER=DEG2RAD(E);

	return(E);
}


// calcul du rayon vecteur et de l'anomalie vraie
double	Orbite::AnomalieVraie(void)
{
	double	x,y;

	x=a*(cos(ER)-e);	// coordonn�es rectangulaire
	y=a*sqrt(1-e*e)*sin(ER);

	r=sqrt(x*x+y*y);	// coordonn�es polaires
	v=RAD2DEG(atan(y/x));

	if(x<0.0)		// lever d'ind�termination suite � atan
		v+=180.0; 	

	v=K360(v);		// 0-360�

	return(v);
}


// Calcul des perturbations sur les �l�ments orbitaux
void	Orbite::Perturbations(double t)
{
	aper=eper=Lper=vkper=rper=lper=bper=0;	// Pas de perturbation par d�faut

	if(pA->astreid==SOL)
		PerturbationsSoleil(t);
	else if(pA->astreid<JUP)
		PerturbationsTellur(t);
	else if(pA->astreid<PLU)
		PerturbationsGeantes(t);
}

void	Orbite::PerturbationsSoleil(double t)
{
	static const VectData	Ve[4]=SOL_PER_VVE;
	static const VectData	Vper[2]={{SOL_PER_VL},{SOL_PER_VB}};
	static const double	Vbc[5]=SOL_PER_VBC;
	static double	Vlc[5]=SOL_PER_VLC;
	double	e[4];


	e[0]=DEG2RAD(Ve[0].DotTprecis(t));
	e[1]=DEG2RAD(Ve[1].DotTprecis(t));
	e[2]=DEG2RAD(Ve[2].DotTprecis(t));
	e[3]=DEG2RAD(Ve[3].DotTprecis(t));

	rper=Vbc[0]*sin(e[0])+Vbc[1]*sin(e[1])+Vbc[2]*sin(e[2])+Vbc[3]*cos(e[3])+Vbc[4]*sin(DEG2RAD(Vper[1].DotTprecis(t)));
	lper=Vlc[0]*cos(e[0])+Vlc[1]*cos(e[1])+Vlc[2]*cos(e[2])+Vlc[3]*sin(e[3])+Vlc[4]*sin(DEG2RAD(Vper[0].DotT(t)));
}

void	Orbite::PerturbationsTellur(double t)
{
	static const VectData	Vmer[8]=MER_PER_VV,Vven[12]=VEN_PER_VV,Vven2=VEN_PER_V;
	static const VectData	Vmar[21]=MAR_PER_VV,Vmar1=MAR_PER_V1,Vmar2=MAR_PER_V2;
	static const double	Cmer[8]=MER_PER_VC,Cven[12]=VEN_PER_VC;
	static const double	Cmar[21]=MAR_PER_VC1,Cmar2[2]=MAR_PER_VC2;
	double	aa;
	unsigned short int	i;
	VectData	VM;


	VM.d[0]=DEG2RAD(pA->pE->Planete[JUP].ElemOrbit.VL.DotTprecis(t)-pA->pE->Planete[JUP].ElemOrbit.Vob.DotT(t));	// anomalie moyenne de Jupiter
	VM.d[2]=MRsp;	// anomalie moyenne de Mercure si MER, de V�nus si VEN, de Mars si MAR
	VM.d[3]=1.0;

	switch(pA->astreid)
	{
		case	MER:
			VM.d[1]=DEG2RAD(pA->pE->Planete[VEN].ElemOrbit.VL.DotTprecis(t)-pA->pE->Planete[VEN].ElemOrbit.Vob.DotT(t));	// anomalie moyenne de Venus
			for(i=0;i<4;i++)
				lper+=Cmer[i]*cos(VM.Dot(Vmer[i]));
			for(i=4;i<8;i++)
				rper+=Cmer[i]*cos(VM.Dot(Vmer[i]));
			break;

		case	VEN:
			VM.d[1]=DEG2RAD(pA->pE->Soleil.ElemOrbit.VL.DotTprecis(t)-pA->pE->Soleil.ElemOrbit.Vob.DotT(t));	// anomalie moyenne du  Soleil
			for(i=0;i<5;i++)
				lper+=Cven[i]*cos(VM.Dot(Vven[i]));
			for(i=5;i<12;i++)
				rper+=Cven[i]*cos(VM.Dot(Vven[i]));

			Lper=VEN_PER_C*sin(Vven2.DotT(t));
			break;

		case	MAR:
			VM.d[1]=DEG2RAD(pA->pE->Soleil.ElemOrbit.VL.DotTprecis(t)-pA->pE->Soleil.ElemOrbit.Vob.DotT(t));	// anomalie moyenne du  Soleil
			for(i=0;i<8;i++)
				lper+=Cmar[i]*cos(VM.Dot(Vmar[i]));

			lper+=MAR_PER_C1*cos(VM.Dot(Vmar1)-DEG2RAD(pA->pE->Planete[VEN].ElemOrbit.VL.DotTprecis(t)-pA->pE->Planete[VEN].ElemOrbit.Vob.DotT(t)));

			for(i=8;i<21;i++)
				rper+=Cmar[i]*cos(VM.Dot(Vmar[i]));

			aa=VM.Dot(Vmar2);
			Lper=Cmar2[0]*sin(aa)+Cmar2[1]*cos(aa);
			break;
	}
}


void	Orbite::PerturbationsGeantes(double t)
{
	static const VectData	VJ[4]=GEA_PER_VJ;
	static const double	jL[33]=JUP_PER_VL,je[39]=JUP_PER_VE;
	static const double	jvk[22]=JUP_PER_VVK,ja[12]=JUP_PER_VA;
	static const double	sL[32]=SAT_PER_VL,se[68]=SAT_PER_VE;
	static const double	svk[27]=SAT_PER_VVK,sa[38]=SAT_PER_VA;
	static const double	sb[6]=SAT_PER_VB;
	static const double	uL[7]=URA_PER_VL,ue[5]=URA_PER_VE;
	static const double	uvk[4]=URA_PER_VVK,ur[13]=URA_PER_VR;
	static const double	ul[15]=URA_PER_VLH,ub[7]=URA_PER_VB;
	static const double	nL[5]=NEP_PER_VL,ne[5]=NEP_PER_VE,na[4]=NEP_PER_VA;
	static const double	nvk[4]=NEP_PER_VVK,nr[6]=NEP_PER_VR;
	static const double	nl[5]=NEP_PER_VLH,nb[2]=NEP_PER_VB;
	static double	j[14];
	static double	u[41];
	double	tmp;


	if(pA->astreid==JUP)	// premier passage = initialisation
	{
		j[0]=t/5.0+.1;		// j1
		j[1]=VJ[0].DotT(t);	// j2
		j[2]=VJ[1].DotT(t);
		j[3]=VJ[2].DotT(t);	// j4

		j[4]=-2*j[1]+5*j[2];	// j5
		j[5]=2*j[1]-6*j[2]+3*j[3];	// j6

		j[6]=j[2]-j[1];		// j7

		j[7]=VJ[3].DotT(t);	// j8 Uranus et Neptune
		j[8]=2*j[7]-j[3];	// j9 Ura Nept

		j[9]=j[3]-j[1];		// ja Uranus
		j[10]=j[3]-j[2];	// jb Uranus, j8 Saturne

		j[11]=j[7]-j[3];	// jc Uranus et Neptune
					
		j[12]=j[7]-j[1];	// ja Neptune
		j[13]=j[7]-j[2];	// jb Neptune

		//---------------------------------------------

		u[0]=sin(j[2]);		// u1
		u[1]=cos(j[2]);

		tmp=j[2]+j[2];
		u[2]=sin(tmp);		// u3
		u[3]=cos(tmp);

		u[4]=sin(j[4]);		// u5
		u[5]=cos(j[4]);
		u[6]=sin(j[4]+j[4]);	// u7
		u[7]=sin(j[5]);
		u[8]=sin(j[6]);		// u9
		u[9]=cos(j[6]);		// ua
		
		tmp=j[6]+j[6];
		u[10]=sin(tmp);		// ub
		u[11]=cos(tmp);		// uc

		tmp+=j[6];
		u[12]=sin(tmp);		// ud
		u[13]=cos(tmp);		// ue

		tmp+=j[6];
		u[14]=sin(tmp);		// uf
		u[15]=cos(tmp);		// ug

		u[16]=cos(tmp+j[6]);	// uh, vh

		tmp=3*j[2];
		u[17]=sin(tmp);		// ui
		u[18]=cos(tmp);		// uj

		tmp+=j[2];
		u[19]=sin(tmp);		// uk
		u[20]=cos(tmp);		// ul

		u[21]=cos(j[4]+j[4]);	// um, vi
		u[22]=sin(5*j[6]);	// un

		tmp=j[10]+j[10];
		u[23]=sin(tmp);		// uo
		u[24]=cos(tmp);		// up

		tmp+=j[10];
		u[25]=sin(tmp);		// uq
		u[26]=cos(tmp);		// ur

		u[27]=sin(j[8]);	// ut, vj
		u[28]=cos(j[8]);	// uu

		tmp=j[8]+j[8];
		u[29]=sin(tmp);		// uv
		u[30]=cos(tmp);		// uw

		u[31]=sin(j[10]);	// ux
		u[32]=cos(j[10]);	// uy
		u[33]=sin(j[3]);	// uz
		u[34]=cos(j[3]);	// va

		tmp=j[3]+j[3];
		u[35]=sin(tmp);		// vb
		u[36]=cos(tmp);		// vc

		tmp=j[11]+j[11];
		u[37]=sin(tmp);		// vd
		u[38]=cos(tmp);		// ve

		u[39]=sin(j[7]);	// vf
		u[40]=cos(j[7]);	// vg
	}

	switch(pA->astreid)
	{
		case	JUP:
			// perturbation en longitude moyenne (L)
			Lper=(jL[0]+(jL[1]+jL[2]*j[0])*j[0])*u[4];
			Lper+=(jL[3]+(jL[4]+jL[5]*j[0])*j[0])*u[5];
			Lper+=(jL[6]+(jL[7]+jL[8]*j[0])*j[0])*u[6];
			Lper+=jL[9]*u[7]+jL[10]*u[8]+jL[11]*u[10];
			Lper+=jL[12]*u[12]+jL[13]*u[14];
			tmp=jL[14]*u[10];
			tmp+=(jL[15]+jL[16]*j[0])*u[8]+jL[17]*u[12];
			tmp+=jL[20]*u[11];
			tmp+=(jL[21]+jL[22]*j[0])*u[9];
			Lper+=tmp*u[0];
			tmp=(jL[18]+jL[19]*j[0])*u[8];
			tmp+=jL[23]*u[10];
			tmp+=(jL[24]*j[0]+jL[25])*u[9]+jL[26];
			tmp+=jL[27]*u[11]+jL[28]*u[13];
			Lper+=tmp*u[1];
			Lper+=(jL[29]*u[8]+jL[30]*u[10])*u[2];
			Lper+=(jL[31]*u[9]+jL[32]*u[11])*u[3];

			// perturbation en excentricit� (e)
			eper=(je[0]+(je[1]+je[2]*j[0])*j[0])*u[4]+(je[3]+je[4]*j[0])*u[5];
			tmp=je[5]*u[8]+je[6]*u[10]+je[7]*u[12]+je[8];
			tmp+=(je[9]+je[10]*j[0])*u[9]+je[11]*u[11];
			eper+=tmp*u[0];
			tmp=(je[12]+je[13]*j[0])*u[8]+je[14]*u[10]+je[15];
			tmp+=je[16]*u[9]+je[17]*u[11]+je[18]*u[13]+je[19]*u[15];
			tmp+=je[20]*u[16];
			eper+=tmp*u[1];
			tmp=(je[21]+je[22]*j[0])*u[8]+je[23]*u[10];
			tmp+=je[24]*u[12]+(je[25]*j[0]+je[26])*u[9]+je[27]*u[11];
			tmp+=je[28]*u[13];
			eper+=tmp*u[2];
			tmp=(je[29]*j[0]+je[30])*u[8]+je[31]*u[10];
			tmp+=je[32]*u[12]+je[33]+(je[34]+je[35]*j[0])*u[9];
			tmp+=je[36]*u[11]+je[37]*u[13];
			eper+=tmp*u[3];
			eper*=je[38];

			// perturbation en longitude du p�rih�lie (omegab)
			vkper=(jvk[0]+jvk[1]*j[0])*u[4];
			vkper+=(j[0]*(jvk[3]*j[0]+jvk[4])+jvk[5])*u[5];
			tmp=jvk[2];
			tmp+=jvk[6]*u[9]+(jvk[7]+jvk[8]*j[0])*u[8];
			tmp+=jvk[9]*u[11]+jvk[10]*u[13];
			vkper+=tmp*u[0];
			vkper+=(jvk[11]*u[8]+jvk[12]*u[10]+jvk[13]*u[9])*u[1];
			vkper+=(jvk[14]*u[8]+jvk[15]*u[10]+jvk[16]*u[9]+jvk[17]*u[11])*u[2];
			vkper+=(jvk[18]*u[8]+jvk[19]*u[10]+jvk[20]*u[9]+jvk[21]*u[11])*u[3];

			// perturbation en demi grand axe (a)
			aper=ja[0]*u[9]+ja[1]*u[5]+ja[2]*u[11]+ja[3]*u[13]+ja[4]*u[15];
			aper+=(ja[5]*u[8]+ja[6]*u[11])*u[0];
			aper+=(ja[7]*u[10]+ja[8]*u[12]+ja[9]*u[9]+ja[10]*u[11])*u[1];
			aper*=ja[11];
			break;

		case	SAT:
			// perturbation en longitude moyenne (L)
			Lper=sL[0]*u[6]+sL[1]*u[7]+sL[2]*u[8];
			Lper+=(sL[3]+(sL[4]+sL[5]*j[0])*j[0])*u[4];
			Lper+=(sL[6]+(sL[7]+sL[8]*j[0])*j[0])*u[5];
			Lper+=sL[9]*u[12]+sL[10]*u[14]+sL[11]*u[0];
			Lper+=sL[13]*u[10];
			tmp=sL[12]*u[10];
			tmp+=(sL[14]+sL[15]*j[0])*u[8]+sL[16]*u[12];
			tmp+=(sL[17]+sL[18]*j[0])*u[9]+sL[19]*u[11];
			Lper+=tmp*u[0];
			tmp=(sL[20]+sL[21]*j[0])*u[8]+sL[22]*u[11];
			tmp+=(sL[23]+sL[24]*j[0])*u[9]+sL[25]*u[13];
			Lper+=tmp*u[1];
			Lper+=(sL[26]*u[8]+sL[27]*u[10]+sL[28]*u[25])*u[2];
			Lper+=(sL[29]*u[9]+sL[30]*u[11]+sL[31]*u[26])*u[3];

			// perturbation en excentricit� (e)
			eper=(se[0]+(se[1]+se[2]*j[0])*j[0])*u[4];
			eper+=(se[3]+(se[4]+se[5]*j[0])*j[0])*u[5]+(se[6]+se[7]*j[0])*u[6];
			eper+=(se[8]+se[9]*j[0])*u[21]+se[10]*u[10];
			tmp=se[11]+(se[12]+se[13]*j[0])*u[8]+(se[14]+se[15]*j[0])*u[10];
			tmp+=se[16]*u[9]+se[17]*u[11]+se[18]*u[13]+se[19]*u[15];
			tmp+=se[20]*u[16]+se[21]*u[24];
			eper+=tmp*u[0];
			tmp=se[22]+se[23]*j[0];
			tmp+=se[24]*u[8]+se[25]*u[10]+se[26]*u[12]+se[27]*u[14];
			tmp+=se[28]*u[22]+(se[29]+se[30]*j[0])*u[9];
			tmp+=(se[31]+se[32]*j[0])*u[11]+se[33]*u[23];
			eper+=tmp*u[1];
			tmp=se[34]+(se[35]+se[36]*j[0])*u[8]+se[37]*u[10]+se[38]*u[12];
			tmp+=se[39]*u[14]+(se[40]+se[41]*j[0])*u[9];
			tmp+=(se[42]+se[43]*j[0])*u[11]+se[44]*u[13]+se[45]*u[25];
			tmp+=se[46]*u[26];
			eper+=tmp*u[2];
			tmp=se[47]+(se[48]+se[49]*j[0])*u[8];
			tmp+=(se[50]+se[51]*j[0])*u[10]+se[52]*u[12];
			tmp+=(se[53]+se[54]*j[0])*u[9]+(se[55]+se[56]*j[0])*u[11];
			tmp+=se[57]*u[13]+se[58]*u[15]+se[59]*u[25]+se[60]*u[26];
			eper+=tmp*u[3];
			eper+=(se[61]*u[8]+se[62]*u[12])*u[17];
			eper+=(se[63]*u[9]+se[64]*u[13])*u[18];
			eper+=se[65]*u[13]*u[19]+se[66]*u[12]*u[20];
			eper*=se[67];
			
			// perturbation en longitude du p�rih�lie (omegab)
			vkper=(svk[0]+(svk[1]+svk[2]*j[0])*j[0])*u[4];
			vkper+=svk[3]*u[8];
			vkper+=(svk[4]+(svk[5]+svk[6]*j[0])*j[0])*u[5];
			vkper+=(svk[8]*u[8]+svk[9]*u[10]+svk[10]*u[12])*u[0];
			vkper+=(svk[7]+svk[11]*u[9]+svk[12]*u[11]+svk[13]*u[13])*u[1];
			tmp=(svk[14]+svk[15]*j[0])*u[8];
			tmp+=(svk[17]+svk[18]*j[0])*u[9];
			tmp+=(svk[19]+svk[20]*j[0])*u[11];
			vkper+=tmp*u[2];
			tmp=svk[16]*u[10]+(svk[21]+svk[22]*j[0])*u[8];
			tmp+=(svk[23]+svk[24]*j[0])*u[9];
			tmp+=(svk[25]+svk[26]*j[0])*u[11];
			vkper+=tmp*u[3];

			// perturbation en demi grand axe (a)
			aper=sa[0]*u[4]+sa[2]*u[5]+sa[4]*u[9]+sa[6]*u[11];
			aper+=sa[8]*u[13]+sa[11]*u[15]+sa[13]*u[16];
			tmp=sa[15]+sa[17]*u[8]+sa[19]*u[10]+sa[21]*u[12];
			tmp+=sa[23]*u[14]+sa[25]*u[9]+sa[27]*u[11];
			tmp+=sa[29]*u[13]+sa[31]*u[15];
			aper+=tmp*u[0];
			tmp=sa[1]*u[10]+sa[3]*u[12]+sa[5]*u[14]+sa[7]*u[9];
			tmp+=(sa[9]+sa[10]*j[0])*u[11]+sa[12]*u[13];
			tmp+=sa[33]+sa[35]*u[8];
			aper+=tmp*u[1];
			aper+=(sa[14]*u[10]+sa[16]*u[9]+sa[18]*u[11]+sa[20]*u[13])*u[2];
			aper+=(sa[22]*u[8]+sa[24]*u[10]+sa[26]*u[11]+sa[28]*u[13])*u[3];
			aper+=(sa[30]*u[8]+sa[32]*u[12])*u[17];
			aper+=(sa[34]*u[9]+sa[36]*u[13])*u[18];
			aper*=sa[37];

			// perturbation en latitude h�liocentrique (b)
			bper=(sb[0]*u[0]+sb[1]*u[1])*u[9];
			bper+=(sb[2]*u[10]+sb[3]*u[11])*u[2];
			bper+=(sb[4]*u[10]+sb[5]*u[11])*u[3];
			break;

		case	URA:
			// perturbation en longitude moyenne (L)
			Lper=(uL[0]+uL[1]*j[0])*u[27];
			Lper+=(uL[2]+uL[3]*j[0])*u[28]+uL[4]*u[29];
			Lper+=uL[5]*u[30]+uL[6]*sin(j[5]);

			// perturbation en longitude du p�rih�lie (omegab)
			vkper=uvk[0]*u[27]+uvk[1]*u[29];
			vkper+=(uvk[2]+uvk[3]*j[0])*u[28];

			// perturbation en excentricit� (e)
			eper=(ue[0]*j[0]+ue[1])*u[27]+ue[2]*u[28]+ue[3]*u[30];
			eper*=ue[4];

			// perturbation en demi grand axe (a)
			aper=URA_PER_A*u[28];

			// perturbation en longitude h�liocentrique (l)
			lper=(ul[0]+(ul[1]+ul[2]*j[0])*j[0])*cos(j[3]+j[10]);
			lper+=(ul[3]+ul[4]*j[0])*sin(j[3]+j[10]);
			lper+=(ul[5]+(ul[6]+ul[7]*j[0])*j[0])*cos(2*j[3]+j[10]);
			lper+=ul[8]*sin(j[3]+3*j[11])+ul[9]*sin(j[9]);
			lper+=ul[10]*u[31]+ul[11]*u[32];
			lper+=ul[12]*sin(j[11])+ul[13]*u[37];
			lper+=ul[14]*sin(3*j[11]);

			// perturbation en latitude h�liocentrique (b)
			bper=(ub[0]*u[31]+ub[1]*u[32]+ub[2]*cos(4*j[11]))*u[33];
			bper+=(ub[3]*u[31]+ub[4]*u[32]+ub[5]*sin(4*j[10]))*u[34];
			bper+=ub[6]*(u[38]*u[35]+u[37]*u[36]);
			
			// perturbation en rayon vecteur (r)
			rper=ur[0]+ur[1]*cos(j[9])+ur[2]*u[34]+ur[3]*u[32];
			rper+=ur[4]*u[38]+ur[5]*(cos(j[11])-cos(3*j[11]));
			rper+=(ur[6]*u[34]+ur[7]*u[33]+ur[8]*u[36])*u[31];
			rper+=(ur[9]*u[34]+ur[10]*u[33]+ur[11]*u[35])*u[32];
			rper*=ur[12];
			break;

		case	NEP:
			// perturbation en longitude moyenne (L)
			Lper=(nL[0]*j[0]+nL[1])*u[27];
			Lper+=(nL[2]*j[0]+nL[3])*u[28]+nL[4]*u[29];

			// perturbation en longitude du p�rih�lie (omegab)
			vkper=nvk[0]*u[27]+nvk[1]*u[28]+nvk[2]*u[29];
			vkper+=nvk[3]*u[30];

			// perturbation en excentricit� (e)
			eper=ne[0]*u[27]+ne[1]*u[29]+ne[2]*u[28]+ne[3]*u[30];
			eper*=ne[4];

			// perturbation en demi grand axe (a)
			aper=na[0]*u[28]+na[1]*u[27]+na[2]*u[30];
			aper*=na[3];

			// perturbation en longitude h�liocentrique (l)
			lper=nl[0]*sin(j[12])+nl[1]*sin(j[13]);
			lper+=nl[2]*u[37]+nl[3]*u[38]*u[39]+nl[4]*u[37]*u[40];

			// perturbation en latitude h�liocentrique (b)
			bper=nb[0]*u[38]*u[39]+nb[1]*u[37]*u[40];

			// perturbation en rayon vecteur (r)
			rper=nr[0]+nr[1]*cos(j[12])+nr[2]*cos(j[13]);
			rper+=nr[3]*cos(j[11])+nr[4]*u[38];
			rper*=nr[5];
			break;
	}
}


// formatage des r�sultats
void	Orbite::Affiche(void) const
{
	cout << "\n" << _("El�ments Orbitaux");
/*	printf("\n%s\t\t a = %9.5f",_("Demi Grand Axe"),a);
	printf("\t%s\t\t e = %9.5f",_("Excentricit�"),e);
	printf("\n%s\t i = %9.5f",_("Inclinaison / Eclipt."),i);

	printf("\t%s\t L = %9.5f",_("Longitude Moyenne"),L);
	printf("\n%s\t O = %9.5f",_("Longitude Noeud Ascend."),Omega);
	printf("\t%s\twb = %9.5f",_("Longitude P�rih�lie"),omegab);


	printf("\n%s\t M = %9.5f",_("Anomalie Moyenne"),M);
	printf("\t%s\t w = %9.5f",_("Argument P�rih�lie"),omega);

	printf("\n%s\t E = %9.5f",_("Anomalie Excentrique"),E);
	
	printf("\n%s\t\t r = %9.5f",_("Rayon Vecteur"),r);
	printf("\t%s\t\t v = %9.5f\n",_("Anomalie Vraie"),v);
*/
}

//////////////////////////////////////////////////////////////////////////////
void	OrbiteLune::Init(AstreLune *pAst)
{
	static const VectData	Vlun[5]={LUN_VO,LUN_VL,LUN_VM,LUN_VF,LUN_VD};
	static const VectData	Vlunper[4]={LUN_PER_VL,LUN_PER_VM,LUN_PER_VF,LUN_PER_VD};
	static const VectData	Vv[4]={LUN_PER_VA,LUN_PER_VB,LUN_PER_VC,LUN_PER_VE};


	pA=pAst;

	// vecteurs de d�finition des �l�ments de l'orbite
	VO=Vlun[0];
	VL=Vlun[1];
	VM=Vlun[2];
	VF=Vlun[3];
	VD=Vlun[4];

	// pour plus de pr�cision
	VLper=Vlunper[0];
	VMper=Vlunper[1];
	VFper=Vlunper[2];
	VDper=Vlunper[3];
	Va=Vv[0];
	Vb=Vv[1];
	Vc=Vv[2];
	Ve=Vv[3];
}


void    OrbiteLune::Elements(const Instant &Inst)
{
	// calcul des �l�ments de l'orbite
	Omega=VO.DotTprecis2(Inst.t,Inst.datejulienne);
	OmegaR=DEG2RAD(Omega);
	L=VL.DotTprecis2(Inst.t,Inst.datejulienne);
	L=K360(L);
	M=VM.DotTprecis2(Inst.t,Inst.datejulienne);
	M=K360(M);
	F=VF.DotTprecis2(Inst.t,Inst.datejulienne);
	F=K360(F);
	D=VD.DotTprecis2(Inst.t,Inst.datejulienne);
	D=K360(D);

	CCR=OmegaR+DEG2RAD(Vc.DotT(Inst.t));
	EE=Ve.DotT(Inst.t);
}
