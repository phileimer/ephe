/* Jean Philippe EIMER - 1999 *
 *  phil.eimer@9online.fr     */

#include <math.h>

#include <iostream>
#include <iomanip>
using namespace std;

#include "include.h"


void	Astre::Init(Ephemerides *pEphe, unsigned int id)
{
	// coordonnées géocentriques pour le Soleil et la Lune
	// éléments orbitaux des planètes


	static const wstring	Vnom[11] = {	MER_NOM, VEN_NOM, \
									MAR_NOM, JUP_NOM, SAT_NOM, \
									URA_NOM, NEP_NOM, PLU_NOM, \
									TER_NOM, SOL_NOM, LUN_NOM};


	pE = pEphe;
	astreid = id;
	nom = Vnom[id];

	ElemOrbit.Init(this);

	if( id == TER )
		Coord[HELIO] = ElemOrbit2EcliptHelio(ElemOrbit);

	Parallaxe();
}


Vecteur Astre::EcliptGeo(void)
{
	Vecteur	Vout;


	ElemOrbit.Elements(pE->Observ.Inst.t);
	ElemOrbit.Perturbations(pE->Observ.Inst.t);
	ElemOrbit.Anomalies();

	Vout=ElemOrbit2EcliptHelio(ElemOrbit);

	if(astreid!=SOL)
	{
		Coord[HELIO]=Vout;
		Vout=EcliptHelio2EcliptGeo(Vout,pE->Soleil.Coord[GEO]);
	}

	Coord[GEO]=Vout;

	return(Vout);
}


Vecteur	Astre::EcliptGeo2Horizon(double t)
{
	Coord[GEOAPP]=EcliptGeo2EcliptGeoApp(Coord[GEO],pE->Soleil.Coord[GEO].Spher[COORD_TETA],t);
	Coord[EQUA]=EcliptGeoApp2Equa(Coord[GEOAPP]);
	Coord[HORAIRE]=Equa2Horaire(Coord[EQUA]);
	Coord[HORIZON]=Horaire2Horizon(Coord[HORAIRE]);

	return(Coord[HORIZON]);
}


Vecteur Astre::ElemOrbit2EcliptHelio(const Orbite &Orb) const
{
	Vecteur	Vout;


	Vout.Spher[COORD_R]=Orb.r;
	Vout.Spher[COORD_TETA]=Orb.omega+Orb.v;
	Vout.Spher[COORD_PHI]=0;

	Vout.Rectangulaire();

	Vout.RotationX(-Orb.i);

	Vout.Spher[COORD_TETA]+=Orb.Omega;

	Vout.Spher[COORD_R]+=Orb.rper;		// Corrections
	Vout.Spher[COORD_TETA]+=Orb.lper;
	Vout.Spher[COORD_PHI]+=Orb.bper;

	Vout.Rectangulaire();

	return(Vout);
}


Vecteur Astre::EcliptHelio2EcliptGeo(const Vecteur &Vin,const Vecteur &GeoSol) const
{
	Vecteur V1=Vin,V2=GeoSol,Vout;

	V1.Spher[COORD_TETA]-=GeoSol.Spher[COORD_TETA];		// l-lamdaS
	V1.Rectangulaire();

	V2.Spher[COORD_TETA]=0;
	V2.Spher[COORD_PHI]=0;
	V2.Rectangulaire();

	Vout=V1+V2;

	Vout.Spher[COORD_TETA]+=GeoSol.Spher[COORD_TETA];	// l+lamdaS = lamda
	Vout.Rectangulaire();

	return(Vout);		// d, lamda et beta
}


Vecteur Astre::EcliptGeo2EcliptHelio(const Vecteur &Vin,const Vecteur &GeoSol) const
{
	Vecteur V1=Vin,V2=GeoSol,Vout;

	V1.Spher[COORD_TETA]-=GeoSol.Spher[COORD_TETA];		// lamda-lamdaS
	V1.Rectangulaire();


	V2.Spher[COORD_TETA]=0;
	V2.Spher[COORD_PHI]=0;
	V2.Rectangulaire();

	Vout=V1-V2;

	Vout.Spher[COORD_TETA]+=GeoSol.Spher[COORD_TETA];	// lamda+lamdaS = l
	Vout.Rectangulaire();

	return(Vout);		// r, l et b
}


Vecteur Astre::EcliptGeo2EcliptGeoApp(const Vecteur &Vin,double lamdaS,double t) const
{
	Orbite	Orb=ElemOrbit;
	Vecteur	Vout;


	if(astreid <= pE->maxastr)
	{
		// trajet de la lumière
		double	tl=TempsLumiere(Vin.Spher[COORD_R]);
		double	dL=tl*Orb.deltaL;
		Orb.L-=dL;		// modification de la longitude moyenne
		Orb.MRsp-=DEG2RAD(dL);	// modification de l'anomalie moyenne

		if(astreid<JUP)
		{
			Orb.L-=Orb.Lper;
			Orb.Perturbations(t-tl);
		}
		else
			Orb.Lper=Orb.vkper=Orb.eper=Orb.aper=0;

		Orb.Anomalies();

		Vout=EcliptHelio2EcliptGeo(ElemOrbit2EcliptHelio(Orb),pE->Soleil.Coord[GEO]);
	}
	else
		Vout=Vin;

	Vout.Spher[COORD_TETA]+=pE->Nutat.lamda;

	if(astreid!=LUN)
	{
		double	lsl=DEG2RAD(lamdaS-Vin.Spher[COORD_TETA]);
		Vout.Spher[COORD_TETA]-=COEF_APP*cos(lsl)/cos(DEG2RAD(Vin.Spher[COORD_PHI]));
		Vout.Spher[COORD_PHI]-=COEF_APP*sin(lsl)*sin(DEG2RAD(Vin.Spher[COORD_PHI]));
	}

	Vout.Rectangulaire();

	return(Vout);			// lamda App et beta App
}

#if 0
Vecteur Astre::EcliptGeoApp2EcliptGeo(const Vecteur &Vin,double lamdaS) const
{
	Vecteur	Vout;
	double	lsl=DEG2RAD(lamdaS-Vin.Spher[COORD_TETA]);


	Vout=Vin;
	Vout.Spher[COORD_TETA]-=pE->Nutat.lamda;

	if(astreid!=LUN)
	{
		Vout.Spher[COORD_TETA]+=COEF_APP*cos(lsl)/cos(Vin.Spher[COORD_PHI]);
		Vout.Spher[COORD_PHI]+=COEF_APP*sin(lsl)*sin(Vin.Spher[COORD_PHI]);
	}

	Vout.Rectangulaire();

	return(Vout);		// lamda et beta
}
#endif

Vecteur Astre::EcliptGeoApp2Equa(const Vecteur &Vin) const
{
	Vecteur Vout;


	Vout=Vin;
	Vout.RotationX(-(pE->epsilon+pE->Nutat.epsilon));

	return(Vout);		// alpha et delta
}

#if 0
Vecteur Astre::Equa2EcliptGeoApp(const Vecteur &Vin) const
{
	Vecteur	Vout;


	Vout=Vin;
	Vout.RotationX(pE->epsilon+pE->Nutat.epsilon);

	return(Vout);		// lamda et beta
}
#endif

Vecteur Astre::Equa2Horaire(const Vecteur &Vin) const
{
	Vecteur	Vout;


	Vout=Vin;

	Vout.Spher[COORD_TETA]=pE->Observ.Inst.ts*KH2D-Vin.Spher[COORD_TETA];	// TS-alpha
	Vout.Rectangulaire();

	return(Vout);		// H et delta
}

#if 0
Vecteur Astre::Horaire2Equa(const Vecteur &Vin) const
{
	Vecteur	Vout;


	Vout=Vin;

	Vout.Spher[COORD_TETA]=pE->Observ.Inst.ts*KH2D-Vin.Spher[COORD_TETA];		// TS-H
	Vout.Rectangulaire();

	return(Vout);			// alpha et delta
}
#endif

Vecteur Astre::Horaire2Horizon(const Vecteur &Vin) const
{
	Vecteur	Vout;


	Vout=Vin;
	Vout.RotationY(pE->Observ.latitude-90.0);

	return(Vout);			// A et h, azimuth et hauteur
}

#if 0
Vecteur Astre::Horizon2Horaire(const Vecteur &Vin) const
{
	Vecteur	Vout;


	Vout=Vin;
	Vout.RotationY(90.0-pE->Observ.latitude);

	return(Vout);			// H et delta
}
#endif

// Calcul du temps de trajet de la lumière
double	Astre::TempsLumiere(double dist) const
{
	return(dist*COEF_LUM);
}


// calcul des levers et couchers
void	Astre::LeverCoucher(void)
{
	Vecteur	Eq12;
	double	h0,tu,tu0,trig[3];


	h0=parallaxe-diametre/7200.0-pE->Observ.refraction-pE->Observ.eta1;

	trig[2]=DEG2RAD(pE->Observ.latitude);
	trig[1]=sin(trig[2]);
	trig[2]=cos(trig[2]);

	// lever
	tu0=12.0;	// départ à 12 h TU
	Eq12=tu2Eq(tu0);
	trig[0]=sin(DEG2RAD(h0+pE->Observ.eta2L));
	tu=Eq2tu(Eq12,trig,LEVER);

	while(TUcond(&tu,tu0,PRECISION_LC))
	{
		tu0=tu;
		tu=Eq2tu(tu2Eq(tu),trig,LEVER);
	}

	Lever.tu2t(pE->Observ.Inst,tu,pE->Observ.zone);

	// coucher
	tu0=12.0;	// départ à 12 h TU
	trig[0]=sin(DEG2RAD(h0+pE->Observ.eta2C));
	tu=Eq2tu(Eq12,trig,COUCHER);

	while(TUcond(&tu,tu0,PRECISION_LC))
	{
		tu0=tu;
		tu=Eq2tu(tu2Eq(tu),trig,COUCHER);
	}

	Coucher.tu2t(pE->Observ.Inst,tu,pE->Observ.zone);

	// passage au méridien
	tu0=12.0;	// départ à 12 h TU
	tu=Eq2tu(Eq12,trig,MERIDIEN);

	while(TUcond(&tu,tu0,PRECISION_M))
	{
		tu0=tu;
		tu=Eq2tu(tu2Eq(tu),trig,MERIDIEN);
	}

	Meridien.tu2t(pE->Observ.Inst,tu,pE->Observ.zone);
}


// calcul du temps universel à partir de la hauteur et des coordonnées équatoriales
double	Astre::Eq2tu(const Vecteur &Eq, double trig[], unsigned int lcm)
{
	double	H,tu=0;


	if(lcm==MERIDIEN)
		H=0;
	else
	{
		H=DEG2RAD(Eq.Spher[COORD_PHI]);
		H=(trig[0]-trig[1]*sin(H))/trig[2]/cos(H);

		if(H>=-1.0 && H<=1.0)
			H=RAD2DEG(acos(H));
		else
			tu=TUSING;

		if(lcm==LEVER)
			H=-H;
	}

	if(tu!=TUSING)
	{
		tu=(H+pE->Observ.longitude+Eq.Spher[COORD_TETA])/KH2D-pE->Observ.Inst.tsg0;
		tu=K24(tu);
		tu=tu/SOL2SID;
	}

	return(tu);
}


// calcul des coordonnées équatoriales à partir du temps universel
Vecteur	Astre::tu2Eq(double tu)
{
	Orbite	Orb=ElemOrbit;
	OrbiteLune	OrbLune;
	Vecteur	Geo;
	Instant	I;


	I.tu2t(pE->Observ.Inst,tu);

	if(astreid==LUN)
	{
		OrbLune=pE->Lune.ElemOrbitLune;
		OrbLune.Elements(I);
		Geo=pE->Lune.ElemOrbit2EcliptGeo(OrbLune,I.t);
	}
	else
	{
		Orb.Elements(I.t);
		Orb.Perturbations(I.t);
		Orb.Anomalies();

		Geo=ElemOrbit2EcliptHelio(Orb);

		if(astreid!=SOL)
		{
			Orb=pE->Soleil.ElemOrbit;
			Orb.Elements(I.t);
			Orb.Perturbations(I.t);
			Orb.Anomalies();
			Geo=EcliptHelio2EcliptGeo(Geo,ElemOrbit2EcliptHelio(Orb));
		}
	}


	return(EcliptGeoApp2Equa(Geo));
}


bool	Astre::TUcond(double *ptu, double tu0, double prec)
{
	bool	cond;


	if(*ptu == TUSING)
		cond = false;
	else
	{
		double	abs=fabs(*ptu-tu0);

		if(abs>12.1)
		{
			*ptu=TUSING;
			cond = false;
		}
		else
			cond = abs>prec;
	}

	return(cond);
}

// calcul de la phase
double	Astre::Phase(void)
{
	phase=acos((Coord[HELIO].Spher[COORD_R]*Coord[HELIO].Spher[COORD_R]+ \
				Coord[GEOAPP].Spher[COORD_R]*Coord[GEOAPP].Spher[COORD_R]- \
				pE->Soleil.Coord[GEO].Spher[COORD_R]*pE->Soleil.Coord[GEO].Spher[COORD_R]) \
                /2/Coord[HELIO].Spher[COORD_R]/Coord[GEOAPP].Spher[COORD_R]);

	phase100=cos(phase/2.0);
	phase=RAD2DEG(phase);

	phase100*=100.0*phase100;

	return(phase);
}


// calcul de l'élongation
double	Astre::Elongation(void)
{
	elongation=acos((pE->Soleil.Coord[GEO].Spher[COORD_R]*pE->Soleil.Coord[GEO].Spher[COORD_R]+ \
                Coord[GEOAPP].Spher[COORD_R]*Coord[GEOAPP].Spher[COORD_R]- \
				Coord[HELIO].Spher[COORD_R]*Coord[HELIO].Spher[COORD_R]) \
                /2/pE->Soleil.Coord[GEO].Spher[COORD_R]/Coord[GEOAPP].Spher[COORD_R]);

	elongation=RAD2DEG(elongation);

	return(elongation);
}


// calcul du diamètre apparent
double	Astre::Diametre(void)
{
	static const double	Vdia[11]={MER_DIA, VEN_DIA, \
								MAR_DIA, JUP_DIA, SAT_DIA, \
								URA_DIA, NEP_DIA, PLU_DIA, 0, SOL_DIA, LUN_DIA};


	if(astreid == LUN)
	{
		diametre=2*asin(Vdia[LUN]*sin(DEG2RAD(parallaxe)));
		diametre=RAD2DEG(diametre)*3600.0;						// en secondes d'arc
	}
	else
		diametre=Vdia[astreid]/Coord[GEOAPP].Spher[COORD_R];	// en secondes d'arc

	return(diametre);
}


// calcul de la parallaxe
double	Astre::Parallaxe(void)
{
	parallaxe=0;

	return(parallaxe);
}


// calcul de la magnitude
double	Astre::Magnitude(void)
{
	static const VectData	Vmag[11]={MER_VMAG,VEN_VMAG, \
				MAR_VMAG,JUP_VMAG,SAT_VMAG, \
				URA_VMAG,NEP_VMAG,PLU_VMAG,{{0,0,0,0}},SOL_VMAG,LUN_VMAG};
	VectData	Vphi;

	Vphi.d[0]=1.0;
	Vphi.d[1]=fabs(phase)/100.0;
	Vphi.d[2]=Vphi.d[1]*Vphi.d[1];
	Vphi.d[3]=Vphi.d[2]*Vphi.d[1];

	magnitude=Vphi.Dot(Vmag[astreid]);

	if(astreid!=SOL)
		magnitude+=5.0*log10(Coord[GEOAPP].Spher[COORD_R]*Coord[HELIO].Spher[COORD_R]);

	return(magnitude);
}

/*
// Gestion de l'affichage des coordonnées
void	Astre::Affiche(unsigned int coord, bool hms) const
{
	cout << "\n" << nom << "\t: ";

	switch(coord)
	{
		case	HELIO:
			cout << _(L"Coordonnées Ecliptiques Héliocentriques");
			Coord[HELIO].Affiche(HELIO, false, hms);
			break;

		case	GEO:
			cout << _(L"Coordonnées Ecliptiques Géocentriques");
			Coord[GEO].Affiche(GEO, false, hms);
			break;

		case	GEOAPP:
			cout << _(L"Coordonnées Ecliptiques Géocentriques Apparentes");
			Coord[GEOAPP].Affiche(GEOAPP, false, hms);
			break;

		case	EQUA:
			cout << _(L"Coordonnées Equatoriales");
			Coord[EQUA].Affiche(EQUA, true, hms);
			break;

		case	HORAIRE:
			cout << _(L"Coordonnées Horaires");
			Coord[HORAIRE].Affiche(HORAIRE, true, hms);
			break;

		case	HORIZON:
			cout << _(L"Coordonnées Horizontales");
			Coord[HORIZON].Affiche(HORIZON, false, hms);
			break;
	}
}
*/


///////////////////////////////////////////////////////////////////////////////
void	AstreLune::Init(Ephemerides *pEphe)
{
	pE = pEphe;
	astreid = LUN;
	nom = LUN_NOM;

	ElemOrbitLune.Init(this);
}


Vecteur	AstreLune::EcliptGeo(void)
{
	ElemOrbitLune.Elements(pE->Observ.Inst);
	Coord[GEO]=ElemOrbit2EcliptGeo(ElemOrbitLune,pE->Observ.Inst.t,&parallaxe);
	Coord[HELIO]=EcliptGeo2EcliptHelio(Coord[GEO],pE->Soleil.Coord[GEO]);

	return(Coord[GEO]);
}


//Vecteur AstreLune::ElemOrbit2EcliptGeo(const OrbiteLune &Orb,double t,double *pparalx=NULL)
Vecteur AstreLune::ElemOrbit2EcliptGeo(const OrbiteLune &Orb,double t,double *pparalx)
{
	static const VectData	Vlamda[50]=LUN_L_VV,Vbeta[45]=LUN_B_VV;
	static const VectData	Vpara[31]=LUN_P_VV;
	static const double	Clamda[50]=LUN_L_VC,Cbeta[45]=LUN_B_VC;
	static const double	Cpara[31]=LUN_P_VC;
	Vecteur	Vout;
	VectData	V,S;
	unsigned short int	i;
	double	paralx;


	S.d[0]=sin(DEG2RAD(Orb.Va.DotT(t)));
	S.d[1]=sin(Orb.OmegaR);
	S.d[2]=LUN_PER_CB*sin(DEG2RAD(Orb.Vb.DotT(t)));
	S.d[3]=sin(Orb.CCR);
	S.d[4]=0.0;

	V.d[0]=DEG2RAD(Orb.M+Orb.VMper.Dot(S));
	V.d[1]=DEG2RAD(Orb.F+Orb.VFper.Dot(S));
	V.d[2]=DEG2RAD(Orb.D+Orb.VDper.Dot(S));
	V.d[3]=pE->Soleil.ElemOrbit.MRsp+DEG2RAD(LUN_PER_CMS*S.d[0]);
	V.d[4]=0.0;

	//************************************
	// calcul de la longitude géocentrique
	Vout.Spher[COORD_TETA]=0;
	for(i=45;i<50;i++)
		Vout.Spher[COORD_TETA]+=Clamda[i]*sin(V.Dot(Vlamda[i]));

	Vout.Spher[COORD_TETA]*=Orb.EE;

	for(i=27;i<45;i++)
		Vout.Spher[COORD_TETA]+=Clamda[i]*sin(V.Dot(Vlamda[i]));

	Vout.Spher[COORD_TETA]*=Orb.EE;

	for(i=0;i<27;i++)
		Vout.Spher[COORD_TETA]+=Clamda[i]*sin(V.Dot(Vlamda[i]));

	Vout.Spher[COORD_TETA]+=Orb.L+Orb.VLper.Dot(S);

	//***********************************
	// calcul de la latitude géocentrique
	Vout.Spher[COORD_PHI]=Orb.EE*Cbeta[44]*sin(V.Dot(Vbeta[44]));

	for(i=29;i<44;i++)
		Vout.Spher[COORD_PHI]+=Cbeta[i]*sin(V.Dot(Vbeta[i]));

	Vout.Spher[COORD_PHI]*=Orb.EE;

	for(i=0;i<29;i++)
		Vout.Spher[COORD_PHI]+=Cbeta[i]*sin(V.Dot(Vbeta[i]));

	Vout.Spher[COORD_PHI]*=(1.0-LUN_PER_CW1*cos(Orb.OmegaR)-LUN_PER_CW2*cos(Orb.CCR));

	//**********************
	// calcul de la parallaxe
	paralx=Orb.EE*Cpara[30]*cos(V.Dot(Vpara[30]));

	for(i=18;i<30;i++)
		paralx+=Cpara[i]*cos(V.Dot(Vpara[i]));

	paralx*=Orb.EE;
	for(i=0;i<18;i++)
		paralx+=Cpara[i]*cos(V.Dot(Vpara[i]));

	Vout.Spher[COORD_R]=TER_RT/sin(DEG2RAD(paralx))/KUA2KM;
	Vout.Rectangulaire();

	if(pparalx!=NULL)
		*pparalx=paralx;

	return(Vout);
}


double	AstreLune::Age(void)
{
	age=fabs(((Coord[GEO].Spher[COORD_TETA] > pE->Soleil.Coord[GEO].Spher[COORD_TETA]?phase:-phase)-180.0)/360.0*LUN_MOIS)/*+1.0*/;

	return(age);
}

//////////////////////////////////////////////////////////////////////////////
void	Nutation::Calcule(Ephemerides *pEphe)
{
	static const double	Vlam[13][5]=NUT_L_VV,Veps[9][5]=NUT_E_VV;
	static const VectData	Vlc1[4]=NUT_L_VC1,Vec1[2]=NUT_E_VC1;
	static const double Vlc2[9]=NUT_L_VC2,Vec2[7]=NUT_E_VC2;
	VectData	V;
	unsigned short int	i;


	pE=pEphe;

	V.d[0]=DEG2RAD(pE->Soleil.ElemOrbit.L);
	V.d[1]=DEG2RAD(pE->Soleil.ElemOrbit.M);
	V.d[2]=DEG2RAD(pE->Lune.ElemOrbitLune.L);
	V.d[3]=DEG2RAD(pE->Lune.ElemOrbitLune.M);
	V.d[4]=pE->Lune.ElemOrbitLune.OmegaR;

	// Nutation en longitude
	lamda=0;

	for(i=0;i<4;i++)
		lamda+=Vlc1[i].DotT(pE->Observ.Inst.t)*sin(V.Dot5(Vlam[i]));

	for(i=4;i<13;i++)
		lamda+=Vlc2[i-4]*sin(V.Dot5(Vlam[i]));

	// Nutation en obliquité
	epsilon=Vec1[0].DotT(pE->Observ.Inst.t)*cos(V.Dot5(Veps[0]));
	epsilon+=Vec1[1].DotT(pE->Observ.Inst.t)*cos(V.Dot5(Veps[1]));

	for(i=2;i<9;i++)
		epsilon+=Vec2[i-2]*cos(V.Dot5(Veps[i]));

	// passage des secondes d'arc en degrés
	lamda/=3600.0;
	epsilon/=3600.0;
}


void	Nutation::Affiche(void)
{
//	wprintf(L"\n%S : %.8f, %S : %.8f", _(L"Nutation en longitude"), lamda, _(L"en obliquité"), epsilon);
	wcout << endl << _(L"Nutation en longitude") << L" : " << fixed << setprecision(8) << lamda << L", " << _(L"en obliquité") << L" : " << epsilon;
}

