/* Jean Philippe EIMER - 1999 *
 *  phil.eimer@9online.fr     */

#include <math.h>
#include <ostream>

#include <stdio.h>

#include "include.h"


// constructeur par défaut
Vecteur::Vecteur()
{
	Spher[COORD_R]=0.0;
	Spher[COORD_TETA]=0.0;
	Spher[COORD_PHI]=0.0;
	Rect[COORD_X]=0.0;
	Rect[COORD_Y]=0.0;
	Rect[COORD_Z]=0.0;
}


// constructeur avec initialisation
Vecteur::Vecteur(double a, double b, double c, bool rect=false)
{
	if(rect)
	{
		Rect[COORD_X]=a;
		Rect[COORD_Y]=b;
		Rect[COORD_Z]=c;
		Spherique();	// mise à jour des coordonnées sphériques
	}
	else
	{
		Spher[COORD_R]=a;
		Spher[COORD_TETA]=b;
		Spher[COORD_PHI]=c;
		Rectangulaire();	// mise à jour des coordonnées rectangulaires
	}
}


// conversion sphérique --> rectangulaire
Vecteur	Vecteur::Rectangulaire(void)
{
	double t=DEG2RAD(Spher[COORD_TETA]),p=DEG2RAD(Spher[COORD_PHI]);
	double cp=cos(p);

	Rect[COORD_X]=Spher[COORD_R]*cos(t)*cp;		// formules
	Rect[COORD_Y]=Spher[COORD_R]*sin(t)*cp;
	Rect[COORD_Z]=Spher[COORD_R]*sin(p);

	Spher[COORD_TETA]=K360(Spher[COORD_TETA]);	// centrage des angles
	Spher[COORD_PHI]=K360(Spher[COORD_PHI]);
	Spher[COORD_PHI]=KNEG(Spher[COORD_PHI]);

	return(*this);
}


// conversion rectangulaire --> sphérique
Vecteur	Vecteur::Spherique(void)
{
	double	p=Rect[COORD_X]*Rect[COORD_X]+Rect[COORD_Y]*Rect[COORD_Y];

	Spher[COORD_R]=sqrt(p+Rect[COORD_Z]*Rect[COORD_Z]);		// formule rayon vecteur

	Spher[COORD_TETA]=RAD2DEG(atan(Rect[COORD_Y]/Rect[COORD_X]));	// formule angle 1

	if(Rect[COORD_X]<0)
		Spher[COORD_TETA]+=180;	// lever d'indétermination suite à atan

	Spher[COORD_PHI]=RAD2DEG(atan(Rect[COORD_Z]/sqrt(p)));	// formule angle 2

	Spher[COORD_TETA]=K360(Spher[COORD_TETA]);	// centrage des angles
	Spher[COORD_PHI]=K360(Spher[COORD_PHI]);
	Spher[COORD_PHI]=KNEG(Spher[COORD_PHI]);

	return(*this);
}


// rotation autour de Ox
Vecteur	Vecteur::RotationX(double angle)
{
	double	t,ar=DEG2RAD(angle),ca=cos(ar),sa=sin(ar);

	t=Rect[COORD_Y]*ca+Rect[COORD_Z]*sa;		// formules
	Rect[COORD_Z]=Rect[COORD_Z]*ca-Rect[COORD_Y]*sa;
	Rect[COORD_Y]=t;

	Spherique();		// mise à jour des coordonnées sphériques

	return(*this);
}


// rotation autour de Oy
Vecteur	Vecteur::RotationY(double angle)
{
	double	t,ar=DEG2RAD(angle),ca=cos(ar),sa=sin(ar);

	t=Rect[COORD_X]*ca+Rect[COORD_Z]*sa;		// formules
	Rect[COORD_Z]=Rect[COORD_Z]*ca-Rect[COORD_X]*sa;
	Rect[COORD_X]=t;

	Spherique();		//mise à jour des coordonnées sphériques

	return(*this);
}


// surcharge +
Vecteur	Vecteur::operator +(const Vecteur &V) const
{
	Vecteur	S;


	S=*this;		// initialisation de la somme

	S.Rect[COORD_X]+=V.Rect[COORD_X];		// somme vectorielle
	S.Rect[COORD_Y]+=V.Rect[COORD_Y];
	S.Rect[COORD_Z]+=V.Rect[COORD_Z];

	S.Spherique();		// mise à jour des coordonnées sphériques
						// de la somme
	return(S);
}


// surcharge -
Vecteur	Vecteur::operator -(const Vecteur &V) const
{
	Vecteur	S;


	S=*this;		// initialisation de la somme

	S.Rect[COORD_X]-=V.Rect[COORD_X];		// soustraction vectorielle
	S.Rect[COORD_Y]-=V.Rect[COORD_Y];
	S.Rect[COORD_Z]-=V.Rect[COORD_Z];

	S.Spherique();		// mise à jour des coordonnées sphériques

	return(S);
}


// surcharge =
Vecteur	&Vecteur::operator =(const Vecteur &V)
{
	Spher[COORD_R]=V.Spher[COORD_R];		// affectation coordonnées sphériques
	Spher[COORD_TETA]=V.Spher[COORD_TETA];
	Spher[COORD_PHI]=V.Spher[COORD_PHI];

	Rect[COORD_X]=V.Rect[COORD_X];			// affectation coordonnées rectangulaires
	Rect[COORD_Y]=V.Rect[COORD_Y];
	Rect[COORD_Z]=V.Rect[COORD_Z];

	return(*this);
}

/*
// formatage des résultats
void	Vecteur::Affiche(unsigned short int coord=255,bool heure=false,bool hms=HMS) const
{
	char u0[3];
	char u1[4]={'�',' ',' ','\0'},u2[4]={'�',' ',' ','\0'};
	char d[2]={'d','\0'};
	char n1[7],n2[7];
	double th=Spher[COORD_TETA],ph=Spher[COORD_PHI];


	strcpy(u0,"ua");


	if(heure)
	{
		u1[0]='h';
		th/=KH2D;
	}


	if(hms)
	{
		th=Deci2Sexa(th);
		ph=Deci2Sexa(ph);

		if(heure)
		{
			u1[1]='m';
			u1[2]='s';
		}
		else
		{
			u1[1]='\'';	// '
			u1[2]='"';
		}

		u2[1]='\'';		// '
		u2[2]='"';
	}

	switch(coord)
	{
		case	HELIO:
			d[0]='r';
			strcpy(n1,"     l");
			strcpy(n2,"     b");
			break;

		case	GEO:
			strcpy(n1,"lambda");
			strcpy(n2,"  beta");
			break;

		case	GEOAPP:
			strcpy(n1,"la App");
			strcpy(n2,"be App");
			break;

		case	EQUA:
			strcpy(n1," alpha");
			strcpy(n2," delta");
			break;

		case	HORAIRE:
			strcpy(n1,"     H");
			strcpy(n2," delta");
			break;

		case	HORIZON:
			strcpy(n1,"     A");
			strcpy(n2,"     h");
			break;

		default:
			d[0]='r';
			strcpy(n1,"  teta");
			strcpy(n2,"   phi");
	}

	printf("\n     %s  %s =  %8.5f\tx = % 8.4f",d,u0,Spher[COORD_R],Rect[COORD_X]);
	printf("\n%s %s = %9.5f\ty = % 8.4f",n1,u1,th,Rect[COORD_Y]);
	printf("\n%s %s = % 9.5f\tz = % 8.4f\n",n2,u2,ph,Rect[COORD_Z]);
}
*/
//////////////////////////////////////////////////////////////////////////////
double	VectData::Dot(const VectData &V) const
{
	return(d[0]*V.d[0]+d[1]*V.d[1]+d[2]*V.d[2]+d[3]*V.d[3]);
}


double	VectData::Dot5(const double V[5]) const
{
	return(d[0]*V[0]+d[1]*V[1]+d[2]*V[2]+d[3]*V[3]+d[4]*V[4]);
}


double	VectData::DotT(double t) const
{
	return(d[0]+(d[1]+(d[2]+d[3]*t)*t)*t);
}


double	VectData::DotTprecis(double t) const
{
	double	dot,dint;

	dot=modf(d[1]*t,&dint)*360.0;
	dot+=d[0]+(d[2]+d[3]*t)*t*t;

	return(dot);
}


double	VectData::DotTprecis2(double t,double dj) const
{
	double	dot,dint;

	dot=modf((dj-DJ0)/d[1],&dint)*360.0;
	dot+=d[0]+(d[2]+d[3]*t)*t*t;

	return(dot);
}


VectData	VectData::operator =(const double Vd[4])
{
	d[0]=Vd[0];
	d[1]=Vd[1];
	d[2]=Vd[2];
	d[3]=Vd[3];
	d[4]=0.0;

	return(*this);
}


void	VectData::Affiche(void)	const
{
	printf("\n{%.5f,%.5f,%.5f,%.5f}",d[0],d[1],d[2],d[3]);

	if(d[4]!=0.0)
		printf("\t%.5f",d[4]);

}

//////////////////////////////////////////////////////////////////////////////
// conversion D�cimal --> Sexag�simal
//double	Deci2Sexa(double angle, int *ph=NULL, int *pm=NULL, int *ps=NULL)
double	Deci2Sexa(double angle, int *ph, int *pm, int *ps)
{
	double	fh,fm,fs;

	fm=modf(fabs(angle),&fh)*60.0;	// extraction des minutes
	fs=modf(fm,&fm)*60.0;		// extraction des secondes (/10000)
					// minutes en nombre entier

	if(fs>59.9)			// gestion de l'arrondi des secondes
	{
		fs=0.0;
		fm+=1.0;
	}

	if(fm==60.0)			// gestion de l'arrondi des minutes
	{
		fm=0.0;
		fh++;
	}

	if(ph!=NULL)
	{
		*ph=(int)fh;

		if(angle < 0.0)
			*ph=-*ph;
	}

	if(pm!=NULL)
		*pm=(unsigned short int)fm;

	if(ps!=NULL)
		*ps=(unsigned short int)fs;

	fm*=0.01;			// minutes /100
	fs*=0.0001;			// secondes /10000

	return(copysign(fh+fm+fs,angle));	// signe r�tabli
}


