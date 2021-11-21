/* Jean Philippe EIMER - 1999 *
 *  phil.eimer@9online.fr     */

#include <math.h>

#include <iostream>
#include <iomanip>

#include "include.h"

#ifndef NO_CONFIGURE
	#include "../config.h"
#endif

using namespace std;

void	Ephemerides::InitAstres(void)
{
	unsigned short int i;


	maxastr=PLU;
	maxastr++;

	Obliquite();		// calcul de l'obliquité (epsilon)

	Soleil.Init(this,SOL);
	Soleil.EcliptGeo();


	Lune.Init(this);
	Lune.EcliptGeo();

	Nutat.Calcule(this);	// calcul de la nutation

	Soleil.EcliptGeo2Horizon(Observ.Inst.t);
	Soleil.Diametre();
	Soleil.Magnitude();
	Soleil.LeverCoucher();


	Lune.EcliptGeo2Horizon(Observ.Inst.t);
	Lune.Diametre();
	Lune.Phase();
	Lune.Elongation();
	Lune.Age();
	Lune.LeverCoucher();


	for(i=0;i<maxastr;i++)
		Planete[i].Init(this,i);

	for(i=0;i<maxastr;i++)
	{
		Planete[i].EcliptGeo();

		Planete[i].EcliptGeo2Horizon(Observ.Inst.t);
		Planete[i].Diametre();
		Planete[i].Phase();
		Planete[i].Elongation();
		Planete[i].Magnitude();
		Planete[i].LeverCoucher();
	}

	Planete[TER].Init(this,TER);
}


void	Ephemerides::Obliquite(void)
{
	VectData	Veps=TER_EPS;


	epsilon=Veps.DotT(Observ.Inst.t-1.0);
}


//void	Ephemerides::Affiche(bool hms=TRUE)
void	Ephemerides::Affiche(bool hms)
{
	double	eps=epsilon;
	wstring	dunit = L"  ";


	if(hms)
	{
		eps=Deci2Sexa(eps);
		dunit = L"\'\"";
	}


	wcout << _(L"Ephémérides Astronomiques");
	wcout << L" (" << PACKAGE << L" " << VERSION << L")" << endl;
	wcout << _(L"-------------------------") << endl;
	Observ.Affiche();
	wcout << endl << _(L"Obliquité") << L" : " /*<< setw(7)*/ << fixed << setprecision(4) << eps << L" °" << dunit;
	AfficheCoord(HELIO,hms);
	AfficheCoord(HELIO,hms);
	AfficheCoord(GEO,hms);
	AfficheCoord(GEOAPP,hms);
	AfficheCoord(EQUA,hms);
	AfficheCoord(HORAIRE,hms);
	AfficheCoord(HORIZON,hms);
	AfficheAPP();
	AfficheLMC();
	wcout << endl;
}


void	separation(unsigned char max)
{
	unsigned char	i;

	wcout << endl;
	wcout << setfill(L'-');

	for(i=0 ; i < max + 3 ; i++)
		wcout << setw(9) << L'-';
}


void	Ephemerides::AfficheCoord(unsigned int coord, bool hms)
{
	wstring	symbol1, symbol2, symbol3;
	wstring	dunit, hunit, unit[3];
	unsigned char	i;



	if(hms)
	{
		dunit = L"°\'\"";
		hunit = L"hms";
	}
	else
	{
		dunit = L"°  ";
		hunit = L"h  ";
	}


	unit[COORD_R] = _(L"ua");
	unit[COORD_TETA] = dunit;
	unit[COORD_PHI] = dunit;
	symbol1 = L"d   ";


	// Separation
	separation(maxastr);

	// Type de coordonnees
	wcout << endl;

	switch(coord)
	{
		case	HELIO:
			wcout << _(L"Coordonnées Ecliptiques Héliocentriques");
			symbol1 = L"r   ";
			symbol2 = L"l   ";
			symbol3 = L"b   ";
			break;

		case	GEO:
			wcout << _(L"Coordonnées Ecliptiques Géocentriques");
			symbol2 = L"lamb";
			symbol3 = L"beta";
			break;

		case	GEOAPP:
			wcout << _(L"Coordonnées Ecliptiques Géocentriques Apparentes");
			symbol2 = L"lamA";
			symbol3 = L"betA";
			break;

		case	EQUA:
			wcout << _(L"Coordonnées Equatoriales");
			unit[COORD_TETA] = hunit;
			symbol2 = L"alph";
			symbol3 = L"delt";
			break;

		case	HORAIRE:
			wcout << _(L"Coordonnées Horaires");
			unit[COORD_TETA] = hunit;
			symbol2 = L"H   ";
			symbol3 = L"delt";
			break;

		case	HORIZON:
			wcout << _(L"Coordonnées Horizontales");
			symbol2 = L"A   ";
			symbol3 = L"h   ";
			break;
	}


	// Separation
	separation(maxastr);

	wcout << setfill(L' ');

	// Noms
	if(coord==HELIO)
		wcout << endl << setw(9) << L'|' << setw(8) << Planete[TER].nom << L'|' << setw(8) << Lune.nom << L'|';
	else
		wcout << endl << setw(9) << L'|' << setw(8) << Soleil.nom << L'|' << setw(8) << Lune.nom << L'|';


	for(i=0 ; i < maxastr ; i++)
		wcout << setw(8) << Planete[i].nom << L'|';


	if(coord<EQUA)
		AfficheLigneCoord(symbol1, unit[COORD_R], 8, 5, coord, COORD_R, false);
	AfficheLigneCoord(symbol2, unit[COORD_TETA], 8, 4, coord, COORD_TETA, hms);
	AfficheLigneCoord(symbol3, unit[COORD_PHI], 8, 4, coord, COORD_PHI, hms);
}


void	Ephemerides::AfficheLigneCoord(const wstring &symbol, const wstring &unit, unsigned char width, unsigned char precision, unsigned int coord, unsigned int coorditem, bool hms)
{
	unsigned short int i;
	double	co,dint;
	bool	heure=(coord==EQUA || coord==HORAIRE) && coorditem==COORD_TETA;


	wcout << setfill(L' ');
	wcout << endl << setw(4) << symbol << L' ' << setw(3) << unit << L'|';


	// Soleil ou Terre
	if(coord==HELIO)
		co=Planete[TER].Coord[coord].Spher[coorditem];
	else
		co=Soleil.Coord[coord].Spher[coorditem];

	if(heure)
		co/=KH2D;
	if(hms)
		co=Deci2Sexa(co);

	wcout << fixed << setw(width) << setprecision(precision) << co << L'|';

	// Lune
	co=Lune.Coord[coord].Spher[coorditem];

	if(heure)
		co/=KH2D;
	if(hms)
		co=Deci2Sexa(co);

	if(coord==GEO && coorditem==COORD_R)
	{
		co=Lune.parallaxe*60.0;
		wcout << L"P " << setw(2) << (unsigned int)co << L'\'' << (unsigned int)(modf(co,&dint)*60.0) << L"\"|";
	}
	else
		wcout << fixed << setw(width) << setprecision(precision) << co << L'|';


	// Planetes
	for(i=0;i<maxastr;i++)
	{
		co=Planete[i].Coord[coord].Spher[coorditem];

		if(heure)
			co/=KH2D;

		if(hms)
			co=Deci2Sexa(co);

		wcout << fixed << setw(width) << setprecision(precision) << co << L'|';
	}
}


void	Ephemerides::AfficheAPP(void)
{
	unsigned short int	i;


	// Séparation
	separation(maxastr);

	wcout << endl << _(L"Phase, Elongation, Diamètre Apparent, Magnitude");

	// Séparation
	separation(maxastr);

	wcout << setfill(L' ');

	// Noms
	wcout << endl << setw(9) << L"|" << setw(8) << Soleil.nom << L'|' << setw(8) << Lune.nom << L'|';

	for(i=0;i<maxastr;i++)
		wcout << setw(8) << Planete[i].nom << L'|';

	// phase
	wcout << endl << _(L"Phase  °|") << setw(9) << "|";
	wcout << setw(8) << Lune.phase << L'|';

	for(i=0;i<maxastr;i++)
		wcout << setw(8) << Planete[i].phase << L'|';

	// phase en %
	wcout << endl << _(L"Phase  %|") << setw(9) << "|";
	wcout << setw(8) << Lune.phase100 << L'|';

	for(i=0;i<maxastr;i++)
		wcout << setw(8) << Planete[i].phase100 << L'|';

	// élongation
	wcout << endl << _(L"Elong  °|") << setw(9) << "|";
	wcout << setw(8) << Lune.elongation << L'|';

	for(i=0;i<maxastr;i++)
		wcout << setw(8) << Planete[i].elongation << L'|';

	// diamètre apparent
	wcout << endl << _(L"Diam   \"|") << setw(4) << (unsigned int)Soleil.diametre/60 << L'\'' << setw(2) << (unsigned int)Soleil.diametre%60 << L"\"|";
	wcout << setw(4) << (unsigned int)Lune.diametre/60 << L'\'' << setw(2) << (unsigned int)Lune.diametre%60 << L"\"|";

	for(i=0;i<maxastr;i++)
		wcout << setw(8) << fixed << setprecision(2) << Planete[i].diametre << L'|';

	// magnitude, âge pour la Lune
	i=(unsigned short int)Lune.age;
	wcout << endl << _(L"Magnit. |") << setw(8) << fixed << setprecision(2) << Soleil.magnitude << L'|';
	wcout << L"A " << setw(2) << i << _(L'j') << setw(2) << (unsigned short int)((Lune.age-i)*24.0) << _(L"h|");

	for(i=0;i<maxastr;i++)
		wcout << setw(8) << fixed << setprecision(2) << Planete[i].magnitude << L'|';
}



void	Ephemerides::AfficheLMC(void)
{
	unsigned int	i;


	// Séparation
	separation(maxastr);

	wcout << endl << _(L"Heures des Lever, Passage au Méridien, Coucher (heures locales)");

	// Séparation
	separation(maxastr);

	wcout << setfill(L' ');

	// Noms
	wcout << endl << setw(9) << L'|' << setw(8) << Soleil.nom << L'|' << setw(8) << Lune.nom << L'|';

	for(i=0;i<maxastr;i++)
		wcout << setw(8) << Planete[i].nom << L'|';

	// levers
	wcout << endl << _(L"Lev  hms|");
	AfficheHeureLMC(Soleil.Lever.tmlocal);
	AfficheHeureLMC(Lune.Lever.tmlocal);

	for(i=0;i<maxastr;i++)
		AfficheHeureLMC(Planete[i].Lever.tmlocal);

	// passages au méridien
	wcout << endl << _(L"Mér  hms|");
	AfficheHeureLMC(Soleil.Meridien.tmlocal);
	AfficheHeureLMC(Lune.Meridien.tmlocal);

	for(i=0;i<maxastr;i++)
		AfficheHeureLMC(Planete[i].Meridien.tmlocal);

	// couchers
	wcout << endl << _(L"Cou  hms|");
	AfficheHeureLMC(Soleil.Coucher.tmlocal);
	AfficheHeureLMC(Lune.Coucher.tmlocal);

	for(i=0;i<maxastr;i++)
		AfficheHeureLMC(Planete[i].Coucher.tmlocal);
}


void	Ephemerides::AfficheHeureLMC(const struct tm &tmt)
{
	if(tmt.tm_hour == TUSING)
		wcout << setfill(L'-') << setw(9) << L'|';
	else
		wcout << setfill(L' ') << setw(2) << tmt.tm_hour << L':' << setfill(L'0') << setw(2) << tmt.tm_min << L':' << setw(2) << tmt.tm_sec << L'|';
}
