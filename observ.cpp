/* Jean Philippe EIMER - 1999 *
 *  phil.eimer@9online.fr     */

#include <math.h>
#include <stdlib.h>
#include <wchar.h>
#include <string.h>

#include <fstream>
#include <iostream>
#include <iomanip>


#include "include.h"



#define LINE_LENGTH 80
#define	TER_RT_M	(TER_RT*1000.0)


float	LitSexa(wchar_t *line)
{
	float l, result;


	l = wcstof(line, NULL);
	result = fabs(l);
	result += wcstof(wcschr(line, L' '), NULL)/60.0;
	result += wcstof(wcsrchr(line, L' '), NULL)/3600.0;

	result = copysign(result, l);

	return(result);
}


bool	Observateur::ChargeLieu(void)
{
	wifstream	*file;
	char *home, filename[LINE_LENGTH];
	bool	ret = false;


    if (!(home = getenv("HOME")) || strlen(home) >= LINE_LENGTH)
		strcpy(filename, ".");
	else
		strcpy(filename, home);

	strcat(filename, "/.epherc");


	file = new wifstream(filename, ios_base::in);

	if( file->fail() )
		wcout << endl << _(L"Fichier de configuration ~/.epherc non trouvé.") << endl << _(L"Utilisation des coordonnées par défaut.") << endl << endl;
	else
	{
		wchar_t line[LINE_LENGTH],word[LINE_LENGTH];

		do
		{
			file->getline(line, LINE_LENGTH);
		}
		while(line[0] == L'#');

		unsigned char	index2 = wcscspn(line, L"\t");
		nom = wstring(line, (index2 > 15 ? 16 : index2));
		//	strncpy(nom, line, (index2>15?16:index2));

		unsigned char 	index1 = index2 + 1;
		index2 = wcscspn(line + index1, L"\t");
		wcsncpy(word, line + index1, index2);
		word[index2]=L'\0';
		longitude = LitSexa(word);

		index1 += index2 + 1;
		index2 = wcscspn(line + index1, L"\t");
		wcsncpy(word, line + index1, index2);
		word[index2]=L'\0';
		latitude = LitSexa(word);

		index1 += index2 + 1;
		wcscpy(word, line + index1);
		altitude = wcstol(word, NULL, 10);

		index1 += index2 + 1;
		wcscpy(word, line + index1);
		zone = wcstol(word, NULL, 10);
		zone *= 3600.0;

		refraction = 0.61;
		eta1 = RAD2DEG(acos(TER_RT_M/(TER_RT_M+altitude)));
		eta2L = RAD2DEG(0);	// atan(h/d);
		eta2C = RAD2DEG(0);	// atan(h/d);

		file->close();
		delete(file);

		ret = true;
	}

	return( ret );
}


void	Observateur::InitLieu(const wchar_t *n, double lng, double lat, int alt, long tz)
{
	nom = wstring(n);
	longitude = lng;
	latitude = lat;
	altitude = alt;
	zone = tz;

	refraction = 0.61;
	eta1 = RAD2DEG(acos(TER_RT_M/(TER_RT_M+altitude)));
	eta2L = RAD2DEG(0);	// atan(h/d);
	eta2C = RAD2DEG(0);	// atan(h/d);
}


double	Observateur::InitTemps(void)
{
	time_t	timet;
	struct tm	*tmtl;

	time(&timet);
	tmtl=localtime(&timet);
	Inst.tmlocal = *tmtl;

	tmtl=gmtime(&timet);
	Inst.tmtu = *tmtl;

#ifndef _WIN32
	zone = Inst.tmlocal.tm_gmtoff;
#endif	// _WIN32

	Inst.Init(longitude, zone, TPSLOC);


	return(Inst.tu);
}


//double	Observateur::InitTemps(unsigned short int jour,unsigned short int mois,int annee,unsigned short int heure=0,unsigned short int min=0, unsigned short int sec=0, bool btu=FALSE)
double	Observateur::InitTemps(unsigned short int jour, unsigned short int mois, int annee, unsigned short int heure, unsigned short int min, unsigned short int sec, bool_tps btu)
{
/*	time_t	timet;
	struct tm	*tmtl;
*/

	if(btu)		// heure entrée en TU
	{
		Inst.tmtu.tm_mday = jour;
		Inst.tmtu.tm_mon = mois-1;
		Inst.tmtu.tm_year = annee-1900;

		Inst.tmtu.tm_hour = heure;
		Inst.tmtu.tm_min = min;
		Inst.tmtu.tm_sec = sec;

/*		timet = timegm(&Inst.tmtu);
		tmtl = localtime(&timet);
		Inst.tmlocal = *tmtl;

		Inst.Init(longitude);
*/	}
	else		// heure entrée comme heure locale
	{
		Inst.tmlocal.tm_mday = jour;
		Inst.tmlocal.tm_mon = mois-1;
		Inst.tmlocal.tm_year = annee-1900;

		Inst.tmlocal.tm_hour = heure;
		Inst.tmlocal.tm_min = min;
		Inst.tmlocal.tm_sec = sec;

/*		Inst.tmlocal.tm_isdst = 0;

		timet = mktime(&Inst.tmlocal);
		tmtl = localtime(&timet);
		Inst.tmlocal = *tmtl;		// pour fixer le DST

		Inst.tmlocal.tm_mday = jour;
		Inst.tmlocal.tm_mon = mois-1;
		Inst.tmlocal.tm_year = annee-1900;

		Inst.tmlocal.tm_hour = heure;
		Inst.tmlocal.tm_min = min;
		Inst.tmlocal.tm_sec = sec;
		timet = mktime(&Inst.tmlocal);
		tmtl = gmtime(&timet);
		Inst.tmtu = *tmtl;

		Inst.Init(longitude);
*/	}

//	zone = Inst.tmlocal.tm_gmtoff;
	Inst.Init(longitude, zone, btu);

	return(Inst.tu);
}


void	Observateur::Affiche(void) const
{
	Inst.Affiche();

	wcout << endl << _(L"Lieu de l'Observation") << L" : " << nom;
	wcout << endl << _(L"Longitude") << L" : " << (longitude>0?_(L"O"):_(L"E")) << L' ' << setw(8) << fixed << setprecision(4) << fabs(Deci2Sexa(longitude)) << L" °\'\"";
	wcout << endl << _(L"Latitude ") << L" : " << (latitude>0?_(L"N"):_(L"S")) << L' ' << setw(8) << fixed << setprecision(4) << fabs(Deci2Sexa(latitude)) << L" °\'\"";
	wcout << endl << _(L"Altitude ") << L" : " /*<< setw(4)*/ << altitude << L" m";
}
