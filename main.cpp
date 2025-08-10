/* Jean Philippe EIMER - 1999 *
 *  phil.eimer@9online.fr     */

#include <iostream>
using namespace std;


#include "include.h"

#ifndef NO_CONFIG_H
	#include "config.h"
#endif	// NO_CONFIG_H

#ifdef HAVE_LIBINTL
	#include <libintl.h>
#endif // HAVE_LIBINTL


int	main(int argc,char* argv[])
{
	(void)argv;
	Ephemerides	E;


#ifdef HAVE_LIBINTL
	setlocale (LC_ALL, "");
	bindtextdomain (PACKAGE, LOCALEDIR);
	textdomain (PACKAGE);
#endif // HAVE_LIBINTL

	wcout.imbue( locale("") );

//	locale::global(locale("fr_FR.UTF-8"));
//	wcout.imbue(locale("fr_FR.UTF-8"));

	if(argc > 1)
	{
		wcout << L" " << PACKAGE << L" " << VERSION << L" - " << AUTHOR;
	}
	else
	{
//		E.Observ.InitLieu("Paris", -2.3375, 48.8363888889, 0, 3600L);
//		E.Observ.InitLieu("Dijon", -5.0333333, 47.31666667, 0, 3600L);
		if( !E.Observ.ChargeLieu() )
			E.Observ.InitLieu(L"Lyon", -4.8333333, 45.76666667, 248, 3600L);

		E.Observ.InitTemps();
//		E.Observ.InitTemps(4,12,2001,11,23,36,TPSLOC);

		E.InitAstres();

		E.Affiche(HMS);
	}

	return(0);
}
