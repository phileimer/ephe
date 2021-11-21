/* Jean Philippe EIMER - 1999 *
 *  phil.eimer@9online.fr     */

#include <math.h>
#include <wchar.h>

#include <iostream>
#include <iomanip>
#include <ctime>
using namespace std;

#include "include.h"



double	Instant::Init(double longitude, long zone, bool_tps btu)
{
	double	dj;


	if (btu == TPSTU)
	{
		dj = Julien(tmtu, NULL);
		dj += (double)zone / 86400.0;
		InvJulien(tmlocal, dj);
	}
	else
	{
		dj = Julien(tmlocal, NULL);
		dj -= (double)zone / 86400.0;
		InvJulien(tmtu, dj);
	}

	datejulienne=Julien(tmtu, &tu);
	dj2t();
	Sideral(longitude);

	return(datejulienne);
}


//double	Instant::tu2t(const Instant &Idate,double tuext,short int zone=127)
double	Instant::tu2t(const Instant &Idate, double tuext, long zone)
{
	if(tuext == TUSING)
	{
		tmtu.tm_mday=0;
		tmtu.tm_mon=0;
		tmtu.tm_year=0;
		tmtu.tm_hour=0;
		tmtu.tm_min=0;
		tmtu.tm_sec=0;
		tmlocal.tm_mday=0;
		tmlocal.tm_mon=0;
		tmlocal.tm_year=0;
		tmlocal.tm_hour=TUSING;
		tmlocal.tm_min=0;
		tmlocal.tm_sec=0;
		datejulienne=0;
		t=0;
	}
	else
	{
		tmtu.tm_mday = Idate.tmtu.tm_mday;
		tmtu.tm_mon = Idate.tmtu.tm_mon;
		tmtu.tm_year = Idate.tmtu.tm_year;
		tmtu.tm_hour=(unsigned short int)tuext;
		tuext-=tmtu.tm_hour;
		tuext*=60.0;
		tmtu.tm_min=(unsigned short int)tuext;
		tuext-=tmtu.tm_min;
		tmtu.tm_sec=(unsigned short int)(tuext*60.0);

		datejulienne=Julien(tmtu,&tu);
		dj2t();

		if(zone != NOZONE)
			DiffHeure(true, zone);
	}

	return(t);
}


//double	Instant::Julien(const struct tm &tmt, double *phdec=NULL)
double	Instant::Julien(const struct tm &tmt, double *phdec)
{
	double	dj,hdec;
	unsigned int	a=tmt.tm_year+1900;
	unsigned int	a1=a;
	unsigned short int	m=tmt.tm_mon+1;
	unsigned short int	aa;
	unsigned short int	m1=m+1;


	if(m<3)
	{
		a1--;
		m1+=12;
	}

	aa=a1/100;

	dj=1720994.5+floor(365.25*a1)+floor(30.6*m1)+tmt.tm_mday;

	if((a*100+m)*100+tmt.tm_mday > 15821004)
		dj+=(double)(2-aa+aa/4);

	hdec=hDeci(tmt);
	dj+=hdec/24.0;

	if(phdec!=NULL)
		*phdec=hdec;

	return(dj);
}


double	Instant::InvJulien(struct tm &tmt, double dj)
{
	double	hdec,dint;
	long int a;
	int	b,c,d,e,g;


	a=(long int)floor(dj+0.5);


	if(a>2299160)
	{
		double f=(int)(((double)a-1867216.25)/36524.25);
		a+=1+f-f/4;
	}

	b=(int)(a+1524);
	c=(int)floor(((double)b-122.1)/365.25);
	d=(int)floor((double)c*365.25);
	e=(int)((double)(b-d)/30.6001);
	g=(int)((double)e*30.6);

	tmt.tm_mday = b-d-g;
	tmt.tm_mon = (e<14?e-1:e-13);
	tmt.tm_year = (tmt.tm_mon>2?c-4716:c-4715);
	tmt.tm_mon -= 1;
	tmt.tm_year -= 1900;

	hdec = modf(dj-0.5,&dint)*24.0;

	Deci2Sexa(hdec,&g,&tmt.tm_min,&tmt.tm_sec);
	tmt.tm_hour = (unsigned short int)g;

	return(hdec);
}


void	Instant::DiffHeure(bool tudefined, double zone)
{
	if(!tudefined)
	{
		InvJulien(tmtu,Julien(tmlocal)-zone/86400.0);
		tmtu.tm_sec = tmlocal.tm_sec;	// conservation des secondes
	}
	else
	{
		InvJulien(tmlocal,datejulienne+zone/86400.0);
		tmlocal.tm_sec = tmtu.tm_sec;	// conservation des secondes
	}
}


double	Instant::Sideral(double longitude)
{
	double	vtsg0[4]=TSG0,t1;
	VectData	V;

	V=vtsg0;

	t1=(floor(datejulienne-0.5)+0.5-DJ0)/36525.0-1.0;	// t-1 � 0hTU

	tsg0=V.DotT(t1);
	tsg0+=(t1-((double)tmtu.tm_year+1900-2000.0)/100.0)*2400.0;
	tsg0=K24(tsg0);

	tsg=tsg0+tu*SOL2SID;
	tsg=K24(tsg);

	ts=tsg-longitude/KH2D;
	ts=K24(ts);

	return(ts);
}



double	Instant::hDeci(const struct tm &tmt) const
{
	double	hdec;

	hdec=tmt.tm_hour+(tmt.tm_min+tmt.tm_sec/60.0)/60.0;
	hdec=K24(hdec);

	return(hdec);
}


double	Instant::dj2t(void)
{
	t=(datejulienne-DJ0)/36525.0;

	return(t);
}


void	Instant::Affiche(void) const
{
	wchar_t str[60];

	wcsftime(str, 60, L"%A %d %B %Y, %T", &tmtu);
	wcout << endl << _(L"Date et Heure TU     ") << L" : " << str;
	wcsftime(str, 60, L"%A %d %B %Y, %T, %z (%Z)", &tmlocal);
	wcout << endl << _(L"Date et Heure Locales") << L" : " << str;

	wcout << endl << _(L"Date Julienne        ") << L" : " << fixed << setprecision(4) << datejulienne;

	wcout << endl << _(L"Temps Sidéral ") << L" (h,ms) : " << _(L"TSG") << L"0=" << Deci2Sexa(tsg0) << L'\t' << _(L"TSG") << L"=" << Deci2Sexa(tsg) << L'\t' << _(L"TS") << L"=" << Deci2Sexa(ts) << endl;
}
