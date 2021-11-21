#define	NO_CONFIGURE
#define	PACKAGE	"ephe"
#define	VERSION	"2.0"
#define	AUTHOR	"Jean Philippe EIMER (phil.eimer@9online.fr)\n"
#define	LOCALEDIR	"/usr/local/share/locale"

//#define FALSE	0
//#define TRUE	!0

#ifndef NULL
#define	NULL	0L
#endif

#ifndef TUSING
#define	TUSING	-100.0
#endif

#define KH2D	15.0		// angles en heures vers angles en degres

#define K360(x) ((x)<0.0?fmod(x,360.0)+360.0:fmod(x,360.0))	// entre 0 et 360
#define KNEG(x) ((x) > 180.0 ? (x)-360.0: (x))			// entre -180 et 180
#define K24(x)	((x)<0.0?fmod(x,24.0)+24.0:fmod(x,24.0))	// entre 0 et 24 h

#define RAD2DEG(x)	((x)*180.0/M_PI)	// conversion radians --> degres
#define DEG2RAD(x)	((x)*M_PI/180.0)	// conversion degres --> radians


#define	PRECISION_KEPLER	0.000005
#define	PRECISION_LC		1.0/3600.0	// précision 1 min pour le lever et le coucher
#define	PRECISION_M		1.0/3600.0	// précision 1 sec pour le passage au méridien


typedef	enum {TPSLOC, TPSTU} bool_tps;
enum {HELIO,GEO,GEOAPP,EQUA,HORAIRE,HORIZON};
enum {COORD_R,COORD_TETA,COORD_PHI};
enum {COORD_X,COORD_Y,COORD_Z};
enum {MER,VEN,MAR,JUP,SAT,URA,NEP,PLU,TER,SOL,LUN};
enum {DECI,HMS};
enum {LEVER,COUCHER,MERIDIEN};

#define	NOZONE	50000L

