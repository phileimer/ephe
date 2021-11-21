#include <string>
using namespace std;


class	Observateur
{
	public:
		// données
		wstring	nom;		// nom du lieu d'observation
		double	longitude;
		double	latitude;
		long 	zone;		// fuseau horaire, en secondes - Paris=+1 en hiver -> +3600
		int		altitude;

		Instant	Inst;		// définition des dates et temps

		double	refraction;	// angle de réfraction

		double	temperature;	// temperature en °C
		double	pression;	// pression en hPa


		// données calculées
		double	eta1,eta2L,eta2C;	// pour le calcul des lever et coucher des astres

		double	InitTemps(void);
		double	InitTemps(unsigned short int, unsigned short int, int, unsigned short int=0, unsigned short int=0, unsigned short int=0, bool_tps=TPSLOC);
		void	InitLieu(const wchar_t *, double, double, int, long);
		bool	ChargeLieu(void);

		void	Affiche(void) const;
};
