#include <string>
using namespace std;


class	Ephemerides
{
	void	AfficheCoord(unsigned int, bool);
	void	AfficheLigneCoord(const wstring &, const wstring &, unsigned char, unsigned char, unsigned int, unsigned int, bool);
	void	AfficheAPP(void);
	void	AfficheLMC(void);
	void	AfficheHeureLMC(const struct tm &);


	public:
		Observateur	Observ;

		Astre	Soleil;
		Astre	Planete[9];
		AstreLune	Lune;

		Nutation	Nutat;

		double	epsilon;	// obliquité : inclinaison de l'équateur terrestre par rapport à l'écliptique

		unsigned short int	maxastr;	// arrêt des calculs à cet astreid

		void	InitAstres(void);
		void	Obliquite(void);
		void	Affiche(bool=true);
};
