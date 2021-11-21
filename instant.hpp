
class	Instant
{
	double	hDeci(const struct tm &) const;
	double	dj2t(void);
	double	Sideral(double);

	public:
		struct tm	tmlocal;	// date et heure locale
		struct tm	tmtu;		// date et heure TU

		double	datejulienne;	// date julienne (avec TU)
		double	t;	 	// nombre de jours julien depuis 31/12/1899 à 12h

		double	tu;		// temps universel
		double	tsg0;		// temps sidéral à Greenwich à 0hTU
		double	tsg;   		// temps sidéral à Greenwich
		double	ts;		// temps sidéral du lieu


		double	Init(double, long, bool_tps);
		double	tu2t(const Instant &, double, long=NOZONE);
		double	Julien(const struct tm &, double * =NULL);
		double	InvJulien(struct tm &, double);
		void	DiffHeure(bool, double);

		void	Affiche(void) const;
};
