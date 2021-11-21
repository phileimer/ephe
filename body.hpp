class Astre
{
	protected:
        Vecteur	ElemOrbit2EcliptHelio(const Orbite &) const;

		// Calcul du temps de trajet de la lumière
		double	TempsLumiere(double) const;

		// conversions héliocentriques <--> géocentriques
		Vecteur	EcliptHelio2EcliptGeo(const Vecteur &,const Vecteur &) const;
		Vecteur	EcliptGeo2EcliptHelio(const Vecteur &,const Vecteur &) const;

		// conversions géocentriques <--> apparentes
		Vecteur	EcliptGeo2EcliptGeoApp(const Vecteur &,double,double) const;
		Vecteur	EcliptGeoApp2EcliptGeo(const Vecteur &,double) const;

		// conversions géocentriques apparentes <--> équatoriales
		Vecteur	EcliptGeoApp2Equa(const Vecteur &) const;
		Vecteur	Equa2EcliptGeoApp(const Vecteur &) const;

		// conversions équatoriales <--> horaires
		Vecteur	Equa2Horaire(const Vecteur &) const;
		Vecteur	Horaire2Equa(const Vecteur &) const;

		// conversions horaires <--> horizontales
		Vecteur	Horaire2Horizon(const Vecteur &) const;
		Vecteur	Horizon2Horaire(const Vecteur &) const;

		// calcul des lever et coucher
		double	Eq2tu(const Vecteur &, double [], unsigned int);
		Vecteur	tu2Eq(double);
		bool	TUcond(double *,double,double);

		double	Parallaxe(void);


	public:
		Ephemerides	*pE;	// pour l'accès à tous les éléments

		std::wstring	nom;		// nom de l'astre
		unsigned short int	astreid;	// numéro d'identification


		Orbite	ElemOrbit;	// éléments de l'orbite

		Vecteur	Coord[6];	// coordonnées

		Instant	Lever;		// instant du lever
		Instant	Coucher;	// instant du coucher
		Instant Meridien;	// instant du passage au méridien

		double	magnitude;	// magnitude relative
		double	diametre;	// diametre apparent
		double	phase;		// phase en degrés
		double	phase100;	// phase en %
		double	elongation;	// élonogation
		double	parallaxe;	// paralaxe


		void	Init(Ephemerides *, unsigned int);

		// coordonnées écliptiques à partir des éléments de l'orbite
		Vecteur	EcliptGeo(void);
		Vecteur	EcliptGeo2Horizon(double);

		double	Phase(void);
		double	Elongation(void);
		double	Diametre(void);
		double	Magnitude(void);
		void	LeverCoucher(void);		// calcul des instants de lever, coucher, passage au méridien

		// formatage des résultats
		void	Affiche(unsigned int, bool=true) const;
};


class	AstreLune : public Astre
{
	public:
		OrbiteLune	ElemOrbitLune;
		double	age;

		void	Init(Ephemerides *);
		double	Age(void);
		Vecteur	EcliptGeo(void);
		Vecteur	ElemOrbit2EcliptGeo(const OrbiteLune &, double, double * =NULL);
};


class	Nutation
{
	public :
		Ephemerides	*pE;

		double	lamda;		// nutation en longitude géocentrique
		double	epsilon;	// nutation en obliquité

		void	Calcule(Ephemerides *);
		void	Affiche(void);
};
