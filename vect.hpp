// conversion D�cimal --> Sexag�simal
double	Deci2Sexa(double,int * =NULL,int * =NULL,int * =NULL);


class Vecteur
{
	public:
		double	Spher[3];	// coordonn�es sph�riques
		double	Rect[3];	// coordonn�es rectangulaires


		// constructeurs et destructeur
		Vecteur(double,double,double,bool);
		Vecteur(const Vecteur &v) { *this=v; };
		Vecteur();
		~Vecteur() {};


		// conversions sph�riques <--> rectangulaires
		Vecteur	Rectangulaire(void);
		Vecteur	Spherique(void);

		// rotations autour des axes Ox et Oy
		Vecteur	RotationX(double);
		Vecteur	RotationY(double);

		// surcharge des op�rateurs +, - et =
		Vecteur	operator+(const Vecteur &) const;
		Vecteur	operator-(const Vecteur &) const;
		Vecteur	&operator=(const Vecteur &);

		// formatage des r�sultats
		void	Affiche(unsigned short int,bool,bool) const;
};


class	VectData
{
	public:
		double	d[5];	// 4 �l�ments en standard, 1 suppl�mentaire

		VectData operator =(const double [4]);

		double	Dot(const VectData &) const;
		double	DotT(double) const;
		double	DotTprecis(double) const;
		double	DotTprecis2(double,double) const;

		double	Dot5(const double [5]) const;	// uniquement pour Nutation

		void	Affiche(void) const;
};
