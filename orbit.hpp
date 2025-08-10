class	Orbite
{
	Orbite	operator =(const double [6][4]);
	void	PerturbationsSoleil(double);	// calcul des perturbations pour le soleil
	void	PerturbationsTellur(double);	// perturbations pour Mercure, V�nus et Mars
	void	PerturbationsGeantes(double);	// perturbations pour les plan�tes g�antes
	double	AnomalieMoyenne(void);	// calcul de M et omega
	double	Kepler(void);		// �quation de Kepler (calcul de E)
	double	AnomalieVraie(void);	// calcul de r et v


	public:
		Astre	*pA;	// pour l'acc�s � tous les �l�ments de l'astre correspondant

		// Vecteurs d'initialisation des �l�ments orbitaux
		VectData	Ve;
		VectData	Vi;
		VectData	VL;
		VectData	VO;
		VectData	Vob;

		// El�ments calcul�s
		double	a;	// demi grand axe
		double	e;	// excentricit�
		double	i;	// inclinaison sur l'�cliptique
		double	L;	// longitude moyenne
		double	Omega;	// longitude du noeud ascendant
		double	omegab;	// longitude du p�rih�lie
		double	deltaL;	// variation �l�mentaire de L, pour la prise en compte du temps de trajet de la lumi�re

		// Perturbations calcul�es
		double	aper;	// demi grand axe
		double	eper;	// excentricit�
		double	Lper;	// longitude moyenne
		double	vkper;	// longitude du p�rih�lie (associ�e avec e)
		double	rper;	// rayon vecteur
		double	lper;	// longitude
		double	bper;	// latitude

		// El�ments calcul�s d'apr�s les pr�cedents
		double	M;	// anomalie moyenne = L-omegab
		double	MRsp;	// anomalie moyenne sans perturbations, en radians
		double	omega;	// argument de latitude du p�rih�lie = omegab-Omega

		double	E,ER;	// anomalie excentrique (obtenue par Kepler) en degr�s et radians

		double	r;	// rayon vecteur
		double	v;	// anomalie vraie


		Orbite() { pA=NULL; a=e=i=L=Omega=omegab=deltaL=aper=eper=Lper=vkper=rper=lper=bper=M=MRsp=omega=E=ER=r=v=0; };
		Orbite(const Orbite&o) { *this = o; };
		void	Init(Astre *);
		Orbite	&operator =(const Orbite &);
		void	Elements(double);	// initialisation des �l�ments
		void	Perturbations(double);	// initialisation des perturbations
		void	Anomalies(void);

		void	Affiche(void) const;		// affichage
};


class	OrbiteLune
{
	public:
		AstreLune	*pA;

		VectData	VO,VL,VM,VF,VD;				// �l�ments de l'orbite
		VectData	VLper,VMper,VFper,VDper;	// pour les perturbations
		VectData	Va,Vb,Vc,Ve;

		double		L,M,F,D;
		double		Omega,OmegaR,CCR,EE;


		void	Init(AstreLune *);
		void	Elements(const Instant &);
};
