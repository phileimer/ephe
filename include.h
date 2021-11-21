#if TIME_WITH_SYS_TIME
# include <sys/time.h>
# include <time.h>
#else
# if HAVE_SYS_TIME_H
#  include <sys/time.h>
# else
#  include <time.h>
# endif
#endif

#include <string>


/* Pour i18n */
//#include <libintl.h>
//#define _(String) gettext(String)
#define _(String) String

#include "defs.h"
#include "data.h"

class	Ephemerides;
class	Astre;
class	AstreLune;

#include "vect.hpp"
#include "instant.hpp"
#include "orbit.hpp"
#include "body.hpp"
#include "observ.hpp"
#include "ephe.hpp"

