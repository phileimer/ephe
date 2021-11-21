prog = $(CROSS_COMPILE)g++
compflags = -c -O2 -Wall $(addflags) #-g
debugflag =
# -ggdb
#linkflags = -lstdc++

exe = ephe
debugexe = ephe

objets = main.o vect.o body.o orbit.o observ.o instant.o ephe.o
inc = include.h data.h defs.h vect.hpp ephe.hpp orbit.hpp instant.hpp observ.hpp body.hpp


default: $(exe)
	$(CROSS_COMPILE)strip $(exe)

$(exe) : $(objets)
	$(prog) $^ $(linkflags) -o $@

$(objets) : %.o : %.cpp $(inc) 
	$(prog) $(compflags) $(debugflag) $< -o $@

distclean: clean
	rm -f $(exe) *~

clean:
	rm -f *.o
