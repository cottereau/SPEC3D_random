FC=ifort
FFLAGS=-O2  -g
LDFLAGS=-O2  -g
AR=ar

LIB=librandom.a

SRC= $(wildcard *.f)
OBJ= $(SRC:.f=.o)

all: $(LIB)

$(LIB): $(OBJ)
	$(AR) crv $@ $(OBJ) 
	ranlib $@     


%.o: %.f
	$(FC) $(FFLAGS) -c -lmpi  $< 

.PHONY: clean 

clean:
	rm -rf *.o $(LIB)



