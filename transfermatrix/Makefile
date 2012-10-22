DRIVER   = driver.f
CORESRC  = linalg.f tmatrix.f tmatrixx.f tmatrixy.f util.f trupdate.f vacancy.f potential.f
CORESRC  += unitarity.f rng.f
PYWRAPPER = _tmatrix.c
HELPERS  = gettime.c

OBJ = $(DRIVER:.f=.o) $(CORESRC:.f=.o) $(HELPERS:.c=.o)

TARGET = tmsim 
PY_TARGET = $(PYWRAPPER:.c=.so)
PY_OBJ = $(CORESRC:.f=.o) $(PYWRAPPER:.c=.o)

C_TARGET = c-example gettime
C_SRC = c-example.c gettime.c
C_OBJ = $(CORESRC:.f=.o) $(C_SRC:.c=.o) 

INC = 
LIB = -lm -lblas -llapack -lgfortran

LDFLAGS  = $(LIB)
CFLAGS   = $(INC) -Wall -fPIC
F77FLAGS = $(INC) -Wall -fPIC

PY_INC = -I/usr/include/python2.7
PY_LDFLAGS = $(LDFLAGS) -lpython2.7

all: $(TARGET) $(PY_TARGET) $(C_TARGET)

$(C_TARGET): $(C_OBJ)
	gcc -o $@ $(C_OBJ) $(LDFLAGS)

$(PY_TARGET): $(PY_OBJ)
	gcc -shared -o $@ $(PY_OBJ) $(PY_LDFLAGS)

_tmatrix.o: _tmatrix.c
	gcc $(CFLAGS) $(PY_INC) -o $@ -c $< 

$(TARGET): $(OBJ)
	gcc -o $@ $(OBJ) $(LDFLAGS)

%.o: %.c
	gcc $(CFLAGS) -o $@ -c $< 

%.o: %.f
	gfortran $(F77FLAGS) -c $< -o $@	

clean::
	-rm $(TARGET) *.o *.so
