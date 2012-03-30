DRIVER   = driver.f
CORESRC  = linalg.f tmatrix.f util.f trupdate.f
HELPERS  = 

OBJ = $(DRIVER:.f=.o) $(CORESRC:.f=.o) $(HELPERS:.c=.o)

TARGET = tmsim

INC = 
LIB = -lm -lblas -llapack

LDFLAGS  = $(LIB)
CFLAGS   = $(INC)
F77FLAGS = $(INC)

$(TARGET): $(OBJ)
	gcc -o $(TARGET) $(OBJ) $(LDFLAGS)

%.o: %.c
	gcc $(CFLAGS) -o $@ -c $< 

%.o: %.f
	gfortran $(F77FLAGS) -o $@ -c $<
clean::
	-rm $(TARGET) *.o
