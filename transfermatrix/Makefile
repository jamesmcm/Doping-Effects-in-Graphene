DRIVER   = driver.f
CORESRC  = linalg.f tmatrix.f tmatrixx.f tmatrixy.f util.f trupdate.f
CORESRC += unitarity.f
HELPERS  = 

OBJ = $(DRIVER:.f=.o) $(CORESRC:.f=.o) $(HELPERS:.c=.o)

TARGET = tmsim

INC = 
LIB = -lm -lblas -llapack

LDFLAGS  = $(LIB)
CFLAGS   = $(INC) -Wall
F77FLAGS = $(INC) -Wall

$(TARGET): $(OBJ)
	gcc -o $(TARGET) $(OBJ) $(LDFLAGS)

%.o: %.c
	gcc $(CFLAGS) -o $@ -c $< 

%.o: %.f
	gfortran $(F77FLAGS) -o $@ -c $<
clean::
	-rm $(TARGET) *.o
