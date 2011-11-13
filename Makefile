CC = gcc
LD = gcc
GSLDIR = /usr/local
GSLINC = $(GSLDIR)/include
GSLLIB = $(GSLDIR)/lib
CFLAGS = -I$(GSLINC) -c
LDFLAGS = -L$(GSLLIB) -lgsl -lm
TARGET = mprojectgsl


$(TARGET): $(TARGET).o
	$(LD)  -o $@ $< $(LDFLAGS)
 
$(TARGET).o: $(TARGET).c
	$(CC) -o $@ $(CFLAGS) $<
