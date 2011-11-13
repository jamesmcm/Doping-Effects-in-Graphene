CC = gcc
LD = gcc
GSLDIR = /opt/local
GSLINC = $(GSLDIR)/include
GSLLIB = $(GSLDIR)/lib
CFLAGS = -I$(GSLINC) -c
LDFLAGS = -L$(GSLLIB) -lgsl
TARGET = mprojectgsl


$(TARGET): $(TARGET).o
	$(LD)  -o $@ $< $(LDFLAGS)
 
$(TARGET).o: $(TARGET).c
	$(CC) -o $@ $(CFLAGS) $<
