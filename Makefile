# Copyright 2014 Ives Rey-Otero <ives.rey-otero@cmla.ens-cachan.fr>

# compilers configuration
CC = gcc
OFLAGS = -O3
LIBS =  -lpng -lfftw3 -L/usr/local/lib -lm

CFLAGS = -Wall -Werror -Wno-write-strings -pedantic -std=c99 -D_POSIX_C_SOURCE=200809L

SRC_ALGO = gaussconv_sampled_kernel \
		   gaussconv_lindeberg \
		   gaussconv_dft \
		   gaussconv_dct \
		   demo_gaussconv_sampled_kernel \
		   demo_gaussconv_lindeberg \
		   demo_gaussconv_dft \
		   demo_gaussconv_dct

SRC2 = lib_gaussconv.c lib_io_gaussconv.c

SRCDIR = src
OBJDIR = src
BINDIR = bin

BIN = $(addprefix $(BINDIR)/,$(SRC_ALGO))
OBJ = $(addprefix $(OBJDIR)/,$(SRC2:.c=.o))

default: $(BIN)

# binary files
$(BIN) : $(BINDIR)/% : $(SRCDIR)/%.c $(OBJ) $(OBJDIR)/io_png.o
	-mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

$(OBJDIR)/%.o : $(SRCDIR)/%.c  $(OBJDIR)/io_png.o
	-mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(OFLAGS) -c -o $@ $<

$(OBJDIR)/io_png.o : $(SRCDIR)/io_png.c
	    -mkdir -p $(OBJDIR)
		    $(CC) $(CFLAGS) $(OFLAGS) -c -o $@ $<

cleanobj:
	-rm -f $(OBJ) $(OBJDIR)/io_png.o


clean: cleanobj
	-rm -f $(BIN)
