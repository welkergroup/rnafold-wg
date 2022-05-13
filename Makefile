VIENNA_RNA_VERSION = 2.4.18

.ONESHELL: # Applies to every targets in the file!

all: build

ViennaRNA-$(VIENNA_RNA_VERSION).tar.gz:
	wget -q https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-$(VIENNA_RNA_VERSION).tar.gz

ViennaRNA-$(VIENNA_RNA_VERSION): ViennaRNA-$(VIENNA_RNA_VERSION).tar.gz
	tar zxvf ViennaRNA-$(VIENNA_RNA_VERSION).tar.gz

vr-build: ViennaRNA-$(VIENNA_RNA_VERSION)
	cd ViennaRNA-$(VIENNA_RNA_VERSION)
	./configure --without-perl --without-python --without-python3 --without-forester --without-rnalocmin
	make

build: vr-build
	gcc -Wall -fopenmp -I ViennaRNA-$(VIENNA_RNA_VERSION)/src -o rnafold-wg rnafold-wg.c -L ViennaRNA-$(VIENNA_RNA_VERSION)/src/ViennaRNA/.libs -lRNA -lm

clean:
	rm -f ViennaRNA-$(VIENNA_RNA_VERSION) ViennaRNA-$(VIENNA_RNA_VERSION).tar.gz
	rm -f rnafold-wg
