SICONOSPATH=/local_home/pissard/Workspace/siconos-user/
NUMERICSPATH=/local_home/pissard/Workspace/siconos/trunk/Numerics/

INC_PATH=$(SICONOSPATH)/include
CFLAGS=-g -fbounds-check -w -Wno-deprecated -fPIC -DLA_COMPLEX_SUPPORT -DWITHOUT_NANA -I$(INC_PATH) -I$(NUMERICSPATH)/include -I/usr/local//include/lapackpp/ `xml2-config --cflags`
DYNFLAGS=-shared

siconos.so:siconos.cpp
	g++  $(CFLAGS) -c siconos.cpp
	g++  $(DYNFLAGS) -o siconos.so siconos.o

clean:
	rm siconos.so siconos.o
