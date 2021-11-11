# https://www.gnu.org/software/make/manual/html_node/Quick-Reference.html
COMP=g++
CFLAGS=-Wall -O2 -std=c++11

MAJOR=0
MINOR=0
LIBNAME=test

OBJ=$(addsuffix .o,$(basename $(wildcard *.cpp)))
HEADERFILES=$(wildcard *.h)

prog: $(OBJ)
	$(COMP) $(CFLAGS) -o prog $^

%.o: %.cpp $(HEADERFILES)
	$(COMP) $(CFLAGS) -c $<

mklib: lib$(LIBNAME).so.$(MAJOR).$(MINOR)


lib$(LIBNAME).so.$(MAJOR).$(MINOR):
	g++ $(CFLAGS) -fPIC -c *.cpp
	g++ $(CFLAGS) -Wl,-soname,lib$(LIBNAME).so.$(MAJOR) -shared -o $@ $^
	rm *.o

install: lib$(LIBNAME).so.$(MAJOR).$(MINOR)
	sudo cp -v $^ /usr/local/lib/
	sudo ln -s /usr/local/lib/lib$(LIBNAME).so.$(MAJOR).$(MINOR) /usr/local/lib/lib$(LIBNAME).so.$(MAJOR)
	sudo ln -s /usr/local/lib/lib$(LIBNAME).so.$(MAJOR) /usr/local/lib/lib$(LIBNAME).so
	sudo cp -v *.h /usr/local/include/

remove:
	sudo rm -v /usr/local/lib/lib$(LIBNAME).so.$(MAJOR).$(MINOR)
	sudo rm -v /usr/local/lib/lib$(LIBNAME).so.$(MAJOR)
	sudo rm -v /usr/local/lib/lib$(LIBNAME).so
	for i in $$(ls *.h); do sudo rm -v /usr/local/include/$$i; done

clean:
	rm *.o *.so* prog 
