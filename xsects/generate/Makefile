
all: xsects

xsects.dat: xsects
	./xsects

xsects: WXSECTS.F90 PHFIT2.F PHOTO.F
	gfortran -Wall -Wextra -Wno-tabs WXSECTS.F90 PHFIT2.F PHOTO.F -o xsects

clean:
	rm xsects


