FF=gfortran
src=map.f90
srci=mapi.f90
src4=map4.f90
target=map
targeti=mapi
target4=map4

default:
	$(FF) $(src) -o $(target) 
imag:
	$(FF) $(srci) -o $(targeti)
4:
	$(FF) $(src4) -o $(target4) 
clean:
	rm *.mod
