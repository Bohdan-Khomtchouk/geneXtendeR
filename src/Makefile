C = clang
LOADER = clang

OBJECTS: $(OBJECTS)
	R CMD SHLIB extract_peaks.c
	R CMD SHLIB extract_number.c
	R CMD SHLIB annotate.c 
	 $(OBJECTS)
 



