In order to use spec2col.py you need:
-cmf_precise.txt
-a "spectra" folder containing your chosen PHOENIX spectra (see line 12 in this file) in a ".fits" format or 
 different spectra in the ".csv" or ".txt" format as a table with (column 1: wavelength, column 2: flux)
	For your own spectra in a ".csv" or ".txt" format:
	- the wavelength resolution should be better than 1AA
	- if the wavelength range in your spectrum is only partly in the range of 3900AA to 8300AA, it needs to be manually extended so that the whole range is covered. Steps of 1AA are accurate enough.
	- if the data is only available as ".txt" files: the standard delimiter is set as ",". It can be changed in line 184 in "spec2col.py".
	  

- PHOENIX spectra are recommended, which can be found under:
	"http://phoenix.astro.physik.uni-goettingen.de/data/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/"
  You will also need a wavelength file, which can be found for the PHOENIX spectra under:
	"http://phoenix.astro.physik.uni-goettingen.de/data/HiResFITS" with the name "WAVE_PHOENIX-ACES-AGSS-COND-2011.fits".
  This file needs to be placed in the same directory as "spec2col.py".


The resulting color(s) will be put into the file "color_result.txt".
