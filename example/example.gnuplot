set datafile separator comma
set terminal pngcairo size 1200,800
set output 'example.png'
set xlabel 'Single-photon wavelength / nm'
set ylabel 'TPA cross section / GM'
set title 'TPA spectrum: linear, lorentzian, width=0.1 eV hwhm, N=4'
set xrange [*:*] reverse
plot 'example.curve.csv' using 1:6 with lines title 'TPA'
