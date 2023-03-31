set terminal post port
set output "tmp.ps"
set nogrid
set tics in
set ticslevel 0.5
set zero 1e-15
set logscale xy 10
set size 1,0.7
set title "PARTICLE SPECTRA: positrons
set xlabel "Kinetic energy, GeV" 0.,0.
set ylabel "Flux, m^-2 s^-1 sr^-1 GeV^-1" 0.,0.
set zlabel "" 0.000000,0.000000
set label "Phi = 500 MV" at 1,0.0004
set nokey
plot [0.1:1e+02][8e-05:*] "42.3ds_DarkSUSY" w l 1,"42.3ds_DarkSUSY" u 1:3 w l 1,"42.3ds_DarkSUSY" u 1:4 w l 1,"42.3ds_DarkSUSY" u 1:5 w l 1 , "positrons0_boezio00_gev.dat" w l 1, "positrons0_boezio00_gev.dat" u 3:4 w p 5, "positrons0_barwick98_gev.dat" w l 1, "positrons0_barwick98_gev.dat" u 3:4 w p 6, "positrons0_grimani02_gev.dat" w l 1, "positrons0_grimani02_gev.dat" u 3:4 w p 7, "positrons0_duvernois01_gev.dat" w l 1, "positrons0_duvernois01_gev.dat" u 3:4 w p 8, "positrons0_alcaraz00_gev.dat" w l 1, "positrons0_alcaraz00_gev.dat" u 3:4 w p 9
