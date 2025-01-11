{d:\lattices\b3_lat\20249999_final\january2024-test2a-425°rb28-wb.opa}


energy = 2.500000;
rotinv = 0;
    betax   = 3.1289737; alphax  = 0.0000000;
    etax    = 0.0000004; etaxp   = 0.0000000;
    betay   = 3.5018782; alphay  = 0.0000000;
    etay    = 0.0000000; etayp   = 0.0000000;

{----- variables ---------------------------------------------------}

rba   = -0.28;
b1a   = 4.25-2.0*rba;
mwba  = 0.2;
mb1a  = 2.75-mwba;

{----- table of elements ---------------------------------------------}

lo     : drift, l = 0.025000, ax = 300.00, ay = 300.00;
l1     : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ml1    : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ml2    : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ml3    : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ml4    : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ul1    : drift, l = 0.150000, ax = 300.00, ay = 300.00;
ul2    : drift, l = 0.150000, ax = 300.00, ay = 300.00;
ul3    : drift, l = 0.130000, ax = 300.00, ay = 300.00;
ul4    : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ul5    : drift, l = 2.800000, ax = 300.00, ay = 300.00;

qd     : quadrupole, l = 0.130000, k = -9.178248, ax = 9.00, ay = 9.00;
uq1    : quadrupole, l = 0.100000, k = 8.677909, ax = 9.00, ay = 9.00;
uq2    : quadrupole, l = 0.200000, k = -9.526291, ax = 9.00, ay = 9.00;
uq3    : quadrupole, l = 0.260000, k = 9.052253, ax = 9.00, ay = 9.00;
uq4    : quadrupole, l = 0.100000, k = -8.509481, ax = 9.00, ay = 9.00;
mb2    : quadrupole, l = 0.180000, k = 9.365333, ax = 50.00, ay = 300.00;

b1     : bending, l = 0.550000, t = b1a/2, k = 0.000000, t1 = 0.000000,
         t2 = b1a/2.0, ax = 50.00, ay = 300.00;
mb1    : bending, l = 0.650000, t = mb1a, k = 0.000000, t1 = mb1a/2.0,
         t2 = mb1a/2.0, ax = 50.00, ay = 300.00;

sd     : sextupole, l = 0.050000, k = -236.499033, n = 1, ax = 300.00,
         ay = 300.00;
sf     : sextupole, l = 0.050000, k = 243.292480, n = 1, ax = 300.00,
         ay = 300.00;
sd5    : sextupole, l = 0.050000, k = -239.967748, n = 1, ax = 300.00,
         ay = 300.00;
sd4    : sextupole, l = 0.050000, k = -211.299921, n = 1, ax = 300.00,
         ay = 300.00;
sd3    : sextupole, l = 0.050000, k = -215.407243, n = 1, ax = 300.00,
         ay = 300.00;
sd2    : sextupole, l = 0.050000, k = -200.093085, n = 1, ax = 300.00,
         ay = 300.00;
sd1    : sextupole, l = 0.050000, k = -226.546244, n = 1, ax = 300.00,
         ay = 300.00;
sf3    : sextupole, l = 0.050000, k = 230.944691, n = 1, ax = 300.00,
         ay = 300.00;
sf2    : sextupole, l = 0.050000, k = 204.087886, n = 1, ax = 300.00,
         ay = 300.00;
sf1    : sextupole, l = 0.050000, k = 190.137366, n = 1, ax = 300.00,
         ay = 300.00;
shx    : sextupole, l = 0.000100, k = 0.000000, n = 1, ax = 50.00,
         ay = 50.00;
shy    : sextupole, l = 0.000100, k = 0.000000, n = 1, ax = 50.00,
         ay = 50.00;

om_sf  : opticsmarker, betax = 5.786310, alphax = 0.000000, betay = 2.828970,
         alphay = 0.000000, etax  = 0.058092, etaxp  = 0.000000,
         etay  = 0.000000, etayp  = 0.000000, ax = 50.00, ay = 50.00;
om_mb1 : opticsmarker, betax = 1.103270, alphax = -1.889940, betay = 5.998490,
         alphay = 0.114930, etax  = 0.000000, etaxp  = 0.000000,
         etay  = 0.000000, etayp  = 0.000000, ax = 50.00, ay = 50.00;
om_c   : opticsmarker, betax = 0.410624, alphax = 0.000000, betay = 5.370696,
         alphay = 0.000000, etax  = 0.009545, etaxp  = 0.000000,
         etay  = 0.000000, etayp  = 0.000000, ax = 50.00, ay = 50.00;

mqd    : combined, l = 0.150000, t = mwba, k = -8.701492, t1 = mwba/2.0,
         t2 = mwba/2.0, ax = 9.00, ay = 9.00;
rb     : combined, l = 0.170000, t = rba, k = 8.642872, t1 = rba/2.0,
         t2 = rba/2.0, ax = 50.00, ay = 300.00;

oxx    : multipole, n = 4,  k = 320.00000000, ax = 50.00, ay = 50.00;
oyy    : multipole, n = 4,  k = 0.00000000, ax = 50.00, ay = 50.00;


{----- table of segments ---------------------------------------------}

ucell : om_sf, sf, l1, rb, l1, qd, l1, sd, sd, l1, -b1, om_c, b1, l1, sd, sd,
        l1, qd, l1, rb, l1, sf, om_sf;
uc21  : om_sf, sf1, l1, rb, l1, qd, l1, sd1, sd1, l1, -b1, om_c, b1, l1, sd2,
        sd2, l1, qd, l1, rb, l1, sf2, om_sf;
uc32  : om_sf, sf2, l1, rb, l1, qd, l1, sd3, sd3, l1, -b1, om_c, b1, l1, sd4,
        sd4, l1, qd, l1, rb, l1, sf3, om_sf;
dcell : om_sf, sf3, sf3, ml1, lo, oxx, lo, ml1, mb2, ml2, mqd, ml3, sd5,
        sd5, sd5, ml4, mb1, om_mb1;
arc   : -dcell, -uc32, -uc21, om_c, uc21, uc32, dcell;
mund  : om_mb1, ul1, uq1, ul2, uq2, shy, ul3, uq3, shx, ul4, uq4, ul5;
sec   : -mund, arc, mund;
sech  : om_c, uc21, uc32, dcell, mund;
sec16 : -mund, arc, mund, nper=16;

{d:\lattices\b3_lat\20249999_final\january2024-test2a-425°rb28-wb.opa}
