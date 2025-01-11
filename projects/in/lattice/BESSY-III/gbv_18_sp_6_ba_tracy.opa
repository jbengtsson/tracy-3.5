{d:\profile\gcf\eigene dateien\opa\tracy\gbv_18_sp_6_ba_tracy.opa}


energy = 2.500000;

    betax   = 2.4602406; alphax  = 0.0000000;
    etax    = -0.0003543; etaxp   = 0.0000000;
    betay   = 3.6407683; alphay  = 0.0000000;
    etay    = 0.0000000; etayp   = 0.0000000;

{----- variables ---------------------------------------------------}


{----- table of elements ---------------------------------------------}

l1    : drift, l = 0.100000, ax = 50.00, ay = 50.00;
l2    : drift, l = 0.100000, ax = 50.00, ay = 50.00;
l3    : drift, l = 0.100000, ax = 50.00, ay = 50.00;
l4    : drift, l = 0.100000, ax = 50.00, ay = 50.00;
l5    : drift, l = 0.100000, ax = 50.00, ay = 50.00;
l6    : drift, l = 0.100000, ax = 50.00, ay = 50.00;
l7    : drift, l = 0.100000, ax = 50.00, ay = 50.00;
l8    : drift, l = 0.100000, ax = 50.00, ay = 50.00;
l9    : drift, l = 2.500000, ax = 50.00, ay = 50.00;

q1    : quadrupole, l = 0.160000, k = -7.786219, ax = 15.00, ay = 15.00;
q2    : quadrupole, l = 0.260000, k = 8.883494, ax = 15.00, ay = 15.00;
q3    : quadrupole, l = 0.100000, k = -6.460781, ax = 15.00, ay = 15.00;

b1e   : bending, l = 0.699987, t = 2.428817, k = -1.201667, t1 = 1.214408,
        t2 = 1.214408, ax = 15.00, ay = 15.00;
b1c   : bending, l = 0.400000, t = 2.207636, k = -2.231373, t1 = 0.000000,
        t2 = 2.207636, ax = 15.00, ay = 15.00;
b2    : bending, l = 0.130000, t = -0.251872, k = 10.557202, t1 = -0.125936,
        t2 = -0.125936, ax = 15.00, ay = 15.00;

s1    : sextupole, l = 0.060000, k = 875.787072, n = 1, ax = 15.00,
        ay = 15.00;
s2    : sextupole, l = 0.080000, k = -818.312905, n = 1, ax = 15.00,
        ay = 15.00;
s3    : sextupole, l = 0.080000, k = 0.000000, n = 1, ax = 15.00,
        ay = 15.00;
s4    : sextupole, l = 0.080000, k = 0.000000, n = 1, ax = 15.00,
        ay = 15.00;

mbac  : opticsmarker, betax = 0.215313, alphax = 0.001151, betay = 6.489271,
        alphay = -0.001595, etax  = 0.006844, etaxp  = 0.000281,
        etay  = 0.000000, etayp  = 0.000000, ax = 50.00, ay = 50.00;
s1c   : opticsmarker, betax = 3.710672, alphax = -0.000731, betay = 1.907478,
        alphay = -0.001071, etax  = 0.032162, etaxp  = 0.000108,
        etay  = 0.000000, etayp  = 0.000000, ax = 50.00, ay = 50.00;
dispe : opticsmarker, betax = 0.738712, alphax = 1.899367, betay = 8.264245,
        alphay = -0.143498, etax  = -0.000075, etaxp  = 0.000435,
        etay  = 0.000000, etayp  = 0.000000, ax = 50.00, ay = 50.00;

o1    : multipole, n = 4,  k = 793.90419417, ax = 15.00, ay = 15.00;
o2    : multipole, n = 4,  k = 0.00000000, ax = 15.00, ay = 15.00;
o3    : multipole, n = 4,  k = 2194.38937212, ax = 15.00, ay = 15.00;
o4    : multipole, n = 4,  k = -3056.12723328, ax = 15.00, ay = 15.00;


{----- table of segments ---------------------------------------------}

cell_hc : mbac, b1c, l1, s2, l2, b2, l3, s1, s1c;
cell_u  : -cell_hc, cell_hc;
cell_he : s1c, s1, l3, b2, l2, s2, l1, b1e, dispe;
cell_m  : l4, s3, o1, l5, q1, l6, o2, q2, l7, s4, o3, l8, q3, o4, l9;
match   : cell_he, cell_m;
achr    : -match, 4*cell_u, match;
ring    : achr, nper=18;

{d:\profile\gcf\eigene dateien\opa\tracy\gbv_18_sp_6_ba_tracy.opa}
