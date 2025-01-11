{d:\profile\gcf\eigene dateien\opa\tracy\jbgv_20_sp_5ba-2.opa}


energy = 2.500000;

    betax   = 2.9759918; alphax  = 0.0000000;
    etax    = -0.0000799; etaxp   = 0.0000000;
    betay   = 2.9171600; alphay  = 0.0000000;
    etay    = 0.0000000; etayp   = 0.0000000;

{----- variables ---------------------------------------------------}


{----- table of elements ---------------------------------------------}

l1   : drift, l = 0.100000, ax = 50.00, ay = 50.00;
l2   : drift, l = 0.100000, ax = 50.00, ay = 50.00;
l3   : drift, l = 0.100000, ax = 50.00, ay = 50.00;
l4   : drift, l = 0.100000, ax = 50.00, ay = 50.00;
l5   : drift, l = 0.100000, ax = 50.00, ay = 50.00;
l6   : drift, l = 0.100000, ax = 50.00, ay = 50.00;
l7   : drift, l = 0.100000, ax = 50.00, ay = 50.00;
l8   : drift, l = 0.100000, ax = 50.00, ay = 50.00;
l9   : drift, l = 2.500000, ax = 50.00, ay = 50.00;

q1   : quadrupole, l = 0.160000, k = -7.813983, ax = 50.00, ay = 50.00;
q2   : quadrupole, l = 0.260000, k = 8.955781, ax = 50.00, ay = 50.00;
q3   : quadrupole, l = 0.100000, k = -6.673491, ax = 50.00, ay = 50.00;

b1c  : bending, l = 0.400000, t = 2.514812, k = -2.233492, t1 = 0.000000,
       t2 = 2.514812, ax = 50.00, ay = 50.00;
b1e  : bending, l = 0.665392, t = 2.764584, k = -1.309446, t1 = 1.382292,
       t2 = 1.382292, ax = 50.00, ay = 50.00;
b2   : bending, l = 0.130000, t = -0.327255, k = 10.689206, t1 = -0.163627,
       t2 = -0.163627, ax = 50.00, ay = 50.00;

s1   : sextupole, l = 0.060000, k = 850.992699, n = 1, ax = 50.00,
       ay = 50.00;
s2   : sextupole, l = 0.080000, k = -805.517729, n = 1, ax = 50.00,
       ay = 50.00;
s3   : sextupole, l = 0.080000, k = 0.000000, n = 1, ax = 50.00, ay = 50.00;
s4   : sextupole, l = 0.080000, k = 0.000000, n = 1, ax = 50.00, ay = 50.00;

mbac : opticsmarker, betax = 0.193000, alphax = 0.000000, betay = 6.800000,
       alphay = 0.000000, etax  = 0.006500, etaxp  = 0.000000,
       etay  = 0.000000, etayp  = 0.000000, ax = 50.00, ay = 50.00;
s1c  : opticsmarker, betax = 4.023190, alphax = 0.000960, betay = 1.989550,
       alphay = 0.001118, etax  = 0.035184, etaxp  = 0.000067,
       etay  = 0.000000, etayp  = 0.000000, ax = 50.00, ay = 50.00;


{----- table of segments ---------------------------------------------}

cell_hc : mbac, b1c, l1, s2, l2, b2, l3, s1, s1c;
cell_u  : -cell_hc, cell_hc;
cell_he : s1c, s1, l3, b2, l2, s2, l1, b1e;
cell_m  : l4, s3, l5, q1, l6, q2, l7, s4, l8, q3, l9;
match   : cell_he, cell_m;
line    : -match, 3*cell_u, match;
achr    : line, nper=20;

{d:\profile\gcf\eigene dateien\opa\tracy\jbgv_20_sp_5ba-2.opa}
