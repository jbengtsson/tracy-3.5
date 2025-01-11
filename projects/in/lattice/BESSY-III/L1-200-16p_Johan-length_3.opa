{d:\profile\gcf\eigene dateien\opa\tracy\l1-200-16p_johan-length_2.opa}


energy = 2.500000;

    betax   = 2.8114542; alphax  = 0.0000000;
    etax    = -0.0000049; etaxp   = 0.0000000;
    betay   = 3.9672688; alphay  = 0.0000000;
    etay    = 0.0000000; etayp   = 0.0000000;

{----- variables ---------------------------------------------------}

orbitdpp  = 0;
b0a       = 2.455*1.0;
b1a       = -0.205*1.0;
mb1a      = -0.165*1.0;
mb0a      = 2.415*1.0;

{----- table of elements ---------------------------------------------}

l01  : drift, l = 0.549613, ax = 50.00, ay = 50.00;
l02  : drift, l = 0.091515, ax = 50.00, ay = 50.00;
l11  : drift, l = 0.150000, ax = 50.00, ay = 50.00;
l12  : drift, l = 0.150000, ax = 50.00, ay = 50.00;
l2   : drift, l = 0.100000, ax = 50.00, ay = 50.00;
l3   : drift, l = 0.100000, ax = 50.00, ay = 50.00;
ml0  : drift, l = 0.440000, ax = 50.00, ay = 50.00;
ml1  : drift, l = 0.239961, ax = 50.00, ay = 50.00;
ml2  : drift, l = 0.100000, ax = 50.00, ay = 50.00;
ml3  : drift, l = 0.100000, ax = 50.00, ay = 50.00;
ul0  : drift, l = 0.100000, ax = 50.00, ay = 50.00;
ul1  : drift, l = 0.059837, ax = 50.00, ay = 50.00;
ul2  : drift, l = 0.100000, ax = 50.00, ay = 50.00;
ul3  : drift, l = 0.100000, ax = 50.00, ay = 50.00;
ul4  : drift, l = 2.800000, ax = 50.00, ay = 50.00;
o1l  : drift, l = 0.040000, ax = 50.00, ay = 50.00;

uq1  : quadrupole, l = 0.100000, k = -7.721786, ax = 50.00, ay = 50.00;
uq2  : quadrupole, l = 0.120000, k = 7.647096, ax = 50.00, ay = 50.00;
uq3  : quadrupole, l = 0.090000, k = -1.990310, ax = 50.00, ay = 50.00;

b0   : bending, l = 0.450000, t = 2.384952, k = -1.604458, t1 = 1.192476,
       t2 = 1.192476, ax = 50.00, ay = 50.00;
b1   : bending, l = 0.150000, t = -0.206882, k = 7.465899, t1 = -0.103441,
       t2 = -0.103441, ax = 50.00, ay = 50.00;
mb0  : bending, l = 0.450000, t = 2.399116, k = 0.313578, t1 = 1.199558,
       t2 = 1.199558, ax = 50.00, ay = 50.00;
mb1  : bending, l = 0.130000, t = 0.138602, k = 7.911741, t1 = 0.069301,
       t2 = 0.069301, ax = 50.00, ay = 50.00;

s0   : sextupole, l = 0.060000, k = 637.205492, n = 1, ax = 50.00,
       ay = 50.00;
s1   : sextupole, l = 0.080000, k = -668.482462, n = 1, ax = 50.00,
       ay = 50.00;
s3   : sextupole, l = 0.080000, k = 0.000000, n = 5, ax = 50.00, ay = 50.00;
s4   : sextupole, l = 0.080000, k = 0.000000, n = 5, ax = 50.00, ay = 50.00;

mbac : opticsmarker, betax = 0.250000, alphax = 0.000000, betay = 7.400000,
       alphay = 0.000000, etax  = 0.009200, etaxp  = 0.000000,
       etay  = 0.000000, etayp  = 0.000000, ax = 50.00, ay = 50.00;
sf   : opticsmarker, betax = 4.616000, alphax = 0.000000, betay = 2.135760,
       alphay = 0.000000, etax  = 0.043577, etaxp  = -0.000033,
       etay  = 0.000000, etayp  = 0.000000, ax = 50.00, ay = 50.00;

o1   : multipole, n = 1,  k = 0.00000000, ax = 50.00, ay = 50.00;
o2   : multipole, n = 1,  k = 0.00000000, ax = 50.00, ay = 50.00;
o3   : multipole, n = 1,  k = 0.00000000, ax = 50.00, ay = 50.00;


{----- table of segments ---------------------------------------------}

cell0    : mbac, b0, l11, s1, l12, b1, l3, s0, sf;
cell1    : mbac, mb0, ml0, s1, ml1, mb1, ml2, s0, sf;
mbacell  : -cell0, cell0;
mund     : ul1, uq1, o1, uq1, ul0, s3, ul0, uq2, o2, uq2, ul2, s4, ul2, uq3,
           o3, uq3, ul4;
straight : -mund;
match    : -cell1, mund;
arc      : -match, 4*mbacell, match;

{d:\profile\gcf\eigene dateien\opa\tracy\l1-200-16p_johan-length_2.opa}
