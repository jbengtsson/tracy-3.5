{d:\profile\gcf\eigene dateien\opa\tracy\l1-200-16p_johan-phase.opa}


energy = 2.500000;

    betax   = 2.4882652; alphax  = 0.0000000;
    etax    = -0.0000024; etaxp   = 0.0000000;
    betay   = 4.9131004; alphay  = 0.0000000;
    etay    = 0.0000000; etayp   = 0.0000000;

{----- variables ---------------------------------------------------}

orbitdpp  = 0;
b0a       = 2.455*1.0;
b1a       = -0.205*1.0;
mb1a      = -0.165*1.0;
mb0a      = 2.415*1.0;

{----- table of elements ---------------------------------------------}

l01  : drift, l = 0.549613, ax = 300.00, ay = 300.00;
l02  : drift, l = 0.091515, ax = 300.00, ay = 300.00;
l11  : drift, l = 0.150000, ax = 300.00, ay = 300.00;
l12  : drift, l = 0.150000, ax = 300.00, ay = 300.00;
l2   : drift, l = 0.100000, ax = 300.00, ay = 300.00;
l3   : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ml0  : drift, l = 0.503876, ax = 300.00, ay = 300.00;
ml1  : drift, l = 0.200000, ax = 300.00, ay = 300.00;
ml2  : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ml3  : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ul0  : drift, l = 0.092888, ax = 300.00, ay = 300.00;
ul1  : drift, l = 0.340142, ax = 300.00, ay = 300.00;
ul2  : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ul3  : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ul4  : drift, l = 2.800000, ax = 300.00, ay = 300.00;

uq1  : quadrupole, l = 0.090000, k = -7.498685, ax = 300.00, ay = 300.00;
uq2  : quadrupole, l = 0.120000, k = 7.286132, ax = 300.00, ay = 300.00;
uq3  : quadrupole, l = 0.090000, k = -1.852174, ax = 300.00, ay = 300.00;

b0   : bending, l = 0.450000, t = b0a, k = -1.608925, t1 = 0.000000, t2 = b0a,
       ax = 50.00, ay = 300.00;
b1   : bending, l = 0.150000, t = b1a, k = 7.470935, t1 = b1a/2., t2 = b1a/2.,
       ax = 50.00, ay = 300.00;
mb0  : bending, l = 0.450000, t = mb0a, k = 0.000000, t1 = mb0a/2.,
       t2 = mb0a/2., ax = 50.00, ay = 300.00;
mb1  : bending, l = 0.120000, t = mb1a, k = 7.403998, t1 = mb1a/2.,
       t2 = mb1a/2., ax = 50.00, ay = 300.00;

s0   : sextupole, l = 0.060000, k = 386.044847, n = 5, ax = 300.00,
       ay = 300.00;
s1   : sextupole, l = 0.080000, k = -477.524900, n = 5, ax = 300.00,
       ay = 300.00;
s3   : sextupole, l = 0.080000, k = 0.000000, n = 5, ax = 50.00, ay = 50.00;
s4   : sextupole, l = 0.080000, k = 683.534611, n = 5, ax = 300.00,
       ay = 300.00;

mbac : opticsmarker, betax = 0.250000, alphax = 0.000000, betay = 7.300000,
       alphay = 0.000000, etax  = 0.009600, etaxp  = 0.000000,
       etay  = 0.000000, etayp  = 0.000000, ax = 50.00, ay = 50.00;
sf   : opticsmarker, betax = 4.622240, alphax = 0.000000, betay = 2.104470,
       alphay = 0.000000, etax  = 0.045076, etaxp  = -0.000017,
       etay  = 0.000000, etayp  = 0.000000, ax = 50.00, ay = 50.00;

o1   : multipole, n = 4,  k = 0.00000000, ax = 50.00, ay = 50.00;
o2   : multipole, n = 4,  k = 0.00000000, ax = 50.00, ay = 50.00;
o3   : multipole, n = 4,  k = 0.00000000, ax = 50.00, ay = 50.00;


{----- table of segments ---------------------------------------------}

cell     : mbac, b0, l11, s1, l12, b1, l3, s0, sf;
cell1    : mbac, mb0, ml0, s1, ml1, mb1, ml2, s0, sf;
mbacell  : -cell, cell;
mund     : ul1, uq1, o1, uq1, ul0, s3, ul0, uq2, o2, uq2, ul2, s4, ul2, uq3,
           o3, uq3, ul4;
straight : -mund;
match    : -cell1, mund;
arc      : -match, 4*mbacell, match, nper=16;

{d:\profile\gcf\eigene dateien\opa\tracy\l1-200-16p_johan-phase.opa}
