{..a\examples\example-cells\20periods\5bend_20p_longbend-tg-tgrb_sym1.opa}


energy = 2.500000;

    betax   = 1.5000000; alphax  = 0.0000000;
    etax    = -0.0000756; etaxp   = 0.0000000;
    betay   = 1.5000000; alphay  = 0.0000000;
    etay    = 0.0000000; etayp   = 0.0000000;

{----- variables ---------------------------------------------------}

orbitdpp  = 0;
b0a       = 2.3*1.0;
b1a       = -0.205*1.0;
mb1a      = -0.0*1.0;
mb0a      = 2.715-mb1a;

{----- table of elements ---------------------------------------------}

l0   : drift, l = 0.100000, ax = 300.00, ay = 300.00;
l1   : drift, l = 0.100000, ax = 300.00, ay = 300.00;
l2   : drift, l = 0.100000, ax = 300.00, ay = 300.00;
l3   : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ml1  : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ml2  : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ml3  : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ul0  : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ul1  : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ul2  : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ul3  : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ul4  : drift, l = 2.500000, ax = 300.00, ay = 300.00;
o1l  : drift, l = 0.040000, ax = 50.00, ay = 50.00;

q0   : quadrupole, l = 0.090000, k = -10.424964, ax = 300.00, ay = 300.00;
q1   : quadrupole, l = 0.130000, k = 11.138731, ax = 300.00, ay = 300.00;
mq1  : quadrupole, l = 0.090000, k = -10.424964, ax = 300.00, ay = 300.00;
mq2  : quadrupole, l = 0.090000, k = 11.138731, ax = 300.00, ay = 300.00;
uq1  : quadrupole, l = 0.090000, k = -7.833987, ax = 300.00, ay = 300.00;
uq2  : quadrupole, l = 0.220000, k = 11.843834, ax = 300.00, ay = 300.00;
uq3  : quadrupole, l = 0.090000, k = -10.534429, ax = 300.00, ay = 300.00;
uq4  : quadrupole, l = 0.320000, k = 0.000000, ax = 300.00, ay = 300.00;

b0   : bending, l = 0.400000, t = b0a, k = -2.225569, t1 = 0.000000, t2 = b0a,
       ax = 50.00, ay = 300.00;
b1   : bending, l = 0.150000, t = b1a, k = 9.524157, t1 = b1a/2., t2 = b1a/2.,
       ax = 50.00, ay = 300.00;
b1a  : bending, l = 0.150000, t = b1a, k = 11.493091, t1 = b1a/2., t2 = b1a/2.,
       ax = 50.00, ay = 300.00;
mb0  : bending, l = 0.400000, t = mb0a, k = -3.923682, t1 = mb0a/2.,
       t2 = mb0a/2., ax = 50.00, ay = 300.00;
mb1  : bending, l = 0.150000, t = mb1a, k = 10.531481, t1 = mb1a/2.,
       t2 = mb1a/2., ax = 50.00, ay = 300.00;

sd1  : sextupole, l = 0.080000, k = 123.125730, n = 5, ax = 300.00,
       ay = 300.00;
sd2  : sextupole, l = 0.080000, k = -1635.611193, n = 5, ax = 300.00,
       ay = 300.00;
sd3  : sextupole, l = 0.080000, k = -759.204665, n = 5, ax = 300.00,
       ay = 300.00;
sd4  : sextupole, l = 0.080000, k = -1333.859875, n = 5, ax = 300.00,
       ay = 300.00;
sf1  : sextupole, l = 0.040000, k = 1242.074776, n = 5, ax = 300.00,
       ay = 300.00;
sf2  : sextupole, l = 0.040000, k = 1555.462511, n = 5, ax = 300.00,
       ay = 300.00;
sh1  : sextupole, l = 0.080000, k = 0.000000, n = 5, ax = 300.00,
       ay = 300.00;
sh2  : sextupole, l = 0.080000, k = 0.000000, n = 5, ax = 300.00,
       ay = 300.00;

mbac : opticsmarker, betax = 0.103000, alphax = 0.000000, betay = 8.000000,
       alphay = 0.000000, etax  = 0.006000, etaxp  = 0.000000,
       etay  = 0.000000, etayp  = 0.000000, ax = 50.00, ay = 50.00;


{----- table of segments ---------------------------------------------}

cell    : mbac, b0, l1, sd1, l1, b1, l3, sf1;
cell1   : mbac, b0, l1, sd2, l1, b1, l3, sf1;
cell2   : mbac, b0, l1, sd3, l1, b1a, l3, sf2;
mbacell : -cell, cell;
mcell   : mb0, ml2, sd4, ml2, mb1, ml1, sf2;
mund    : ul0, sh1, ul0, uq1, ul1, uq2, ul2, sh2, ul2, uq3, ul3, ul4;
match   : cell2, -mcell, mund;
arc     : -match, cell1, mbacell, -cell1, match, nper=20;

{..a\examples\example-cells\20periods\5bend_20p_longbend-tg-tgrb_sym1.opa}
