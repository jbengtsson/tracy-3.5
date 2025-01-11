{..le-cells\lattices talk 2\test4_dipl_variation-2p19deg-tg-pure-tgrb.opa}


energy = 2.500000;

    betax   = 1.6295376; alphax  = 0.0000000;
    etax    = -0.0004301; etaxp   = 0.0000000;
    betay   = 1.6088879; alphay  = 0.0000000;
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
ml0  : drift, l = 0.050000, ax = 300.00, ay = 300.00;
ml1  : drift, l = 0.050000, ax = 300.00, ay = 300.00;
ml2  : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ml3  : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ul0  : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ul1  : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ul2  : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ul3  : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ul4  : drift, l = 2.500000, ax = 300.00, ay = 300.00;
o1l  : drift, l = 0.040000, ax = 50.00, ay = 50.00;

q0   : quadrupole, l = 0.090000, k = -10.424964, ax = 300.00, ay = 300.00;
q0a  : quadrupole, l = 0.090000, k = -11.000000, ax = 300.00, ay = 300.00;
q1   : quadrupole, l = 0.130000, k = 11.138731, ax = 300.00, ay = 300.00;
q1a  : quadrupole, l = 0.090000, k = 3.348660, ax = 300.00, ay = 300.00;
mq0  : quadrupole, l = 0.100000, k = 0.000000, ax = 300.00, ay = 300.00;
mq1  : quadrupole, l = 0.090000, k = -4.393642, ax = 300.00, ay = 300.00;
mq2  : quadrupole, l = 0.090000, k = 11.616342, ax = 300.00, ay = 300.00;
uq1  : quadrupole, l = 0.160000, k = -11.036002, ax = 300.00, ay = 300.00;
uq2  : quadrupole, l = 0.260000, k = 11.463513, ax = 300.00, ay = 300.00;
uq3  : quadrupole, l = 0.100000, k = -11.313212, ax = 300.00, ay = 300.00;

b0   : bending, l = 0.400000, t = b0a, k = -2.251588, t1 = 0.000000, t2 = b0a,
       ax = 50.00, ay = 300.00;
b1   : bending, l = 0.130000, t = b1a, k = 11.105635, t1 = b1a/2., t2 = b1a/2.,
       ax = 50.00, ay = 300.00;
mb0  : bending, l = 0.230000, t = mb0a, k = 0.000000, t1 = mb0a/2.,
       t2 = mb0a/2., ax = 50.00, ay = 300.00;
mb1  : bending, l = 0.050000, t = mb1a, k = 0.000000, t1 = mb1a/2.,
       t2 = mb1a/2., ax = 50.00, ay = 300.00;

s0   : sextupole, l = 0.060000, k = 1233.318357, n = 1, ax = 300.00,
       ay = 300.00;
s0a  : sextupole, l = 0.040000, k = 753.588285, n = 1, ax = 300.00,
       ay = 300.00;
s1   : sextupole, l = 0.080000, k = -1182.973681, n = 1, ax = 300.00,
       ay = 300.00;
s1a  : sextupole, l = 0.080000, k = 785.769128, n = 1, ax = 300.00,
       ay = 300.00;
s3   : sextupole, l = 0.080000, k = -454.323306, n = 1, ax = 300.00,
       ay = 300.00;
s4   : sextupole, l = 0.080000, k = -252.946142, n = 1, ax = 300.00,
       ay = 300.00;
s5   : sextupole, l = 0.080000, k = 226.205241, n = 1, ax = 300.00,
       ay = 300.00;

mbac : opticsmarker, betax = 0.103000, alphax = 0.000000, betay = 8.000000,
       alphay = 0.000000, etax  = 0.006000, etaxp  = 0.000000,
       etay  = 0.000000, etayp  = 0.000000, ax = 50.00, ay = 50.00;


{----- table of segments ---------------------------------------------}

cell    : mbac, b0, l1, s1, l1, b1, l3, s0;
cell1   : mbac, b0, l1, s1, l1, q1a, l2, b1, l3, s0a;
mbacell : -cell, cell;
mcell   : mb0, ml3, mq1, ml2, s1a, ml2, mq2, ml1, ml0, s0a;
mund    : ul0, s3, ul0, uq1, ul1, uq2, ul2, s5, ul2, uq3, ul4;
match   : cell1, -mcell, mund;
arc     : -match, cell, mbacell, -cell, match, nper=20;

{..le-cells\lattices talk 2\test4_dipl_variation-2p19deg-tg-pure-tgrb.opa}
