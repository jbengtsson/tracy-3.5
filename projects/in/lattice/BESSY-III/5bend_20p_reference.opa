{..en\opa\examples\example-cells\lattices talk 2\test4_5bend_20period.opa}


energy = 2.500000;

    betax   = 0.0594000; alphax  = 0.0000000;
    etax    = 0.0016000; etaxp   = 0.0000000;
    betay   = 15.0000000; alphay  = 0.0000000;
    etay    = 0.0000000; etayp   = 0.0000000;

{----- variables ---------------------------------------------------}

orbitdpp  = 0;
b0a       = 2.39*1.0;
b1a       = -0.205*1.0;
mb1a      = -0.0*1.0;
mb0a      = 2.445-mb1a;

{----- table of elements ---------------------------------------------}

l0  : drift, l = 0.100000, ax = 20.00, ay = 20.00;
l1  : drift, l = 0.100000, ax = 20.00, ay = 20.00;
l2  : drift, l = 0.100000, ax = 20.00, ay = 20.00;
l3  : drift, l = 0.100000, ax = 20.00, ay = 20.00;
ml0 : drift, l = 0.050000, ax = 20.00, ay = 20.00;
ml1 : drift, l = 0.050000, ax = 20.00, ay = 20.00;
ml2 : drift, l = 0.100000, ax = 20.00, ay = 20.00;
ml3 : drift, l = 0.100000, ax = 20.00, ay = 20.00;
ul0 : drift, l = 0.100000, ax = 20.00, ay = 20.00;
ul1 : drift, l = 0.100000, ax = 20.00, ay = 20.00;
ul2 : drift, l = 0.100000, ax = 20.00, ay = 20.00;
ul3 : drift, l = 0.100000, ax = 20.00, ay = 20.00;
ul4 : drift, l = 2.500000, ax = 20.00, ay = 20.00;
o1l : drift, l = 0.040000, ax = 20.00, ay = 20.00;

q0  : quadrupole, l = 0.090000, k = -11.102108, ax = 20.00, ay = 20.00;
q1  : quadrupole, l = 0.140000, k = 11.052113, ax = 20.00, ay = 20.00;
q1a : quadrupole, l = 0.120000, k = 11.335714, ax = 20.00, ay = 20.00;
mq0 : quadrupole, l = 0.100000, k = 0.000000, ax = 20.00, ay = 20.00;
mq1 : quadrupole, l = 0.130000, k = -11.304641, ax = 20.00, ay = 20.00;
mq2 : quadrupole, l = 0.160000, k = 11.684475, ax = 20.00, ay = 20.00;
uq1 : quadrupole, l = 0.090000, k = -9.208116, ax = 20.00, ay = 20.00;
uq2 : quadrupole, l = 0.240000, k = 11.669158, ax = 20.00, ay = 20.00;
uq3 : quadrupole, l = 0.090000, k = -11.061242, ax = 20.00, ay = 20.00;

b0  : bending, l = 0.230000, t = b0a, k = 0.000000, t1 = 0.000000, t2 = b0a,
      ax = 20.00, ay = 20.00;
b1  : bending, l = 0.050000, t = b1a, k = 0.000000, t1 = b1a/2., t2 = b1a/2.,
      ax = 20.00, ay = 20.00;
mb0 : bending, l = 0.230000, t = mb0a, k = 0.000000, t1 = mb0a/2.,
      t2 = mb0a/2., ax = 20.00, ay = 20.00;
mb1 : bending, l = 0.050000, t = mb1a, k = 0.000000, t1 = mb1a/2.,
      t2 = mb1a/2., ax = 20.00, ay = 20.00;

s0  : sextupole, l = 0.060000, k = 721.734906, n = 1, ax = 20.00,
      ay = 20.00;
s0a : sextupole, l = 0.040000, k = 782.633295, n = 1, ax = 20.00,
      ay = 20.00;
s1  : sextupole, l = 0.080000, k = -762.247450, n = 1, ax = 20.00,
      ay = 20.00;
s1a : sextupole, l = 0.080000, k = -87.578394, n = 1, ax = 20.00,
      ay = 20.00;
s3  : sextupole, l = 0.080000, k = -407.911570, n = 1, ax = 20.00,
      ay = 20.00;
s4  : sextupole, l = 0.080000, k = -252.946142, n = 1, ax = 20.00,
      ay = 20.00;
s5  : sextupole, l = 0.080000, k = 503.956264, n = 1, ax = 20.00,
      ay = 20.00;

mbc : opticsmarker, betax = 0.059400, alphax = 0.000000, betay = 15.000000,
      alphay = 0.000000, etax  = 0.001600, etaxp  = 0.000000, etay  = 0.000000,
      etayp  = 0.000000, ax = 20.00, ay = 20.00;


{----- table of segments ---------------------------------------------}

cell    : mbc, b0, l0, q0, l1, s1, l1, q1, l2, b1, l3, s0;
cell1   : mbc, b0, l0, q0, l1, s1, l1, q1a, l2, b1, l3, s0a;
mbacell : -cell, cell;
mcell   : mb0, ml3, mq1, ml2, s1a, ml2, mq2, ml1, ml0, s0a;
mund    : ul0, s3, ul0, uq1, ul1, uq2, ul2, s5, ul2, uq3, ul4;
match   : cell1, -mcell, mund;
arc     : -match, cell, mbacell, -cell, match, nper=20;

{..en\opa\examples\example-cells\lattices talk 2\test4_5bend_20period.opa}
