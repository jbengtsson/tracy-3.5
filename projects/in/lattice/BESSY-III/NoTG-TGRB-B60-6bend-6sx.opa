{..examples\example-cells\talk 1 low gradient\notg-tgrb-b60-6bend-6sx.opa}


energy = 2.500000;

    betax   = 2.6142733; alphax  = 0.0000000;
    etax    = 0.0002142; etaxp   = 0.0000000;
    betay   = 2.4368208; alphay  = 0.0000000;
    etay    = 0.0000000; etayp   = 0.0000000;

{----- variables ---------------------------------------------------}

orbitdpp  = 0;
b0a       = 2.205*1.0;
b1a       = -0.25*1.0;
mb0a      = 3.43*1.0;
mb1a      = -0.0*1.0;

{----- table of elements ---------------------------------------------}

l0  : drift, l = 0.100000, ax = 300.00, ay = 300.00;
l1  : drift, l = 0.100000, ax = 300.00, ay = 300.00;
l2  : drift, l = 0.100000, ax = 300.00, ay = 300.00;
l3  : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ml0 : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ml1 : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ml2 : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ml3 : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ul0 : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ul1 : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ul2 : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ul3 : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ul4 : drift, l = 2.800000, ax = 300.00, ay = 300.00;

q0  : quadrupole, l = 0.110000, k = -7.702058, ax = 9.00, ay = 9.00;
q1  : quadrupole, l = 0.180000, k = 7.965777, ax = 9.00, ay = 9.00;
q1a : quadrupole, l = 0.090000, k = 4.351168, ax = 9.00, ay = 9.00;
mq0 : quadrupole, l = 0.100000, k = 0.000000, ax = 9.00, ay = 9.00;
mq1 : quadrupole, l = 0.080000, k = -8.878971, ax = 9.00, ay = 9.00;
mq2 : quadrupole, l = 0.110000, k = 8.617091, ax = 9.00, ay = 9.00;
uq1 : quadrupole, l = 0.100000, k = -6.174755, ax = 9.00, ay = 9.00;
uq2 : quadrupole, l = 0.280000, k = 7.731768, ax = 9.00, ay = 9.00;
uq3 : quadrupole, l = 0.090000, k = -7.696539, ax = 9.00, ay = 9.00;

b0  : bending, l = 0.600000, t = b0a, k = 0.000000, t1 = 0.000000, t2 = b0a,
      ax = 50.00, ay = 300.00;
b1  : bending, l = 0.140000, t = b1a, k = 8.081387, t1 = b1a/2., t2 = b1a/2.,
      ax = 50.00, ay = 300.00;
b1a : bending, l = 0.220000, t = b1a, k = 7.630257, t1 = b1a/2., t2 = b1a/2.,
      ax = 20.00, ay = 9.00;
mb0 : bending, l = 0.600000, t = mb0a, k = -2.237166, t1 = mb0a/2.,
      t2 = mb0a/2., ax = 50.00, ay = 300.00;
mb1 : bending, l = 0.000000, t = mb1a, k = 0.000000, t1 = mb1a/2.,
      t2 = mb1a/2., ax = 50.00, ay = 300.00;

s0  : sextupole, l = 0.040000, k = 505.109148, n = 1, ax = 300.00,
      ay = 300.00;
s0a : sextupole, l = 0.040000, k = 1023.621439, n = 1, ax = 300.00,
      ay = 300.00;
s1  : sextupole, l = 0.080000, k = -346.018550, n = 1, ax = 300.00,
      ay = 300.00;
s1a : sextupole, l = 0.080000, k = -734.288733, n = 1, ax = 300.00,
      ay = 300.00;
s3  : sextupole, l = 0.080000, k = 398.901079, n = 1, ax = 300.00,
      ay = 300.00;
s4  : sextupole, l = 0.080000, k = 37.389977, n = 1, ax = 300.00,
      ay = 300.00;

b0c : opticsmarker, betax = 0.395000, alphax = 0.000000, betay = 6.100000,
      alphay = 0.000000, etax  = 0.009800, etaxp  = 0.000000, etay  = 0.000000,
      etayp  = 0.000000, ax = 50.00, ay = 50.00;
s0c : opticsmarker, betax = 7.666260, alphax = 0.000000, betay = 3.526500,
      alphay = 0.000000, etax  = 0.052420, etaxp  = -0.000011,
      etay  = 0.000000, etayp  = 0.000000, ax = 50.00, ay = 50.00;


{----- table of segments ---------------------------------------------}

cell    : b0c, b0, l0, q0, l1, s1, l1, b1, l3, s0;
cell1   : b0c, b0, l0, q0, l1, s1, l1, b1a, l3, s0a, s0c;
mbacell : -cell, cell;
mcell   : mb0, ml2, s1a, ml2, mq2, ml1, s0a;
mund    : ul0, s3, ul0, uq1, ul1, uq2, ul2, s4, ul3, uq3, ul4;
match   : cell1, -mcell, mund;
arc     : -match, cell, mbacell, mbacell, -cell, match, nper=16;

{..examples\example-cells\talk 1 low gradient\notg-tgrb-b60-6bend-6sx.opa}
