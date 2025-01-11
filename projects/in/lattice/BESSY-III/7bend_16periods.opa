{..igene dateien\opa\examples\example-cells\16periods\7bend_16periods.opa}


energy = 2.500000;

    betax   = 1.5000000; alphax  = 0.0000000;
    etax    = 0.0000695; etaxp   = 0.0000000;
    betay   = 1.5000000; alphay  = 0.0000000;
    etay    = 0.0000000; etayp   = 0.0000000;

{----- variables ---------------------------------------------------}

orbitdpp  = 0;
b0a       = 2.0*1.0;
b1a       = -0.165*1.0;
mb0a      = 2.075*1.0;
mb1a      = -0.0*1.0;

{----- table of elements ---------------------------------------------}

l0   : drift, l = 0.100000, ax = 300.00, ay = 300.00;
l1   : drift, l = 0.100000, ax = 300.00, ay = 300.00;
l2   : drift, l = 0.100000, ax = 300.00, ay = 300.00;
l3   : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ml0  : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ml1  : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ml2  : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ml3  : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ul0  : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ul1  : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ul2  : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ul3  : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ul4  : drift, l = 2.500000, ax = 300.00, ay = 300.00;

q0   : quadrupole, l = 0.090000, k = -11.195751, ax = 300.00, ay = 300.00;
q0a  : quadrupole, l = 0.090000, k = -11.188811, ax = 300.00, ay = 300.00;
q1   : quadrupole, l = 0.140000, k = 11.167885, ax = 300.00, ay = 300.00;
q1a  : quadrupole, l = 0.130000, k = 10.480515, ax = 300.00, ay = 300.00;
mq1  : quadrupole, l = 0.130000, k = -11.123916, ax = 300.00, ay = 300.00;
mq2  : quadrupole, l = 0.160000, k = 11.748542, ax = 300.00, ay = 300.00;
uq1  : quadrupole, l = 0.090000, k = -9.793062, ax = 300.00, ay = 300.00;
uq2  : quadrupole, l = 0.240000, k = 11.814885, ax = 300.00, ay = 300.00;
uq3  : quadrupole, l = 0.100000, k = -10.046458, ax = 300.00, ay = 300.00;

b0   : bending, l = 0.220000, t = b0a, k = 0.000000, t1 = 0.000000, t2 = b0a,
       ax = 50.00, ay = 300.00;
b1   : bending, l = 0.050000, t = b1a, k = 0.000000, t1 = b1a/2., t2 = b1a/2.,
       ax = 50.00, ay = 300.00;
b1a  : bending, l = 0.050000, t = b1a, k = 0.000000, t1 = b1a/2., t2 = b1a/2.,
       ax = 50.00, ay = 300.00;
mb0  : bending, l = 0.220000, t = mb0a, k = 0.000000, t1 = mb0a/2.,
       t2 = mb0a/2., ax = 50.00, ay = 300.00;
mb1  : bending, l = 0.000000, t = mb1a, k = 0.000000, t1 = mb1a/2.,
       t2 = mb1a/2., ax = 50.00, ay = 300.00;

sf0  : sextupole, l = 0.040000, k = 774.682337, n = 1, ax = 300.00,
       ay = 300.00;
sf1  : sextupole, l = 0.040000, k = 1802.110305, n = 1, ax = 300.00,
       ay = 300.00;
sf2  : sextupole, l = 0.040000, k = 667.060645, n = 1, ax = 300.00,
       ay = 300.00;
sf3  : sextupole, l = 0.040000, k = 488.766869, n = 1, ax = 300.00,
       ay = 300.00;
sd0  : sextupole, l = 0.080000, k = 63.145599, n = 1, ax = 300.00,
       ay = 300.00;
sd1  : sextupole, l = 0.080000, k = -1345.621471, n = 1, ax = 300.00,
       ay = 300.00;
sd2  : sextupole, l = 0.080000, k = -1117.362064, n = 1, ax = 300.00,
       ay = 300.00;
sd3  : sextupole, l = 0.080000, k = -1067.142576, n = 1, ax = 300.00,
       ay = 300.00;
sd4  : sextupole, l = 0.080000, k = -488.929576, n = 1, ax = 300.00,
       ay = 300.00;
sd5  : sextupole, l = 0.080000, k = -108.533722, n = 1, ax = 300.00,
       ay = 300.00;
sh1  : sextupole, l = 0.080000, k = -305.670767, n = 1, ax = 300.00,
       ay = 300.00;
sh2  : sextupole, l = 0.080000, k = 28.346071, n = 1, ax = 300.00,
       ay = 300.00;

mbac : opticsmarker, betax = 0.059400, alphax = 0.000000, betay = 15.000000,
       alphay = 0.000000, etax  = 0.001300, etaxp  = 0.000000,
       etay  = 0.000000, etayp  = 0.000000, ax = 50.00, ay = 50.00;


{----- table of segments ---------------------------------------------}

cella   : mbac, b0, l0, q0, l1, sd0, l1, q1, l2, b1, l3, sf0;
cellb1  : mbac, b0, l0, q0, l1, sd1, l1, q1, l2, b1, l3, sf0;
cellb2  : mbac, b0, l0, q0, l1, sd2, l1, q1, l2, b1, l3, sf1;
cellc   : mbac, b0, l0, q0, l1, sd3, l1, q1, l2, b1, l3, sf1;
cell1   : mbac, b0, l0, q0a, l1, sd4, l1, q1a, l2, b1a, l3, sf2;
mbacell : -cellb2, cellb1, -cella, cella, -cellb1, cellb2;
mcell   : mb0, ml3, mq1, ml2, sd5, ml2, mq2, ml1, sf2;
mund    : ul0, sh1, ul0, uq1, ul1, uq2, ul2, sh2, ul3, uq3, ul4;
match   : cell1, -mcell, mund;
arc     : -match, cellc, mbacell, -cellc, match, nper=16;

{..igene dateien\opa\examples\example-cells\16periods\7bend_16periods.opa}
