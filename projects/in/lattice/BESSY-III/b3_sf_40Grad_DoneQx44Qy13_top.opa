{..xtcloud\biii_meins\lattices\b3_lat\b3_sf_40grad_doneqx44qy13_top_z.opa}


energy = 2.500000;

    betax   = 2.6555738; alphax  = 0.0000000;
    etax    = 0.0000654; etaxp   = 0.0000000;
    betay   = 2.5202533; alphay  = 0.0000000;
    etay    = 0.0000000; etayp   = 0.0000000;

{----- variables ---------------------------------------------------}

b2a   = -0.23;
b1a   = 4.0-2.0*b2a;
mb2a  = -0.244832764868045;
mb1a  = 3.25-mb2a;

{----- table of elements ---------------------------------------------}

l1     : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ml1    : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ml1a   : drift, l = 0.140000, ax = 300.00, ay = 300.00;
ml1b   : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ul1    : drift, l = 0.280000, ax = 300.00, ay = 300.00;
ul2    : drift, l = 0.253150, ax = 300.00, ay = 300.00;
ul3    : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ul4    : drift, l = 2.800000, ax = 300.00, ay = 300.00;

q1     : quadrupole, l = 0.130000, k = -9.366917, ax = 9.00, ay = 9.00;
mq1    : quadrupole, l = 0.130000, k = -5.000000, ax = 9.00, ay = 9.00;
uq1    : quadrupole, l = 0.100000, k = -4.009569, ax = 9.00, ay = 9.00;
uq2    : quadrupole, l = 0.240000, k = 9.290481, ax = 9.00, ay = 9.00;
uq3    : quadrupole, l = 0.120000, k = -9.513104, ax = 9.00, ay = 9.00;

b1     : bending, l = 1.000000, t = b1a, k = 0.000000, t1 = b1a/2.0,
         t2 = b1a/2.0, ax = 50.00, ay = 300.00;
b2     : bending, l = 0.160000, t = b2a, k = 9.389515, t1 = b2a/2.0,
         t2 = b2a/2.0, ax = 50.00, ay = 300.00;
mb1    : bending, l = 0.700000, t = mb1a, k = -1.600000, t1 = mb1a/2.0,
         t2 = mb1a/2.0, ax = 50.00, ay = 300.00;
mb2    : bending, l = 0.160000, t = mb2a, k = 7.896711, t1 = mb2a/2.0,
         t2 = mb2a/2.0, ax = 50.00, ay = 300.00;

s1     : sextupole, l = 0.120000, k = -233.935253, n = 1, ax = 300.00,
         ay = 300.00;
s2     : sextupole, l = 0.070000, k = 246.209678, n = 1, ax = 300.00,
         ay = 300.00;

om_s2  : opticsmarker, betax = 5.550480, alphax = 0.000000, betay = 2.810490,
         alphay = 0.000000, etax  = 0.053363, etaxp  = 0.000000,
         etay  = 0.000000, etayp  = 0.000000, ax = 50.00, ay = 50.00;
om_mb1 : opticsmarker, betax = 0.881052, alphax = -2.423020, betay = 9.665924,
         alphay = 3.622847, etax  = 0.000000, etaxp  = 0.000000,
         etay  = 0.000000, etayp  = 0.000000, ax = 50.00, ay = 50.00;
om_c   : opticsmarker, betax = 5.681470, alphax = 0.000000, betay = 2.779680,
         alphay = 0.000000, etax  = 0.060370, etaxp  = 0.000000,
         etay  = 0.000000, etayp  = 0.000000, ax = 50.00, ay = 50.00;


{----- table of segments ---------------------------------------------}

ucell : om_s2, s2, l1, b2, l1, q1, l1, s1, l1, b1, l1, s1, l1, q1, l1, b2, l1,
        s2, om_s2;
dcell : om_s2, s2, ml1, mb2, ml1b, s1, ml1a, mb1, om_mb1;
arc   : -dcell, 2*ucell, om_c, 2*ucell, dcell;
mund  : ul1, uq1, ul2, uq2, ul3, uq3, ul4;
sec   : -mund, arc, mund;
sech  : om_c, 2*ucell, dcell, mund;
sec16 : -mund, arc, mund, nper=16;
ring  : 16*sec;

{..xtcloud\biii_meins\lattices\b3_lat\b3_sf_40grad_doneqx44qy13_top_z.opa}
