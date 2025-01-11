{..es\b3_lat\b3_cf_2poles_45grad_v3_tunesinteger_shortquads_longbends.opa}


energy = 2.500000;

    betax   = 2.3160326; alphax  = 0.0000000;
    etax    = 0.0000000; etaxp   = 0.0000000;
    betay   = 2.1449538; alphay  = 0.0000000;
    etay    = 0.0000000; etayp   = 0.0000000;

{----- variables ---------------------------------------------------}

b2a   = -0.4;
b1a   = 4.5-2.0*b2a;
mb2a  = 0;
mb1a  = 2.25-mb2a;

{----- table of elements ---------------------------------------------}

l1     : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ml1    : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ml1a   : drift, l = 0.226126, ax = 300.00, ay = 300.00;
ml1b   : drift, l = 0.568183, ax = 300.00, ay = 300.00;
ul1    : drift, l = 0.280150, ax = 300.00, ay = 300.00;
ul2    : drift, l = 0.186225, ax = 300.00, ay = 300.00;
ul3    : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ul4    : drift, l = 2.800000, ax = 300.00, ay = 300.00;

uq1    : quadrupole, l = 0.150000, k = -7.785685, ax = 9.00, ay = 9.00;
uq2    : quadrupole, l = 0.200000, k = 8.451754, ax = 9.00, ay = 9.00;
uq3    : quadrupole, l = 0.100000, k = -5.985905, ax = 9.00, ay = 9.00;

b1     : bending, l = 1.000000, t = b1a, k = -1.571294, t1 = b1a/2.0,
         t2 = b1a/2.0, ax = 50.00, ay = 300.00;
b2     : bending, l = 0.140000, t = b2a, k = 8.479348, t1 = b2a/2.0,
         t2 = b2a/2.0, ax = 50.00, ay = 300.00;
mb1    : bending, l = 0.280000, t = mb1a, k = 0.000000, t1 = mb1a/2.0,
         t2 = mb1a/2.0, ax = 50.00, ay = 300.00;
mb2    : bending, l = 0.120000, t = mb2a, k = 7.755013, t1 = mb2a/2.0,
         t2 = mb2a/2.0, ax = 50.00, ay = 300.00;

s1     : sextupole, l = 0.100000, k = -461.665510, n = 1, ax = 300.00,
         ay = 300.00;
s2     : sextupole, l = 0.050000, k = 653.439253, n = 1, ax = 300.00,
         ay = 300.00;
s3     : sextupole, l = 0.100000, k = 0.000000, n = 1, ax = 300.00,
         ay = 300.00;
s4     : sextupole, l = 0.100000, k = 0.000000, n = 1, ax = 300.00,
         ay = 300.00;

om_s2  : opticsmarker, betax = 4.446660, alphax = 0.000000, betay = 2.027240,
         alphay = 0.000000, etax  = 0.043000, etaxp  = 0.000000,
         etay  = 0.000000, etayp  = 0.000000, ax = 50.00, ay = 50.00;
om_mb1 : opticsmarker, betax = 0.379946, alphax = -0.897109, betay = 10.385652,
         alphay = -6.173339, etax  = 0.000000, etaxp  = 0.000000,
         etay  = 0.000000, etayp  = 0.000000, ax = 50.00, ay = 50.00;
om_c   : opticsmarker, betax = 3.280220, alphax = 0.000000, betay = 1.778260,
         alphay = 0.000000, etax  = 0.033189, etaxp  = 0.000000,
         etay  = 0.000000, etayp  = 0.000000, ax = 50.00, ay = 50.00;


{----- table of segments ---------------------------------------------}

ucell : om_s2, s2, l1, b2, l1, s1, l1, b1, l1, s1, l1, b2, l1, s2, om_s2;
dcell : om_s2, s2, ml1, mb2, ml1b, s1, ml1a, mb1, om_mb1;
arc   : -dcell, 2*ucell, om_c, 2*ucell, dcell;
mund  : ul1, uq1, ul2, s3, ul2, uq2, ul3, s4, ul3, uq3, ul4;
sec   : -mund, arc, mund;
secm  : dcell, mund;
sech  : om_c, 2*ucell, dcell, mund;
sec16 : -mund, arc, mund, nper=16;
ring  : 16*sec;

{..es\b3_lat\b3_cf_2poles_45grad_v3_tunesinteger_shortquads_longbends.opa}
