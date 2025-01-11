{..3_lat\b3_cf_4poles_45grad_v3_tunesinteger_shortquads_increasealpha.opa}


energy = 2.500000;

    betax   = 2.5053043; alphax  = 0.0000000;
    etax    = 0.0000367; etaxp   = 0.0000000;
    betay   = 2.2809490; alphay  = 0.0000000;
    etay    = 0.0000000; etayp   = 0.0000000;

{----- variables ---------------------------------------------------}

b2a   = -0,28;
b1a   = 4.5-2.0*b2a;
mb2a  = 0;
mb1a  = 2.25-mb2a;

{----- table of elements ---------------------------------------------}

l1     : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ml1    : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ml1a   : drift, l = 0.218000, ax = 300.00, ay = 300.00;
ml1b   : drift, l = 0.300000, ax = 300.00, ay = 300.00;
ul1    : drift, l = 0.104038, ax = 300.00, ay = 300.00;
ul2    : drift, l = 0.300000, ax = 300.00, ay = 300.00;
ul3    : drift, l = 0.200000, ax = 300.00, ay = 300.00;
ul4    : drift, l = 2.800000, ax = 300.00, ay = 300.00;

uq1    : quadrupole, l = 0.140000, k = -8.133350, ax = 9.00, ay = 9.00;
uq2    : quadrupole, l = 0.180000, k = 8.158217, ax = 9.00, ay = 9.00;
uq3    : quadrupole, l = 0.080000, k = -7.270339, ax = 9.00, ay = 9.00;

b1     : bending, l = 0.600000, t = b1a, k = -3.085460, t1 = b1a/2.0,
         t2 = b1a/2.0, ax = 50.00, ay = 300.00;
b2     : bending, l = 0.175000, t = b2a, k = 8.378757, t1 = b2a/2.0,
         t2 = b2a/2.0, ax = 50.00, ay = 300.00;
mb1    : bending, l = 0.300000, t = mb1a, k = 0.000000, t1 = mb1a/2.0,
         t2 = mb1a/2.0, ax = 50.00, ay = 300.00;
mb2    : bending, l = 0.175000, t = mb2a, k = 6.910000, t1 = mb2a/2.0,
         t2 = mb2a/2.0, ax = 50.00, ay = 300.00;

s1     : sextupole, l = 0.100000, k = 0.000000, n = 1, ax = 300.00,
         ay = 300.00;
s2     : sextupole, l = 0.050000, k = 0.000000, n = 1, ax = 300.00,
         ay = 300.00;
s3     : sextupole, l = 0.100000, k = 0.000000, n = 1, ax = 300.00,
         ay = 300.00;
s4     : sextupole, l = 0.100000, k = 0.000000, n = 1, ax = 300.00,
         ay = 300.00;

om_s2  : opticsmarker, betax = 3.280880, alphax = 0.000000, betay = 1.778220,
         alphay = 0.000000, etax  = 0.033618, etaxp  = 0.000000,
         etay  = 0.000000, etayp  = 0.000000, ax = 50.00, ay = 50.00;
om_mb1 : opticsmarker, betax = 0.380408, alphax = -0.887775, betay = 10.372159,
         alphay = -6.138640, etax  = 0.000028, etaxp  = 0.000001,
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

{..3_lat\b3_cf_4poles_45grad_v3_tunesinteger_shortquads_increasealpha.opa}
