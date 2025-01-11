{..tcloud\biii_meins\lattices\b3_lat\b3_cf_425grad_doneqx44qy13_top_z.opa}


energy = 2.500000;

    betax   = 2.7870479; alphax  = 0.0000000;
    etax    = 0.0000003; etaxp   = 0.0000000;
    betay   = 2.5614845; alphay  = 0.0000000;
    etay    = 0.0000000; etayp   = 0.0000000;

{----- variables ---------------------------------------------------}

b2a   = -0.38;
b1a   = 4.25-2.0*b2a;
mb2a  = -0.160989817118871;
mb1a  = 2.75-mb2a;

{----- table of elements ---------------------------------------------}

l1     : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ml1    : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ml1a   : drift, l = 0.260000, ax = 300.00, ay = 300.00;
ml1b   : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ul1    : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ul2    : drift, l = 0.258148, ax = 300.00, ay = 300.00;
ul3    : drift, l = 0.300000, ax = 300.00, ay = 300.00;
ul4    : drift, l = 2.800000, ax = 300.00, ay = 300.00;

uq1    : quadrupole, l = 0.165000, k = -9.476508, ax = 9.00, ay = 9.00;
uq2    : quadrupole, l = 0.200000, k = 9.591353, ax = 9.00, ay = 9.00;
uq3    : quadrupole, l = 0.100000, k = -3.051951, ax = 9.00, ay = 9.00;

b1     : bending, l = 1.100000, t = b1a, k = -1.199036, t1 = b1a/2.0,
         t2 = b1a/2.0, ax = 50.00, ay = 300.00;
b2     : bending, l = 0.163000, t = b2a, k = 6.190683, t1 = b2a/2.0,
         t2 = b2a/2.0, ax = 50.00, ay = 300.00;
mb1    : bending, l = 0.600000, t = mb1a, k = 0.000000, t1 = mb1a/2.0,
         t2 = mb1a/2.0, ax = 50.00, ay = 300.00;
mb2    : bending, l = 0.163000, t = mb2a, k = 6.476866, t1 = mb2a/2.0,
         t2 = mb2a/2.0, ax = 50.00, ay = 300.00;

s1     : sextupole, l = 0.180000, k = -223.646713, n = 1, ax = 300.00,
         ay = 300.00;
s2     : sextupole, l = 0.130000, k = 224.408612, n = 1, ax = 300.00,
         ay = 300.00;

om_s2  : opticsmarker, betax = 4.981030, alphax = 0.000000, betay = 2.530010,
         alphay = 0.000000, etax  = 0.046882, etaxp  = 0.000000,
         etay  = 0.000000, etayp  = 0.000000, ax = 50.00, ay = 50.00;
om_mb1 : opticsmarker, betax = 0.660220, alphax = -1.587761, betay = 15.524128,
         alphay = -7.452246, etax  = 0.000000, etaxp  = 0.000000,
         etay  = 0.000000, etayp  = 0.000000, ax = 50.00, ay = 50.00;
om_c   : opticsmarker, betax = 4.975580, alphax = 0.000000, betay = 2.608240,
         alphay = 0.000000, etax  = 0.047311, etaxp  = 0.000000,
         etay  = 0.000000, etayp  = 0.000000, ax = 50.00, ay = 50.00;


{----- table of segments ---------------------------------------------}

ucell : om_s2, s2, l1, b2, l1, s1, l1, b1, l1, s1, l1, b2, l1, s2, om_s2;
dcell : om_s2, s2, ml1, mb2, ml1b, s1, ml1a, mb1, om_mb1;
arc   : -dcell, 2*ucell, om_c, 2*ucell, dcell;
mund  : om_mb1, ul1, uq1, ul2, uq2, ul3, uq3, ul4;
sec   : -mund, arc, mund;
sech  : om_c, 2*ucell, dcell, mund;
sec16 : -mund, arc, mund, nper=16;
ring  : 16*sec;

{..tcloud\biii_meins\lattices\b3_lat\b3_cf_425grad_doneqx44qy13_top_z.opa}
