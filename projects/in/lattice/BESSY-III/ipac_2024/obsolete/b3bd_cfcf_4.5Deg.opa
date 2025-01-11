{..lat\20231124_4.5deg_allesaufnull\b3bd_cfcf_bb1.25m_a0a1_1.56_6.25_.opa}


energy = 2.500000;
rotinv = 0;
    betax   = 2.1018363; alphax  = 0.0000000;
    etax    = 0.0000000; etaxp   = 0.0000000;
    betay   = 3.3802562; alphay  = 0.0000000;
    etay    = 0.0000000; etayp   = 0.0000000;

{----- variables ---------------------------------------------------}

brv   = -0.44;
bbv   = 4.5-2.0*brv;
mbrv  = -0.126755112936131;
mbbv  = 2.25-mbrv;

{----- table of elements ---------------------------------------------}

l1     : drift, l = 0.100000, ax = 9.00, ay = 9.00;
ml1    : drift, l = 0.100000, ax = 9.00, ay = 9.00;
ml1a   : drift, l = 0.200000, ax = 9.00, ay = 9.00;
ml1b   : drift, l = 0.100000, ax = 9.00, ay = 9.00;
ul1    : drift, l = 0.354784, ax = 9.00, ay = 9.00;
ul2    : drift, l = 0.120000, ax = 9.00, ay = 9.00;
ul3    : drift, l = 0.100000, ax = 9.00, ay = 9.00;
ul4    : drift, l = 2.800000, ax = 9.00, ay = 9.00;

uq1    : quadrupole, l = 0.100000, k = -7.429366, ax = 9.00, ay = 9.00;
uq2    : quadrupole, l = 0.250000, k = 9.061852, ax = 9.00, ay = 9.00;
uq3    : quadrupole, l = 0.100000, k = -9.984271, ax = 9.00, ay = 9.00;

sd     : sextupole, l = 0.070000, k = -231.652249, n = 1, ax = 9.00,
         ay = 9.00;
sf     : sextupole, l = 0.110000, k = 242.166329, n = 1, ax = 9.00,
         ay = 9.00;

om_sd  : opticsmarker, betax = 3.347847, alphax = 3.730128, betay = 3.756874,
         alphay = -3.102013, etax  = 0.040228, etaxp  = -0.048817,
         etay  = 0.000000, etayp  = 0.000000, ax = 9.00, ay = 9.00;
om_sf  : opticsmarker, betax = 5.242950, alphax = 0.000000, betay = 2.478150,
         alphay = 0.000000, etax  = 0.051235, etaxp  = 0.000000,
         etay  = 0.000000, etayp  = 0.000000, ax = 9.00, ay = 9.00;
om_mbb : opticsmarker, betax = 1.805766, alphax = -2.688999, betay = 8.675730,
         alphay = 2.250039, etax  = 0.000000, etaxp  = 0.000000,
         etay  = 0.000000, etayp  = 0.000000, ax = 9.00, ay = 9.00;
om_c   : opticsmarker, betax = 5.261120, alphax = 0.000000, betay = 2.498280,
         alphay = 0.000000, etax  = 0.051482, etaxp  = 0.000000,
         etay  = 0.000000, etayp  = 0.000000, ax = 9.00, ay = 9.00;

bb     : combined, l = 1.250000, t = bbv, k = -1.067836, t1 = bbv/2.0,
         t2 = bbv/2.0, ax = 9.00, ay = 9.00;
br     : combined, l = 0.140000, t = brv, k = 7.151128, t1 = brv/2.0,
         t2 = brv/2.0, ax = 9.00, ay = 9.00;
mbb    : combined, l = 1.200000, t = mbbv, k = -0.700000, t1 = mbbv/2.0,
         t2 = mbbv/2.0, ax = 9.00, ay = 9.00;
mbr    : combined, l = 0.110000, t = mbrv, k = 8.418136, t1 = mbrv/2.0,
         t2 = mbrv/2.0, ax = 9.00, ay = 9.00;


{----- table of segments ---------------------------------------------}

ucell : om_sf, sf, l1, br, l1, sd, om_sd, sd, l1, bb, l1, sd, om_sd, sd, l1,
        br, l1, sf, om_sf;
dcell : om_sf, sf, ml1, mbr, ml1b, sd, om_sd, sd, ml1a, mbb, om_mbb;
arc   : -dcell, ucell, ucell, om_c, ucell, ucell, dcell;
mund  : om_mbb, ul1, uq1, ul2, uq2, ul3, uq3, ul4;
sec   : -mund, arc, mund;
sech  : om_c, ucell, ucell, dcell, mund;
sec16 : -mund, arc, mund, nper=16;
ring  : 16*sec;

{..lat\20231124_4.5deg_allesaufnull\b3bd_cfcf_bb1.25m_a0a1_1.56_6.25_.opa}
