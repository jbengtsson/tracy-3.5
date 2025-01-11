{..ices\b3_lat\20231124_4.5deg_allesaufnull\b3bd_sfsf4q_bb1.20m_a0a1_.opa}


energy = 2.500000;
rotinv = 0;
    betax   = 2.3284313; alphax  = 0.0000000;
    etax    = 0.0000000; etaxp   = 0.0000000;
    betay   = 4.1234308; alphay  = 0.0000000;
    etay    = 0.0000000; etayp   = 0.0000000;

{----- variables ---------------------------------------------------}

brv   = -0.34;
bbv   = 4.5-2.0*brv;
mbrv  = -0.0842697566666138;
mbbv  = 2.25-mbrv;

{----- table of elements ---------------------------------------------}

l1     : drift, l = 0.100000, ax = 9.00, ay = 9.00;
ml1    : drift, l = 0.100000, ax = 9.00, ay = 9.00;
ml1a   : drift, l = 0.250000, ax = 9.00, ay = 9.00;
ml1b   : drift, l = 0.100000, ax = 9.00, ay = 9.00;
ml1c   : drift, l = 0.100000, ax = 9.00, ay = 9.00;
ul1    : drift, l = 0.200000, ax = 9.00, ay = 9.00;
ul2    : drift, l = 0.430000, ax = 9.00, ay = 9.00;
ul3    : drift, l = 0.150000, ax = 9.00, ay = 9.00;
ul4    : drift, l = 0.100000, ax = 9.00, ay = 9.00;
ul5    : drift, l = 2.800000, ax = 9.00, ay = 9.00;

qd     : quadrupole, l = 0.125000, k = -9.361708, ax = 9.00, ay = 9.00;
mqd    : quadrupole, l = 0.130000, k = -9.200906, ax = 9.00, ay = 9.00;
uq1    : quadrupole, l = 0.100000, k = 4.705477, ax = 9.00, ay = 9.00;
uq2    : quadrupole, l = 0.150000, k = -9.491984, ax = 9.00, ay = 9.00;
uq3    : quadrupole, l = 0.240000, k = 9.148943, ax = 9.00, ay = 9.00;
uq4    : quadrupole, l = 0.100000, k = -8.789846, ax = 9.00, ay = 9.00;

bb     : bending, l = 1.250000, t = bbv, k = 0.000000, t1 = bbv/2.0,
         t2 = bbv/2.0, ax = 9.00, ay = 9.00;
mbb    : bending, l = 1.000000, t = mbbv, k = 0.000000, t1 = mbbv/2.0,
         t2 = mbbv/2.0, ax = 9.00, ay = 9.00;

sd     : sextupole, l = 0.050000, k = -187.984660, n = 1, ax = 9.00,
         ay = 9.00;
sf     : sextupole, l = 0.050000, k = 218.556304, n = 1, ax = 9.00,
         ay = 9.00;

om_sd  : opticsmarker, betax = 1.691706, alphax = 1.762554, betay = 5.746677,
         alphay = 0.229201, etax  = 0.032591, etaxp  = -0.040746,
         etay  = 0.000000, etayp  = 0.000000, ax = 9.00, ay = 9.00;
om_sf  : opticsmarker, betax = 6.106480, alphax = 0.000000, betay = 2.990900,
         alphay = 0.000000, etax  = 0.064748, etaxp  = 0.000000,
         etay  = 0.000000, etayp  = 0.000000, ax = 9.00, ay = 9.00;
om_mbb : opticsmarker, betax = 1.210750, alphax = -1.392509, betay = 5.451241,
         alphay = 0.000000, etax  = 0.000000, etaxp  = 0.000000,
         etay  = 0.000000, etayp  = 0.000000, ax = 9.00, ay = 9.00;
om_c   : opticsmarker, betax = 5.972870, alphax = 0.000000, betay = 2.949620,
         alphay = 0.000000, etax  = 0.063420, etaxp  = 0.000000,
         etay  = 0.000000, etayp  = 0.000000, ax = 9.00, ay = 9.00;

br     : combined, l = 0.170000, t = brv, k = 8.419588, t1 = brv/2.0,
         t2 = brv/2.0, ax = 9.00, ay = 9.00;
mbr    : combined, l = 0.180000, t = mbrv, k = 8.074015, t1 = mbrv/2.0,
         t2 = mbrv/2.0, ax = 9.00, ay = 9.00;


{----- table of segments ---------------------------------------------}

ucell : om_sf, sf, l1, br, l1, qd, l1, sd, om_sd, sd, l1, bb, l1, sd, om_sd,
        sd, l1, qd, l1, br, l1, sf, om_sf;
dcell : om_sf, sf, ml1, mbr, ml1c, mqd, ml1b, sd, om_sd, sd, ml1a, mbb,
        om_mbb;
arc   : -dcell, ucell, ucell, om_c, ucell, ucell, dcell;
mund  : om_mbb, ul1, uq1, ul2, uq2, ul3, uq3, ul4, uq4, ul5;
sec   : -mund, arc, mund;
sech  : om_c, ucell, ucell, dcell, mund;
sec16 : -mund, arc, mund, nper=16;
ring  : 16*sec;

{..ices\b3_lat\20231124_4.5deg_allesaufnull\b3bd_sfsf4q_bb1.20m_a0a1_.opa}
