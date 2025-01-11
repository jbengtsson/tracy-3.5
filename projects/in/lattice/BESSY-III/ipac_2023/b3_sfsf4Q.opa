{d:\lattices\b3_lat\b3_sfsf_4quads_xsextadjust.opa}


energy = 2.500000;

    betax   = 2.6692567; alphax  = 0.0000000;
    etax    = 0.0000000; etaxp   = 0.0000000;
    betay   = 3.5165057; alphay  = 0.0000000;
    etay    = 0.0000000; etayp   = 0.0000000;

{----- variables ---------------------------------------------------}

brv   = -0.34;
bbv   = 4.375-2.0*brv;
mbrv  = -0.0889026597923138;
mbbv  = 2.5-mbrv;

{----- table of elements ---------------------------------------------}

l1     : drift, l = 0.100000, ax = 9.00, ay = 9.00;
ml1    : drift, l = 0.100000, ax = 9.00, ay = 9.00;
ml1a   : drift, l = 0.100000, ax = 9.00, ay = 9.00;
ml1b   : drift, l = 0.100000, ax = 9.00, ay = 9.00;
ml1c   : drift, l = 0.100000, ax = 9.00, ay = 9.00;
ul1    : drift, l = 0.150000, ax = 9.00, ay = 9.00;
ul2    : drift, l = 0.300000, ax = 9.00, ay = 9.00;
ul3    : drift, l = 0.150000, ax = 9.00, ay = 9.00;
ul4    : drift, l = 0.100000, ax = 9.00, ay = 9.00;
ul5    : drift, l = 2.800000, ax = 9.00, ay = 9.00;

qd     : quadrupole, l = 0.125000, k = -9.602787, ax = 9.00, ay = 9.00;
mqd    : quadrupole, l = 0.130000, k = -9.782278, ax = 9.00, ay = 9.00;
uq1    : quadrupole, l = 0.100000, k = 6.144712, ax = 9.00, ay = 9.00;
uq2    : quadrupole, l = 0.170000, k = -9.390129, ax = 9.00, ay = 9.00;
uq3    : quadrupole, l = 0.250000, k = 9.183668, ax = 9.00, ay = 9.00;
uq4    : quadrupole, l = 0.100000, k = -8.997373, ax = 9.00, ay = 9.00;

bb     : bending, l = 1.100000, t = bbv, k = 0.000000, t1 = bbv/2.0,
         t2 = bbv/2.0, ax = 9.00, ay = 9.00;
mbb    : bending, l = 0.850000, t = mbbv, k = 0.000000, t1 = mbbv/2.0,
         t2 = mbbv/2.0, ax = 9.00, ay = 9.00;

sd     : sextupole, l = 0.050000, k = -231.702860, n = 1, ax = 9.00,
         ay = 9.00;
sf     : sextupole, l = 0.055000, k = 231.266556, n = 1, ax = 9.00,
         ay = 9.00;

om_sd  : opticsmarker, betax = 1.602863, alphax = 1.705960, betay = 5.464352,
         alphay = 0.110793, etax  = 0.027576, etaxp  = -0.044142,
         etay  = 0.000000, etayp  = 0.000000, ax = 9.00, ay = 9.00;
om_sf  : opticsmarker, betax = 5.735060, alphax = 0.000000, betay = 2.837550,
         alphay = 0.000000, etax  = 0.059187, etaxp  = 0.000000,
         etay  = 0.000000, etayp  = 0.000000, ax = 9.00, ay = 9.00;
om_mbb : opticsmarker, betax = 0.762357, alphax = -1.177548, betay = 5.682860,
         alphay = 0.000000, etax  = 0.000000, etaxp  = 0.000000,
         etay  = 0.000000, etayp  = 0.000000, ax = 9.00, ay = 9.00;
om_c   : opticsmarker, betax = 6.106180, alphax = 0.000000, betay = 2.990910,
         alphay = 0.000000, etax  = 0.064601, etaxp  = 0.000000,
         etay  = 0.000000, etayp  = 0.000000, ax = 9.00, ay = 9.00;

br     : combined, l = 0.170000, t = brv, k = 8.697105, t1 = brv/2.0,
         t2 = brv/2.0, ax = 9.00, ay = 9.00;
mbr    : combined, l = 0.180000, t = mbrv, k = 8.899091, t1 = mbrv/2.0,
         t2 = mbrv/2.0, ax = 9.00, ay = 9.00;


{----- table of segments ---------------------------------------------}

ucell : om_sf, sf, l1, br, l1, qd, l1, sd, om_sd, sd, l1, bb, l1, sd, om_sd,
        sd, l1, qd, l1, br, l1, sf, om_sf;
dcell : om_sf, sf, ml1, mbr, ml1c, mqd, ml1b, sd, om_sd, sd, ml1a, mbb,
        om_mbb;
arc   : -dcell, 2*ucell, om_c, 2*ucell, dcell;
mund  : om_mbb, ul1, uq1, ul2, uq2, ul3, uq3, ul4, uq4, ul5;
sec   : -mund, arc, mund;
sech  : om_c, 2*ucell, dcell, mund;
sec16 : -mund, arc, mund, nper=16;
ring  : 16*sec;

{d:\lattices\b3_lat\b3_sfsf_4quads_xsextadjust.opa}
