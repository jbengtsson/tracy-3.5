{d:\lattices\b3_lat\b3_sf_uc4375_dsc_sf_4372_1276.opa}


energy = 2.500000;

    betax   = 3.0306951; alphax  = 0.0000000;
    etax    = 0.0000000; etaxp   = 0.0000000;
    betay   = 2.9470044; alphay  = 0.0000000;
    etay    = 0.0000000; etayp   = 0.0000000;

{----- variables ---------------------------------------------------}

brv   = -0.34;
bbv   = 4.375-2.0*brv;
mbrv  = -0.0852053644577534;
mbbv  = 2.5-mbrv;

{----- table of elements ---------------------------------------------}

l1     : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ml1    : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ml1a   : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ml1b   : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ml1c   : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ul1    : drift, l = 0.120000, ax = 300.00, ay = 300.00;
ul2    : drift, l = 0.250000, ax = 300.00, ay = 300.00;
ul3    : drift, l = 0.100000, ax = 300.00, ay = 300.00;
ul4    : drift, l = 2.800000, ax = 300.00, ay = 300.00;

qd     : quadrupole, l = 0.125000, k = -9.573181, ax = 9.00, ay = 9.00;
mqd    : quadrupole, l = 0.130000, k = -9.753806, ax = 9.00, ay = 9.00;
uq1    : quadrupole, l = 0.150000, k = 9.025159, ax = 9.00, ay = 9.00;
uq2    : quadrupole, l = 0.250000, k = -9.108800, ax = 9.00, ay = 9.00;
uq3    : quadrupole, l = 0.180000, k = 9.434980, ax = 9.00, ay = 9.00;

bb     : bending, l = 1.100000, t = bbv, k = 0.000000, t1 = bbv/2.0,
         t2 = bbv/2.0, ax = 50.00, ay = 300.00;
mbb    : bending, l = 0.850000, t = mbbv, k = 0.000000, t1 = mbbv/2.0,
         t2 = mbbv/2.0, ax = 50.00, ay = 300.00;

sd     : sextupole, l = 0.110000, k = -232.753589, n = 1, ax = 300.00,
         ay = 300.00;
sf     : sextupole, l = 0.050000, k = 240.259413, n = 1, ax = 300.00,
         ay = 300.00;

om_sf  : opticsmarker, betax = 5.797560, alphax = 0.000000, betay = 2.844570,
         alphay = 0.000000, etax  = 0.059669, etaxp  = 0.000000,
         etay  = 0.000000, etayp  = 0.000000, ax = 50.00, ay = 50.00;
om_mbb : opticsmarker, betax = 0.762357, alphax = -1.177548, betay = 5.682860,
         alphay = 0.000000, etax  = 0.000000, etaxp  = 0.000000,
         etay  = 0.000000, etayp  = 0.000000, ax = 50.00, ay = 50.00;
om_c   : opticsmarker, betax = 6.106180, alphax = 0.000000, betay = 2.990910,
         alphay = 0.000000, etax  = 0.064601, etaxp  = 0.000000,
         etay  = 0.000000, etayp  = 0.000000, ax = 50.00, ay = 50.00;

br     : combined, l = 0.170000, t = brv, k = 8.660561, t1 = brv/2.0,
         t2 = brv/2.0, ax = 50.00, ay = 300.00;
mbr    : combined, l = 0.180000, t = mbrv, k = 8.862581, t1 = mbrv/2.0,
         t2 = mbrv/2.0, ax = 50.00, ay = 300.00;


{----- table of segments ---------------------------------------------}

ucell : om_sf, sf, l1, br, l1, qd, l1, sd, l1, bb, l1, sd, l1, qd, l1, br, l1,
        sf, om_sf;
dcell : om_sf, sf, ml1, mbr, ml1c, mqd, ml1b, sd, ml1a, mbb, om_mbb;
arc   : -dcell, 2*ucell, om_c, 2*ucell, dcell;
mund  : om_mbb, ul1, uq1, ul2, uq2, ul3, uq3, ul4;
sec   : -mund, arc, mund;
sech  : om_c, 2*ucell, dcell, mund;
sec16 : -mund, arc, mund, nper=16;
ring  : 16*sec;

{d:\lattices\b3_lat\b3_sf_uc4375_dsc_sf_4372_1276.opa}
