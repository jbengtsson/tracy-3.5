{..0231124_4.5deg_allesaufnull\b3bd_sfsf4q_bb1bb2_vary_morealpha_top1.opa}


energy = 2.500000;
rotinv = 0;
    betax   = 2.3742338; alphax  = 0.0000000;
    etax    = 0.0000000; etaxp   = 0.0000000;
    betay   = 4.2860024; alphay  = 0.0000000;
    etay    = 0.0000000; etayp   = 0.0000000;

{----- variables ---------------------------------------------------}

brb    = -0.3;
bb1b   = 0.45;
bb2b   = (4.5-2*bb1b)-2.0*brb;
mbrb   = -0.31015962150778;
mbb1b  = 0.45;
mbb2b  = 2.25-mbb1b-mbrb;

{----- table of elements ---------------------------------------------}

l1     : drift, l = 0.100000, ax = 9.00, ay = 9.00;
ml1    : drift, l = 0.100000, ax = 9.00, ay = 9.00;
ml1a   : drift, l = 0.100000, ax = 9.00, ay = 9.00;
ml1b   : drift, l = 0.100000, ax = 9.00, ay = 9.00;
ml1c   : drift, l = 0.100000, ax = 9.00, ay = 9.00;
ul1    : drift, l = 0.250000, ax = 9.00, ay = 9.00;
ul2    : drift, l = 0.250000, ax = 9.00, ay = 9.00;
ul3    : drift, l = 0.200000, ax = 9.00, ay = 9.00;
ul4    : drift, l = 0.100000, ax = 9.00, ay = 9.00;
ul5    : drift, l = 2.800000, ax = 9.00, ay = 9.00;

qd     : quadrupole, l = 0.125000, k = -9.361708, ax = 9.00, ay = 9.00;
mqd    : quadrupole, l = 0.130000, k = -9.200906, ax = 9.00, ay = 9.00;
uq1    : quadrupole, l = 0.100000, k = 6.792969, ax = 9.00, ay = 9.00;
uq2    : quadrupole, l = 0.180000, k = -8.180370, ax = 9.00, ay = 9.00;
uq3    : quadrupole, l = 0.240000, k = 8.732444, ax = 9.00, ay = 9.00;
uq4    : quadrupole, l = 0.100000, k = -9.282014, ax = 9.00, ay = 9.00;

bb2    : bending, l = 0.800000, t = bb2b, k = 0.000000, t1 = bb2b/2.0,
         t2 = bb2b/2.0, ax = 9.00, ay = 9.00;
mbb2   : bending, l = 1.000000, t = mbb2b, k = 0.000000, t1 = mbb2b/2.0,
         t2 = mbb2b/2.0, ax = 9.00, ay = 9.00;

sd     : sextupole, l = 0.050000, k = -249.287977, n = 1, ax = 9.00,
         ay = 9.00;
sf     : sextupole, l = 0.050000, k = 265.950890, n = 1, ax = 9.00,
         ay = 9.00;

om_sd  : opticsmarker, betax = 1.154072, alphax = 1.455183, betay = 5.302664,
         alphay = 0.112571, etax  = 0.023938, etaxp  = -0.036833,
         etay  = 0.000000, etayp  = 0.000000, ax = 9.00, ay = 9.00;
om_sf  : opticsmarker, betax = 5.293770, alphax = 0.000000, betay = 2.557610,
         alphay = 0.000000, etax  = 0.056002, etaxp  = 0.000000,
         etay  = 0.000000, etayp  = 0.000000, ax = 9.00, ay = 9.00;
om_mbb : opticsmarker, betax = 1.378946, alphax = -1.650765, betay = 5.289277,
         alphay = -0.100001, etax  = 0.000000, etaxp  = 0.000000,
         etay  = 0.000000, etayp  = 0.000000, ax = 9.00, ay = 9.00;
om_c   : opticsmarker, betax = 5.293840, alphax = 0.000000, betay = 2.557700,
         alphay = 0.000000, etax  = 0.055901, etaxp  = 0.000000,
         etay  = 0.000000, etayp  = 0.000000, ax = 9.00, ay = 9.00;

bb1    : combined, l = 0.160000, t = bb1b, k = -7.485289, t1 = bb1b/2.0,
         t2 = bb1b/2.0, ax = 9.00, ay = 9.00;
mbb1   : combined, l = 0.170000, t = mbb1b, k = -7.017857, t1 = mbb1b/2.0,
         t2 = mbb1b/2.0, ax = 9.00, ay = 9.00;
br     : combined, l = 0.180000, t = brb, k = 8.348017, t1 = brb/2.0,
         t2 = brb/2.0, ax = 9.00, ay = 9.00;
mbr    : combined, l = 0.180000, t = mbrb, k = 8.296055, t1 = mbrb/2.0,
         t2 = mbrb/2.0, ax = 9.00, ay = 9.00;


{----- table of segments ---------------------------------------------}

ucell : om_sf, sf, l1, br, l1, bb1, l1, sd, om_sd, sd, l1, bb2, l1, sd, om_sd,
        sd, l1, bb1, l1, br, l1, sf, om_sf;
dcell : om_sf, sf, ml1, mbr, ml1c, mbb1, ml1b, sd, om_sd, sd, ml1a, mbb2,
        om_mbb;
arc   : -dcell, ucell, ucell, om_c, ucell, ucell, dcell;
mund  : om_mbb, ul1, uq1, ul2, uq2, ul3, uq3, ul4, uq4, ul5;
sec   : -mund, arc, mund;
sech  : om_c, ucell, ucell, dcell, mund;
sec16 : -mund, arc, mund, nper=16;
ring  : 16*sec;

{..0231124_4.5deg_allesaufnull\b3bd_sfsf4q_bb1bb2_vary_morealpha_top1.opa}
