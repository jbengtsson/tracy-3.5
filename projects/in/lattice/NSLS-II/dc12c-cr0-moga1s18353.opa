{c:\users\streun\opadat\sls-2\dc12c.opa}
{com period 12 version of dc01a, only m straights com}

allocation  = sls2names_per12.dat;

energy = 2.400000;

    betax   = 2.6739905; alphax  = 0.0000000;
    etax    = 0.0000226; etaxp   = 0.0000005;
    betay   = 2.7969540; alphay  = 0.0000000;
    etay    = 0.0000000; etayp   = 0.0000000;

{----- variables ----------------------------------------------------}

aban  = 0.78;

{----- table of elements ----------------------------------------------------}

dmon   : drift, l = 0.050000, ax = 10.00, ay = 10.00;
dnvb   : drift, l = 0.010000, ax = 10.00, ay = 10.00;
dsvb   : drift, l = 0.084000, ax = 10.00, ay = 10.00;
dnm    : drift, l = 0.300000, ax = 10.00, ay = 10.00;
dnm1   : drift, l = 0.050000, ax = 10.00, ay = 10.00;
dnm2   : drift, l = 0.150000, ax = 10.00, ay = 10.00;
dxs    : drift, l = 0.200000, ax = 10.00, ay = 10.00;
dxm    : drift, l = 0.218000, ax = 10.00, ay = 10.00;
dmp    : drift, l = 0.150000, ax = 10.00, ay = 10.00;
dsx    : drift, l = 0.100000, ax = 10.00, ay = 10.00;
dme    : drift, l = 0.300000, ax = 10.00, ay = 10.00;
ds1    : drift, l = 0.050000, ax = 10.00, ay = 10.00;
ds2    : drift, l = 0.070000, ax = 10.00, ay = 10.00;
mgp    : drift, l = 1.000000, ax = 10.00, ay = 10.00;
dm1    : drift, l = 0.100000, ax = 10.00, ay = 10.00;
dm2    : drift, l = 0.100000, ax = 10.00, ay = 10.00;
dl1    : drift, l = 0.100000, ax = 10.00, ay = 10.00;
dl2    : drift, l = 0.100000, ax = 10.00, ay = 10.00;
dms    : drift, l = 2.692000, ax = 10.00, ay = 10.00;
dxd    : drift, l = 0.050000, ax = 10.00, ay = 10.00;
doc    : drift, l = 0.050000, ax = 10.00, ay = 10.00;
dbs    : drift, l = 0.074000, ax = 10.00, ay = 10.00;

center : marker, ax = 10.00, ay = 10.00;

qs1    : quadrupole, l = 0.150000, k = -6.442563, ax = 10.00, ay = 10.00;
qs2    : quadrupole, l = 0.200000, k = 8.781815, ax = 10.00, ay = 10.00;
qs3    : quadrupole, l = 0.100000, k = -4.907590, ax = 10.00, ay = 10.00;
qm3      : quadrupole, l =   0.100000, k =  -4.822851, ax = 10.0, ay = 10.0;
qm2      : quadrupole, l =   0.250000, k =   6.216424, ax = 10.0, ay = 10.0;
qm1      : quadrupole, l =   0.150000, k =  -6.552535, ax = 10.0, ay = 10.0;
ql1    : quadrupole, l = 0.150000, k = -6.402326, ax = 10.00, ay = 10.00;
ql2    : quadrupole, l = 0.250000, k = 6.114719, ax = 10.00, ay = 10.00;
ql3    : quadrupole, l = 0.100000, k = -4.455466, ax = 10.00, ay = 10.00;
ql4    : quadrupole, l = 0.100000, k = -6.247519, ax = 10.00, ay = 10.00;
ql5    : quadrupole, l = 0.100000, k = 6.284436, ax = 10.00, ay = 10.00;

bn00   : bending, l = 0.021832, t = 0.312500, k = 0.000000,
         t1 = 0.000000, t2 = 0.312500, ax = 10.00, ay = 10.00;
bn01   : bending, l = 0.022967, t = 0.312559, k = 0.000000,
         t1 = -0.312500, t2 = 0.625059, ax = 10.00, ay = 10.00;
bn02   : bending, l = 0.031095, t = 0.312474, k = 0.000000,
         t1 = -0.625059, t2 = 0.937533, ax = 10.00, ay = 10.00;
bn03   : bending, l = 0.038596, t = 0.312485, k = 0.000000,
         t1 = -0.937533, t2 = 1.250018, ax = 10.00, ay = 10.00;
bn04   : bending, l = 0.045882, t = 0.312491, k = 0.000000,
         t1 = -1.250018, t2 = 1.562509, ax = 10.00, ay = 10.00;
bn05   : bending, l = 0.053107, t = 0.312495, k = 0.000000,
         t1 = -1.562509, t2 = 1.875004, ax = 10.00, ay = 10.00;
bn06   : bending, l = 0.060418, t = 0.312497, k = 0.000000,
         t1 = -1.875004, t2 = 2.187501, ax = 10.00, ay = 10.00;
bs00   : bending, l = 0.004613, t = 0.182296, k = 0.000000,
         t1 = 0.000000, t2 = 0.182296, ax = 10.00, ay = 10.00;
bs01   : bending, l = 0.007205, t = 0.182289, k = 0.000000,
         t1 = -0.182296, t2 = 0.364585, ax = 10.00, ay = 10.00;
bs02   : bending, l = 0.009727, t = 0.182291, k = 0.000000,
         t1 = -0.364585, t2 = 0.546876, ax = 10.00, ay = 10.00;
bs03   : bending, l = 0.012042, t = 0.182292, k = 0.000000,
         t1 = -0.546876, t2 = 0.729168, ax = 10.00, ay = 10.00;
bs04   : bending, l = 0.014194, t = 0.182292, k = 0.000000,
         t1 = -0.729168, t2 = 0.911460, ax = 10.00, ay = 10.00;
bs05   : bending, l = 0.016223, t = 0.182292, k = 0.000000,
         t1 = -0.911460, t2 = 1.093751, ax = 10.00, ay = 10.00;
bs06   : bending, l = 0.018160, t = 0.182292, k = 0.000000,
         t1 = -1.093751, t2 = 1.276043, ax = 10.00, ay = 10.00;
bs07   : bending, l = 0.020024, t = 0.182292, k = 0.000000,
         t1 = -1.276043, t2 = 1.458335, ax = 10.00, ay = 10.00;
bs08   : bending, l = 0.021833, t = 0.182292, k = 0.000000,
         t1 = -1.458335, t2 = 1.640627, ax = 10.00, ay = 10.00;
bs09   : bending, l = 0.023599, t = 0.182292, k = 0.000000,
         t1 = -1.640627, t2 = 1.822918, ax = 10.00, ay = 10.00;
bs10   : bending, l = 0.025288, t = 0.181982, k = 0.000000,
         t1 = -1.822918, t2 = 2.004900, ax = 10.00, ay = 10.00;
bs11   : bending, l = 0.027092, t = 0.182602, k = 0.000000,
         t1 = -2.004900, t2 = 2.187502, ax = 10.00, ay = 10.00;

sxx_mh   : sextupole, l =   0.050000, k =    164.377175,
                  n =4, ax = 10.0, ay = 10.0;
sxy_mh   : sextupole, l =   0.050000, k =   -188.749729,
                  n =4, ax = 10.0, ay = 10.0;
syy_mh   : sextupole, l =   0.050000, k =    185.724436,
                  n =4, ax = 10.0, ay = 10.0;
sdx      : sextupole, l =   0.100000, k =   -243.051964,
                  n =4, ax = 10.0, ay = 10.0;
sfxh     : sextupole, l =   0.050000, k =    333.224049,
                  n =4, ax = 10.0, ay = 10.0;
sd       : sextupole, l =   0.100000, k =   -297.731747,
                  n =4, ax = 10.0, ay = 10.0;
sfh      : sextupole, l =   0.050000, k =    387.461814,
                  n =4, ax = 10.0, ay = 10.0;


bnom   : opticsmarker, betax = 0.305892, alphax = 0.000000,
         betay = 7.108460, alphay = 0.000000, etax  = 0.001279,
         etaxp  = 0.000000, etay  = 0.000000, etayp  = 0.000000, ax = 10.00,
         ay = 10.00;
bsom   : opticsmarker, betax = 0.305320, alphax = -0.000001,
         betay = 7.126130, alphay = 0.000000, etax  = -0.000012,
         etaxp  = 0.000024, etay  = 0.000000, etayp  = 0.000000, ax = 10.00,
         ay = 10.00;

vb_bs  : combined, l = 0.206100, t = 0.312499+aban, k = -3.856634,
         t1 = -2.187501, t2 = 2.5+aban, ax = 10.00, ay = 10.00;
an_bs  : combined, l = 0.300000, t = -aban, k = 3.918671, t1 = -aban,
         t2 = 0.000000, ax = 10.00, ay = 10.00;

an       : combined, l =   0.300000, t =  -0.780000, k =   3.905239,
       t1 =  -0.780000, t2 =   0.000000, ax = 10.0, ay = 10.0;
anm      : combined, l =   0.300000, t =  -0.780000, k =   3.512432,
       t1 =  -0.780000, t2 =   0.000000, ax = 10.0, ay = 10.0;
vbm      : combined, l =   0.206100, t =   1.092499, k =  -2.162146,
       t1 =  -2.187496, t2 =   3.280000, ax = 10.0, ay = 10.0;
vb       : combined, l =   0.206100, t =   1.092499, k =  -3.795373,
       t1 =  -2.187496, t2 =   3.280000, ax = 10.0, ay = 10.0;


xm     : photonbeam, xl = 2.20, style = 1, snap = 1, ax = 10.00, ay = 10.00;
xo     : photonbeam, xl = 1.50, style = 2, snap = 1, ax = 10.00, ay = 10.00;
xs     : photonbeam, xl = 3.00, style = 0, snap = 1, ax = 50.00, ay = 50.00;

lgbs   : girder, typ = 3, shift = 0.00000, ax = 10.00, ay = 10.00;
lgbn   : girder, typ = 3, shift = 0.00000, ax = 10.00, ay = 10.00;
lgbe   : girder, typ = 3, shift = 0.00000, ax = 10.00, ay = 10.00;

oxx_m    : multipole, n = 4, k =     54.535976,
                  ax = 10.0, ay = 10.0;
oxy_m    : multipole, n = 4, k =    122.809992,
                  ax = 10.0, ay = 10.0;
oyy_m    : multipole, n = 4, k =   -175.832660,
                  ax = 10.0, ay = 10.0;
ocxx     : multipole, n = 4, k =   -150.028357,
                  ax = 10.0, ay = 10.0;
ocxx2    : multipole, n = 4, k =     23.544000,
                  ax = 10.0, ay = 10.0;


mon    : monitor, ax = 10.00, ay = 10.00;

ch     : h-corrector, ax = 10.00, ay = 10.00;

cv     : v-corrector, ax = 10.00, ay = 10.00;


{----- table of segments ----------------------------------------------------}

dss    : dss0, mgp;
bs     : bsom, bs00, bs01, bs02, bs03, bs04, bs05, bs06, bs07, bs08,
         bs09, bs10, bs11;
bsl    : bs, lgbs;
bn     : bnom, bn00, bn01, bn02, bn03, bn04, bn05, bn06;
bnl    : bn, lgbn;
bnsup  : lgbe, bn, lgbe;
cc     : ch, cv;
sxx_m  : sxx_mh, cc, sxx_mh;
sxy_m  : sxy_mh, cc, sxy_mh;
syy_m  : syy_mh, cc, syy_mh;
hnc    : bnl, dnvb, vb, dnm1, sd, dnm2, an;
hsc    : bsl, dbs, dnvb, vb_bs, dnm1, sd, dnm2, an_bs;
mpfho  : doc, ocxx2, doc, sfh;
mpfh_  : doc, mon, doc, sfh;
mpf    : mpfh_, cc, -mpfho;
ncell  : -mpfho, -hnc, hnc, mpfh_;
ncellr : hnc, mpfh_, -mpfho, -hnc;
scell  : -mpfho, -hsc, hsc, mpfh_;
mpfx   : doc, mon, doc, sfxh, cc, sfxh, doc, ocxx, doc;
hcor   : xm, hnc, mpf, -hnc, xo, hnc, mpf, -hnc, xo;
hscor  : hsc, mpf, -hnc, hnc, mpf, -hnc;
moxy_m : doc, oxy_m, doc, sxy_m, doc;
moxx_m : doc, oxx_m, doc, sxx_m, doc;
moyy_m : doc, oyy_m, doc, syy_m, doc;
dsupm  : hnc, mpfx, anm, dxm, sdx, dnm1, -vbm, dnvb, -bnsup;
matm   : moyy_m, dmon, mon, dme, qm1, dm1, dmon, mon, moxy_m, qm2,
         moxx_m, dmon, mon, dm2, qm3, dms;
tm     : dsupm, matm;
hm     : hcor, tm, nper=24;
hsm    : hscor, tm;
marc   : hm, -hm;
per    : center, xs, -hm, hm, xs, -hsm, hsm, xs, -hm, hm, xs, -hm, hm,
         nper=3;
arcn   : center, xs, -hm, hm, nper=12;
ring   : 12*arcn;

{c:\users\streun\opadat\sls-2\dc12c.opa}
