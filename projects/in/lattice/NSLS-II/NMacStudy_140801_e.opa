{ ..cuments\accelerators\newmachinestudies\lattices\nmacstudy_140801_e.opa }
 
{----- global parameters (units: gev, m, rad) -------------------------------}
 
 
energy = 3.000000;
 
    betax   = 3.3958627; alphax  = 0.0000000;
    etax    = -0.0000046; etaxp   = 0.0000000;
    betay   = 2.9417282; alphay  = 0.0000000;
    etay    = 0.0000000; etayp   = 0.0000000;
    orbitx  = 0.0000000000; orbitxp = 0.0000000000;
    orbity  = 0.0000000000; orbityp = 0.0000000000;
    orbitdpp= 0.0000000000;
 
{----- table of elements (units: m, m^-2, deg, t; mm, mrad) ---------------- }
{      conventions: quadrupole: k>0 horizontally focusing                    }
{                   sextupole : k=m*l, m:=bpoletip/r^2/(b*rho)               }
 
l1       : drift, l = 0.050000, ax = 4.00, ay = 4.00;
l2       : drift, l = 0.100000, ax = 4.00, ay = 4.00;
l4       : drift, l = 0.154311, ax = 4.00, ay = 4.00;
l5       : drift, l = 0.756190, ax = 4.00, ay = 4.00;
ls       : drift, l = 2.500000, ax = 4.00, ay = 4.00;
l6       : drift, l = 0.173080, ax = 4.00, ay = 4.00;
l3h      : drift, l = 0.030029, ax = 4.00, ay = 4.00;
qf       : quadrupole, l = 0.075000, k = 21.889700, ax = 4.00, ay = 4.00;
qfe      : quadrupole, l = 0.100000, k = 23.402022, ax = 4.00, ay = 4.00;
qde      : quadrupole, l = 0.100000, k = -19.788989, ax = 4.00, ay = 4.00;
qm       : quadrupole, l = 0.150000, k = 18.290790, ax = 4.00, ay = 4.00;
bh       : bending, l = 0.166667, t = 0.500000, k = -7.009950, 
           t1 = 0.000000, t2 = 0.000000, gap = 0.00, 
           k1in = 0.0000, k1ex = 0.0000, k2in = 0.0000, 
           k2ex = 0.0000, ax = 4.00, ay = 4.00;
bm       : bending, l = 0.166667, t = 0.500000, k = -3.000000, 
           t1 = 0.000000, t2 = 0.000000, gap = 0.00, 
           k1in = 0.0000, k1ex = 0.0000, k2in = 0.0000, 
           k2ex = 0.0000, ax = 4.00, ay = 4.00;
sd       : sextupole, l = 0.100000, k = -1971.546180, n =5, 
           ax = 4.00, ay = 4.00;
sf       : sextupole, l = 0.050000, k = 3356.914090, n =5, 
           ax = 4.00, ay = 4.00;
sde      : sextupole, l = 0.100000, k = 0.000000, n =1, 
           ax = 4.00, ay = 4.00;
sfe      : sextupole, l = 0.100000, k = 0.000000, n =1, 
           ax = 4.00, ay = 4.00;
sga      : sextupole, l = 0.000000, k = 0.000000, n =1, 
           ax = 4.00, ay = 4.00;
sgb      : sextupole, l = 0.000000, k = 0.000000, n =1, 
           ax = 4.00, ay = 4.00;
sgc      : sextupole, l = 0.000000, k = 0.000000, n =1, 
           ax = 4.00, ay = 4.00;
sgd      : sextupole, l = 0.000000, k = 0.000000, n =1, 
           ax = 4.00, ay = 4.00;
sge      : sextupole, l = 0.000000, k = 0.000000, n =1, 
           ax = 4.00, ay = 4.00;
sd1      : sextupole, l = 0.100000, k = -6611.404500, n =1, 
           ax = 4.00, ay = 4.00;
sf1      : sextupole, l = 0.050000, k = 3326.180000, n =1, 
           ax = 4.00, ay = 4.00;
sd2      : sextupole, l = 0.100000, k = -6611.404500, n =1, 
           ax = 4.00, ay = 4.00;
sf2      : sextupole, l = 0.050000, k = 3326.180000, n =1, 
           ax = 4.00, ay = 4.00;
sd3      : sextupole, l = 0.100000, k = -6611.404500, n =1, 
           ax = 4.00, ay = 4.00;
sf3      : sextupole, l = 0.050000, k = 3326.180000, n =1, 
           ax = 4.00, ay = 4.00;
sd4      : sextupole, l = 0.100000, k = -6611.404500, n =1, 
           ax = 4.00, ay = 4.00;
sf4      : sextupole, l = 0.050000, k = 5027.292930, n =1, 
           ax = 4.00, ay = 4.00;
sd5      : sextupole, l = 0.100000, k = -6611.404500, n =1, 
           ax = 4.00, ay = 4.00;
sf5      : sextupole, l = 0.050000, k = 3326.180000, n =1, 
           ax = 4.00, ay = 4.00;
sd6      : sextupole, l = 0.100000, k = -6611.404500, n =1, 
           ax = 4.00, ay = 4.00;
sf6      : sextupole, l = 0.050000, k = 3326.180000, n =1, 
           ax = 4.00, ay = 4.00;
sd7      : sextupole, l = 0.100000, k = -6611.404500, n =1, 
           ax = 4.00, ay = 4.00;
sf7      : sextupole, l = 0.050000, k = 3326.180000, n =1, 
           ax = 4.00, ay = 4.00;
sd8      : sextupole, l = 0.100000, k = -6611.404500, n =1, 
           ax = 4.00, ay = 4.00;
sf8      : sextupole, l = 0.050000, k = 3326.180000, n =1, 
           ax = 4.00, ay = 4.00;
o1       : multipole, n = 4, k = 0.000, ax = 4.00, ay = 4.00;
o2       : multipole, n = 4, k = 0.000, ax = 4.00, ay = 4.00;
o3       : multipole, n = 4, k = 0.000, ax = 4.00, ay = 4.00;
 
{----- table of segments ----------------------------------------------------}
 
uch       : sf, l1, qf, l2, sd, bh;
uc        : uch, -uch;
ucb       : -uch, uch;
mcb       : ls, o1, sga, qfe, sgb, l3h, o2, l3h, qde, sgc, l4, sfe,
            o3, bm, sde, l5, sgd, qm, l6, sge, bh;
mcblinear : ls, qfe, l3h, l3h, qde, l4, sfe, bm, sde, l5, qm, l6,
            bh;
uch1      : sf1, l1, qf, l2, sd, bh;
ucb1      : -uch1, uch1;
uch2      : sf2, l1, qf, l2, sd, bh;
ucb2      : -uch2, uch2;
uch3      : sf3, l1, qf, l2, sd, bh;
ucb3      : -uch3, uch3;
uch4      : sf4, l1, qf, l2, sd, bh;
ucb4      : -uch1, uch1;
uch5      : sf5, l1, qf, l2, sd, bh;
ucb5      : -uch5, uch5;
uch6      : sf6, l1, qf, l2, sd, bh;
ucb6      : -uch6, uch6;
uch7      : sf7, l1, qf, l2, sd, bh;
ucb7      : -uch7, uch7;
uch8      : sf8, l1, qf, l2, sd, bh;
ucb8      : -uch8, uch8;
hper      : mcb, 8*ucb;
per       : hper, -hper, nper=20;
hpera     : mcb, ucb1, ucb2, ucb3, ucb4, ucb5, ucb6, ucb7, ucb8;
pera      : hpera, -hpera, nper=20;
perucs    : ucb, nper=360;
 
{ ..cuments\accelerators\newmachinestudies\lattices\nmacstudy_140801_e.opa }
