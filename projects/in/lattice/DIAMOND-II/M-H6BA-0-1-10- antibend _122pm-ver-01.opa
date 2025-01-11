{..bs36\desktop\m-h6ba-3-1-1\m-h6ba-0-1-3 - antibend -ver-05_150pm-03.opa}


energy = 3.500000;

    betax   = 2.4571816; alphax  = 0.0000000;
    etax    = 0.0203992; etaxp   = 0.0000000;
    betay   = 2.0175939; alphay  = 0.0000000;
    etay    = 0.0000000; etayp   = 0.0000000;

{----- variables ----------------------------------------------------}

aban1  = 0.035;
aban2  = 0;
aban3  = 0.315;
aban4  = 0.32;
rk1    = 5.26028347;
rk2    =-3.62182625;
rk3    = 6.12183082;
rk4    = 6.22723242;

{----- table of elements ----------------------------------------------------}

dr_01      : drift, l = 2.450000, ax = 50.00, ay = 50.00;
dr_02      : drift, l = 0.075000, ax = 50.00, ay = 50.00;
dr_03      : drift, l = 0.260000, ax = 50.00, ay = 50.00;
dr_04      : drift, l = 0.075000, ax = 50.00, ay = 50.00;
dr_05      : drift, l = 0.075000, ax = 50.00, ay = 50.00;
dr_06      : drift, l = 0.075000, ax = 50.00, ay = 50.00;
dr_07      : drift, l = 0.450000, ax = 50.00, ay = 50.00;
dr_08      : drift, l = 0.075000, ax = 50.00, ay = 50.00;
dr_09      : drift, l = 0.075000, ax = 50.00, ay = 50.00;
dr_10      : drift, l = 0.100000, ax = 50.00, ay = 50.00;
dr_11      : drift, l = 0.340000, ax = 50.00, ay = 50.00;
dr_12      : drift, l = 0.075000, ax = 50.00, ay = 50.00;
dr_13      : drift, l = 0.075000, ax = 50.00, ay = 50.00;
dr_14      : drift, l = 0.475000, ax = 50.00, ay = 50.00;
dr_15      : drift, l = 0.090000, ax = 50.00, ay = 50.00;
dr_16      : drift, l = 0.090000, ax = 50.00, ay = 50.00;
dr_17      : drift, l = 0.075000, ax = 50.00, ay = 50.00;
dr_18      : drift, l = 1.450000, ax = 50.00, ay = 50.00;
dl_01      : drift, l = 1.450000, ax = 50.00, ay = 50.00;
dl_02      : drift, l = 0.075000, ax = 50.00, ay = 50.00;
dl_03      : drift, l = 0.090000, ax = 50.00, ay = 50.00;
dl_04      : drift, l = 0.090000, ax = 50.00, ay = 50.00;
dl_05      : drift, l = 0.475000, ax = 50.00, ay = 50.00;
dl_06      : drift, l = 0.075000, ax = 50.00, ay = 50.00;
dl_07      : drift, l = 0.075000, ax = 50.00, ay = 50.00;
dl_08      : drift, l = 0.340000, ax = 50.00, ay = 50.00;
dl_09      : drift, l = 0.100000, ax = 50.00, ay = 50.00;
dl_10      : drift, l = 0.075000, ax = 50.00, ay = 50.00;
dl_11      : drift, l = 0.075000, ax = 50.00, ay = 50.00;
dl_12      : drift, l = 0.450000, ax = 50.00, ay = 50.00;
dl_13      : drift, l = 0.075000, ax = 50.00, ay = 50.00;
dl_14      : drift, l = 0.075000, ax = 50.00, ay = 50.00;
dl_15      : drift, l = 0.075000, ax = 50.00, ay = 50.00;
dl_16      : drift, l = 0.230000, ax = 50.00, ay = 50.00;
dl_17      : drift, l = 0.075000, ax = 50.00, ay = 50.00;
dl_18      : drift, l = 0.156000, ax = 50.00, ay = 50.00;
dl_19      : drift, l = 3.900000, ax = 50.00, ay = 50.00;

cellcenter : marker, ax = 50.00, ay = 50.00;
sextmark   : marker, ax = 50.00, ay = 50.00;
ms         : marker, ax = 50.00, ay = 50.00;
ss         : marker, ax = 50.00, ay = 50.00;
ls         : marker, ax = 50.00, ay = 50.00;



qf1        : quadrupole, l = 0.155000, k = 7.21705738, ax = 10.00, ay = 10.00;
qd2        : quadrupole, l = 0.102000, k =-9.62859996, ax = 10.00, ay = 10.00;
qd5        : quadrupole, l = 0.136000, k =-6.23861195, ax = 10.00, ay = 10.00;
qd2_c1     : quadrupole, l = 0.162000, k = 0.04375929, ax = 10.00, ay = 10.00;
qf1_c1     : quadrupole, l = 0.215000, k =-5.34266058, ax = 10.00, ay = 10.00;
quad_add   : quadrupole, l = 0.162000, k = 7.60144713 , ax = 10.00, ay = 10.00;

sf1        : sextupole, l = 0.070000, k = 258.09567592 , n = 4,
             ax = 10.00, ay = 10.00;
sd1        : sextupole, l = 0.140000, k =-275.23969438, n = 4,
             ax = 10.00, ay = 10.00;
sd2        : sextupole, l = 0.140000, k =-198.42458308, n = 4,
             ax = 50.00, ay = 10.00;
sh1        : sextupole, l = 0.050000, k = 0.000000, n = 4,
             ax = 10.00, ay = 10.00;
sh2        : sextupole, l = 0.050000, k = 0.000000, n = 4,
             ax = 10.00, ay = 10.00;

bnom1      : opticsmarker, betax = 0.205000, alphax = -0.125500,
             betay = 7.245000, alphay = -2.277000, etax  = 0.000000,
             etaxp  = 0.000000, etay  = 0.000000, etayp  = 0.000000,
             ax = 10.00, ay = 10.00;

qf4        : combined, l = 0.142000, t = -aban1, k = rk1,
             t1 = -aban1/2.0, t2 = -aban1/2.0, ax = 10.00, ay = 10.00;
qd3        : combined, l = 0.100000, t = -aban2, k = rk2,
             t1 = -aban2/2.0, t2 = -aban2/2.0, ax = 10.00, ay = 10.00;
qf8        : combined, l = 0.278000, t = -aban3, k = rk3,
             t1 = -aban3/2.0, t2 = -aban3/2.0, ax = 10.00, ay = 10.00;
qf6        : combined, l = 0.362000, t = -aban4, k = rk4,
             t1 = -aban4/2.0, t2 = -aban4/2.0, ax = 10.00, ay = 10.00;
dq1        : combined, l = 0.869804, t =3.11617749, k =-2.75325131, t1 =1.55808875, t2 =1.55808875, ax = 10.00, ay = 10.00;



bl1_1      : combined, l = 0.200983, t =0.33481650, k =0.02346832,
             t1 =0.000000, t2 =0.33481650, ax = 10.00, ay = 10.00;
bl1_2      : combined, l = 0.200983, t = 0.38540370, k =0.02346832,
             t1 =-0.33481650, t2 =0.72022020, ax = 10.00, ay = 10.00;
bl1_3      : combined, l = 0.200983, t =0.45119986, k =0.02346832,
             t1 =-0.38540370, t2 =1.17142006, ax = 10.00, ay = 10.00;
bl1_4      : combined, l = 0.200983, t =0.55781748, k =0.02346832,
             t1 =-0.45119986, t2 =1.72923753, ax = 10.00, ay = 10.00;
bl1_5      : combined, l = 0.200983, t =0.85292543, k =0.02346832,
             t1 =-0.55781748, t2 =2.58216296, ax = 10.00, ay = 10.00;

bl2_1      : combined, l = 0.200983, t =0.32502612, k =0.22166965,
             t1 =0.00000000, t2 =0.32502612, ax = 10.00, ay = 10.00;
bl2_2      : combined, l = 0.200983, t =0.37413411, k =0.22166965,
             t1 =-0.32502612, t2 =0.69916023, ax = 10.00, ay = 10.00;
bl2_3      : combined, l = 0.200983, t =0.43800632, k =0.22166965,
             t1 =-0.37413411, t2 = 1.13716656, ax = 10.00, ay = 10.00;
bl2_4      : combined, l = 0.200983, t =0.54150633, k =0.22166965,
             t1 =-0.43800632, t2 =1.67867289, ax = 10.00, ay = 10.00;
bl2_5      : combined, l = 0.200983, t =0.82798503, k =0.22166965,
             t1 =-0.54150633, t2 =2.50665792, ax = 10.00, ay = 10.00;

of1s       : multipole, n = 4,  k = 0.00000000, ax = 10.00, ay = 10.00;

bpm        : monitor, ax = 50.00, ay = 50.00;

ch         : h-corrector, ax = 50.00, ay = 50.00;

cv         : v-corrector, ax = 50.00, ay = 50.00;


{----- table of segments ----------------------------------------------------}

bl1      : bl1_5, bl1_4, bl1_3, bl1_2, bl1_1;
bl2      : bl2_5, bl2_4, bl2_3, bl2_2, bl2_1;
ss_mcell : dr_01, qf1, dr_02, sh1, dr_03, qd2, dr_04;
dip_cell : bl1, dr_05, qd3, dr_06, sd1, dr_07, qf4, dr_08, sf1,
           sextmark, sf1, dr_09, qf4, dr_10, of1s, dr_11, sd2, dr_12, qd5,
           dr_13, -bl2;
ms_mcell : dr_14, qf6, dr_15, dq1, dr_16, qf8, dr_17, sh2, dr_18;
ls_mcell : dl_15, qd2_c1, dl_16, sh1, dl_17, qf1_c1, dl_18, quad_add, dl_19;
std_cell : ms, -ms_mcell, -dip_cell, -ss_mcell, ss, ss_mcell,
           dip_cell, ms_mcell, ms;
ls_hcell : -ms_mcell, -dip_cell, ls_mcell, ls;
sp_short : ls, -ls_hcell, std_cell, ls_hcell, ls;
sp       : ls, -ls_hcell, 3*std_cell, ls_hcell, ls;

{..bs36\desktop\m-h6ba-3-1-1\m-h6ba-0-1-3 - antibend -ver-05_150pm-03.opa}
