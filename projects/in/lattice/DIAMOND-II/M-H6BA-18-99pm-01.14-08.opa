energy = 3.500000;

c0 = 2.99792458e8; gamma = 3.5e9/c0;
h_rf = 934;

{----- table of elements ----------------------------------------------------}

{ CDR: Bare Lattice 1.66 MV, with IDs 2.69.
       V_RF = 1.5 MV => delta^ = 4.0%, sigma_s = 3.0 mm. }
cav: Cavity, L = 0.0, Frequency = c0/C*h_rf, Voltage = 1.5e6, harnum = h_rf,
     phi = 0.0;

dcor     : drift, l = 0.000000;
dcor1    : drift, l = 0.040000;
dL1 = 0.05;
dL2 = 0.1;
dr_01    : drift, l = 2.521350 - dL1 - dL2;
dr_02    : drift, l = 0.075000;
dr_03    : drift, l = 0.075000;
dr_04    : drift, l = 0.150000 + dL2;
dr_05    : drift, l = 0.075000;

d_14     : drift, l = 0.075000;
d_13     : drift, l = 0.075000;
d_12     : drift, l = 0.165500;
d_11     : drift, l = 0.075000;
d_10     : drift, l = 0.075000;
d_09     : drift, l = 0.075000;
d_08     : drift, l = 0.100000;
d_07     : drift, l = 0.080000;
d_06     : drift, l = 0.197000;
d_05     : drift, l = 0.090000;
d_04     : drift, l = 0.090000;
d_03     : drift, l = 0.075000;
d_02     : drift, l = 0.075000;
d_01     : drift, l = 1.384200;

dl_27    : drift, l = 3.696350;
dl_26    : drift, l = 0.075000;
dl_25    : drift, l = 0.035000;
dl_24    : drift, l = 0.035000;
dl_23    : drift, l = 0.075000;
dl_22    : drift, l = 0.075000;
dl_21    : drift, l = 0.075000;
dl_20    : drift, l = 0.075000;

b_t      : marker;
ms       : marker;
ss       : marker;
ls       : marker;

dl1a_1:   bending, l =  0.20000000,
    t =  0.40523456, t1 =  0.00000000, t2 =  0.40523456, k =  0.03990854;
dl1a_2:   bending, l =  0.20000000,
    t =  0.46646119, t1 = -0.40523456, t2 =  0.87169575, k =  0.03990854;
dl1a_3:   bending, l =  0.20000000,
    t =  0.54609549, t1 = -0.46646119, t2 =  1.41779124, k =  0.03990854;
dl1a_4:   bending, l =  0.20000000,
    t =  0.67513675, t1 = -0.54609549, t2 =  2.09292799, k =  0.03990854;
dl1a_5:   bending, l =  0.20000000,
    t =  1.03231134, t1 = -0.67513675, t2 =  3.12523933, k =  0.03990854;
dl2a_1:   bending, l =  0.16667000,
    t =  0.35814357, t1 =  0.00000000, t2 =  0.35814357, k = -0.36838095;
dl2a_2:   bending, l =  0.16667000,
    t =  0.41225525, t1 = -0.35814357, t2 =  0.77039882, k = -0.36838095;
dl2a_3:   bending, l =  0.33334000,
    t =  0.48263551, t1 = -0.41225525, t2 =  1.25303434, k = -0.36838095;
dl2a_4:   bending, l =  0.16667000,
    t =  0.59668131, t1 = -0.48263551, t2 =  1.84971564, k = -0.36838095;
dl2a_5:   bending, l =  0.16667000,
    t =  0.91234980, t1 = -0.59668131, t2 =  2.76206544, k = -0.36838095;

qf4:      bending, l =  0.15000000,
    t = -0.26178834, t1 = -0.13089417, t2 = -0.13089417, k =  4.84724428;
qf8:      bending, l =  0.25000000,
    t = -0.36746657, t1 = -0.18373329, t2 = -0.18373329, k =  6.81570447;
dq1:      bending, l =  0.87000000,
    t =  2.50373848, t1 =  1.25186924, t2 =  1.25186924, k = -2.77224175;

qd3:      quadrupole, l = 0.15000000, k =   -3.62025211;
qd5:      quadrupole, l = 0.10500000, k =   -4.89174834;
qf6:      quadrupole, l = 0.36000000, k =    7.04754125;

{qf1:      quadrupole, l = 0.20000000, k =    8.43607057;
qd2:      quadrupole, l = 0.15000000, k =   -7.56247065;
qf1_c1:   quadrupole, l = 0.18500000, k =   -6.03911106;
qd2_c1:   quadrupole, l = 0.10500000, k =    0.61175899;
quad_add: quadrupole, l = 0.18500000, k =    7.34101051;}

{ nu = [64.87, 19.70]. }
{qf1:      quadrupole, l = 0.20000000, k =    8.36359597;
qd2:      quadrupole, l = 0.15000000, k =   -7.52234783;
qf1_c1:   quadrupole, l = 0.18500000, k =   -6.03679910;
qd2_c1:   quadrupole, l = 0.10500000, k =    0.63005007;
quad_add: quadrupole, l = 0.18500000, k =    7.30778719;}

{ nu = [64.88, 19.68]. }
qf1:      quadrupole, l = 0.20000000, k =    8.38301056;
qd2:      quadrupole, l = 0.15000000, k =   -7.50203775;
qf1_c1:   quadrupole, l = 0.18500000, k =   -6.02996775;
qd2_c1:   quadrupole, l = 0.10500000, k =    0.64010006;
quad_add: quadrupole, l = 0.18500000, k =    7.31734861;

{ nu = [64.90, 19.67]. }
{qf1:      quadrupole, l = 0.20000000, k =    8.38702440;
qd2:      quadrupole, l = 0.15000000, k =   -7.50259118;
qf1_c1:   quadrupole, l = 0.18500000, k =   -6.02968657;
qd2_c1:   quadrupole, l = 0.10500000, k =    0.63989567;
quad_add: quadrupole, l = 0.18500000, k =    7.31925836;}

{ nu = [64.91, 19.66]. }
{qf1:      quadrupole, l = 0.20000000, k =    8.38900830;
qd2:      quadrupole, l = 0.15000000, k =   -7.50244659;
qf1_c1:   quadrupole, l = 0.18500000, k =   -6.02944742;
qd2_c1:   quadrupole, l = 0.10500000, k =    0.63999751;
quad_add: quadrupole, l = 0.18500000, k =    7.32020605;}

sf1:      sextupole,  l = 0.07000000, k =  133.51274409;
sd1:      sextupole,  l = 0.14000000, k = -113.85554836;
sd2:      sextupole,  l = 0.14000000, k = -131.81065897;
s:        sextupole,  l = 0.10000000, k = -121.81032466;
sh1:      sextupole,  l = 0.10000000, k =   64.28900179;
sh2:      sextupole,  l = 0.10000000, k = -201.07250629;
of_dh:    drift, l = 0.09/2.0;
of1_m:    multipole, n = 4, k = 2.21098729e+03*0.09;
of1: of_dh, of1_m, of_dh;

of2_m    : multipole, n = 4, k = 0.0;
of2: of_dh, of2_m, of_dh;

bpm: monitor;
ch:  h-corrector;
cv:  v-corrector;
chv: ch, cv;


{----- table of segments ----------------------------------------------------}

dl1a     : dl1a_5, dl1a_4, dl1a_3, dl1a_2, dl1a_1;
dl2a     : dl2a_5, dl2a_4, dl2a_3, dl2a_2, dl2a_1;

match_ss: dr_05, sh1, dr_03, chv, qd2, dr_04, qf1, dr_02, bpm, dr_01;

match_ms: d_01, bpm, d_02, chv, sh2, d_03, qf8, d_04;

match_ls: dl_20, sh1, dl_21, bpm, dl_22, chv, qd2_c1, dl_23,
          qf1_c1, dl_24, dcor1, chv, dcor1, dl_25, quad_add,
          dl_26, bpm, dl_27;

dip_cell: dq1, d_05, qf6, d_06, bpm, d_07, dcor, s, chv, dcor,
          d_08, dl2a, d_09, qd5, d_10, sd2, chv, d_11,
          bpm, d_12, of1, d_13, qf4, d_14, sf1,
          sf1, d_14, qf4, d_13, of2, d_12, bpm,
          d_11, chv, sd1, d_10, qd3, d_09, -dl1a;

ms_to_ss : match_ms, dip_cell, match_ss;

ms_to_ls : match_ms, dip_cell, match_ls;

std_cell : ms, ms_to_ss, ss, b_t, -ms_to_ss;
sp       : ls, b_t, -ms_to_ls, 3*std_cell, ms, ms_to_ls, ls;

ring: 1*sp, b_t{, cav};

cell: ring, nper=6;

end;