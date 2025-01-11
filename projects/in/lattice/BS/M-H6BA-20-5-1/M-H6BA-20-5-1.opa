{s:\technical\accelerator physics\ghasem-beni\m-h6ba-20-5-1.opa}


energy = 3.500000;

    betax   = 14.5713226; alphax  = 0.0000000;
    etax    = -0.0000215; etaxp   = 0.0000000;
    betay   = 2.3612895; alphay  = 0.0000000;
    etay    = 0.0000000; etayp   = 0.0000000;

{----- variables ----------------------------------------------------}

nsext  = 2;
aban1  = 0.20585295;
aban2  = -0.3853;
aban3  = 0.3579;
rk1    = 5.26028347;
rk2    = -3.62182625;
rk3    = 6.12183082;
rk4    = 6.22723242;

{----- table of elements ----------------------------------------------------}

dmult    : drift, l = 0.045000, ax = 12.50, ay = 12.50;
dcor     : drift, l = 0.000000, ax = 12.50, ay = 12.50;
dcor1    : drift, l = 0.040000, ax = 12.50, ay = 12.50;
dr_01    : drift, l = 2.521350, ax = 12.00, ay = 12.00;
dr_02    : drift, l = 0.075000, ax = 12.00, ay = 12.00;
dr_03    : drift, l = 0.075000, ax = 12.00, ay = 12.00;
dr_04    : drift, l = 0.150000, ax = 12.00, ay = 12.00;
dr_05    : drift, l = 0.075000, ax = 12.00, ay = 12.00;
dr_06    : drift, l = 0.075000, ax = 12.00, ay = 12.00;
dr_07    : drift, l = 0.075000, ax = 12.00, ay = 12.00;
dr_08    : drift, l = 0.075000, ax = 12.00, ay = 12.00;
dr_09    : drift, l = 0.330500, ax = 12.00, ay = 12.00;
dr_091   : drift, l = 0.075000, ax = 12.00, ay = 12.00;
dr_092   : drift, l = 0.165500, ax = 12.00, ay = 12.00;
dr_10    : drift, l = 0.075000, ax = 12.00, ay = 12.00;
dr_11    : drift, l = 0.075000, ax = 12.00, ay = 12.00;
dr_12    : drift, l = 0.075000, ax = 12.00, ay = 12.00;
dr_13    : drift, l = 0.165500, ax = 12.00, ay = 12.00;
dr_14    : drift, l = 0.075000, ax = 12.00, ay = 12.00;
dr_15    : drift, l = 0.075000, ax = 12.00, ay = 12.00;
dr_16    : drift, l = 0.075000, ax = 12.00, ay = 12.00;
dr_17    : drift, l = 0.100000, ax = 12.00, ay = 12.00;
dr_18    : drift, l = 0.080000, ax = 12.00, ay = 12.00;
dr_19    : drift, l = 0.197000, ax = 12.00, ay = 12.00;
dr_20    : drift, l = 0.090000, ax = 12.00, ay = 12.00;
dr_21    : drift, l = 0.090000, ax = 12.00, ay = 12.00;
dr_22    : drift, l = 0.075000, ax = 12.00, ay = 12.00;
dr_23    : drift, l = 0.075000, ax = 12.00, ay = 12.00;
dr_24    : drift, l = 1.384200, ax = 12.00, ay = 12.00;
dl_27    : drift, l = 3.696350, ax = 12.00, ay = 12.00;
dl_26    : drift, l = 0.075000, ax = 12.00, ay = 12.00;
dl_25    : drift, l = 0.035000, ax = 12.00, ay = 12.00;
dl_24    : drift, l = 0.035000, ax = 12.00, ay = 12.00;
dl_23    : drift, l = 0.075000, ax = 12.00, ay = 12.00;
dl_22    : drift, l = 0.075000, ax = 12.00, ay = 12.00;
dl_21    : drift, l = 0.075000, ax = 12.00, ay = 12.00;
dl_20    : drift, l = 0.075000, ax = 12.00, ay = 12.00;
dl_19    : drift, l = 0.075000, ax = 12.00, ay = 12.00;
dl_18    : drift, l = 0.075000, ax = 12.00, ay = 12.00;
dl_17    : drift, l = 0.075000, ax = 12.00, ay = 12.00;
dl_16    : drift, l = 0.330500, ax = 12.00, ay = 12.00;
dl_161   : drift, l = 0.075000, ax = 12.00, ay = 12.00;
dl_162   : drift, l = 0.165500, ax = 12.00, ay = 12.00;
dl_15    : drift, l = 0.075000, ax = 12.00, ay = 12.00;
dl_14    : drift, l = 0.075000, ax = 12.00, ay = 12.00;
dl_13    : drift, l = 0.075000, ax = 12.00, ay = 12.00;
dl_12    : drift, l = 0.165500, ax = 12.00, ay = 12.00;
dl_11    : drift, l = 0.075000, ax = 12.00, ay = 12.00;
dl_10    : drift, l = 0.075000, ax = 12.00, ay = 12.00;
dl_09    : drift, l = 0.075000, ax = 12.00, ay = 12.00;
dl_08    : drift, l = 0.100000, ax = 12.00, ay = 12.00;
dl_07    : drift, l = 0.080000, ax = 12.00, ay = 12.00;
dl_06    : drift, l = 0.197000, ax = 12.00, ay = 12.00;
dl_05    : drift, l = 0.090000, ax = 12.00, ay = 12.00;
dl_04    : drift, l = 0.090000, ax = 12.00, ay = 12.00;
dl_03    : drift, l = 0.075000, ax = 12.00, ay = 12.00;
dl_02    : drift, l = 0.075000, ax = 12.00, ay = 12.00;
dl_01    : drift, l = 1.384200, ax = 12.00, ay = 12.00;

ms       : marker, ax = 12.00, ay = 12.00;
ss       : marker, ax = 12.00, ay = 12.00;
ls       : marker, ax = 12.00, ay = 12.00;

qf1      : quadrupole, l = 0.150000, k = 10.860458, ax = 12.00, ay = 12.00;
qd2      : quadrupole, l = 0.150000, k = -7.292291, ax = 12.00, ay = 12.00;
qd3      : quadrupole, l = 0.150000, k = -2.411602, ax = 12.00, ay = 12.00;
qd5      : quadrupole, l = 0.105000, k = -3.317328, ax = 12.00, ay = 12.00;
qd3_c1   : quadrupole, l = 0.150000, k = -2.143560, ax = 12.00, ay = 12.00;
qd2_c1   : quadrupole, l = 0.105000, k = -3.137978, ax = 12.00, ay = 12.00;
qf1_c1   : quadrupole, l = 0.185000, k = -3.699175, ax = 12.00, ay = 12.00;
quad_add : quadrupole, l = 0.185000, k = 6.398149, ax = 12.00, ay = 12.00;

dl1a_5   : bending, l = 0.200000, t = 0.933264, k = -0.108016,
           t1 = -0.610359, t2 = 2.825380, ax = 12.00, ay = 12.00;
dl1a_4   : bending, l = 0.200000, t = 0.610359, k = -0.108016,
           t1 = -0.493699, t2 = 1.892117, ax = 12.00, ay = 12.00;
dl1a_3   : bending, l = 0.200000, t = 0.493699, k = -0.108016,
           t1 = -0.421705, t2 = 1.281758, ax = 12.00, ay = 12.00;
dl1a_2   : bending, l = 0.200000, t = 0.421705, k = -0.108016,
           t1 = -0.366353, t2 = 0.788059, ax = 12.00, ay = 12.00;
dl1a_1   : bending, l = 0.200000, t = 0.366353, k = -0.108016,
           t1 = 0.000000, t2 = 0.366353, ax = 12.00, ay = 12.00;
dl2a_5   : bending, l = 0.166670, t = 0.842316, k = -0.573710,
           t1 = -0.550879, t2 = 2.550044, ax = 12.00, ay = 12.00;
dl2a_4   : bending, l = 0.166670, t = 0.550879, k = -0.573710,
           t1 = -0.445588, t2 = 1.707728, ax = 12.00, ay = 12.00;
dl2a_3   : bending, l = 0.333340, t = 0.445588, k = -0.573710,
           t1 = -0.380610, t2 = 1.156849, ax = 12.00, ay = 12.00;
dl2a_2   : bending, l = 0.166670, t = 0.380610, k = -0.573710,
           t1 = -0.330652, t2 = 0.711262, ax = 12.00, ay = 12.00;
dl2a_1   : bending, l = 0.166670, t = 0.330652, k = -0.573710,
           t1 = 0.000000, t2 = 0.330652, ax = 12.00, ay = 12.00;

sd1      : sextupole, l = 0.140000, k = -187.275767, n = 1,
           ax = 12.00, ay = 12.00;
sd2      : sextupole, l = 0.140000, k = -141.444215, n = 1,
           ax = 12.00, ay = 12.00;
sf1      : sextupole, l = 0.140000, k = 168.051198, n = 1,
           ax = 12.00, ay = 12.00;
sh1      : sextupole, l = 0.100000, k = 22.700000, n = 1, ax = 12.00,
           ay = 12.00;
sh2      : sextupole, l = 0.100000, k = 19.900000, n = 1, ax = 12.00,
           ay = 12.00;
s        : sextupole, l = 0.100000, k = -14.200000, n = 1,
           ax = 12.00, ay = 12.00;

qf4      : combined, l = 0.150000, t = -0.205853, k = 4.698648,
           t1 = -0.102926, t2 = -0.102926, ax = 12.00, ay = 12.00;
qf4l     : combined, l = 0.150000, t = -0.205853, k = 4.698648,
           t1 = -0.102926, t2 = -0.102926, ax = 12.00, ay = 12.00;
qf4_c1   : combined, l = 0.150000, t = -0.205853, k = 4.637100,
           t1 = -0.102926, t2 = -0.102926, ax = 12.00, ay = 12.00;
qf6      : combined, l = 0.360000, t = 0.385300, k = 8.217016,
           t1 = 0.192650, t2 = 0.192650, ax = 12.00, ay = 12.00;
qf8      : combined, l = 0.250000, t = -0.357900, k = 6.262614,
           t1 = -0.178950, t2 = -0.178950, ax = 12.00, ay = 12.00;
dq1      : combined, l = 0.870000, t = 2.508881, k = -2.814485,
           t1 = 1.254441, t2 = 1.254441, ax = 12.00, ay = 12.00;

of1      : multipole, n = 4,  k = 95.00000000, ax = 12.00, ay = 12.00;

bpm_01   : monitor, ax = 12.50, ay = 12.50;
bpm_02   : monitor, ax = 12.50, ay = 12.50;
bpm_03   : monitor, ax = 12.50, ay = 12.50;
bpm_04   : monitor, ax = 12.50, ay = 12.50;
bpm_05   : monitor, ax = 12.50, ay = 12.50;
bpm_06   : monitor, ax = 12.50, ay = 12.50;
bpm_07   : monitor, ax = 12.50, ay = 12.50;
bpm_08   : monitor, ax = 12.50, ay = 12.50;
bpm_09   : monitor, ax = 12.50, ay = 12.50;
bpm_10   : monitor, ax = 12.50, ay = 12.50;
bpm_11   : monitor, ax = 12.50, ay = 12.50;

ch_01    : h-corrector, ax = 12.50, ay = 12.50;
ch_02    : h-corrector, ax = 12.50, ay = 12.50;
ch_03    : h-corrector, ax = 12.50, ay = 12.50;
ch_04    : h-corrector, ax = 12.50, ay = 12.50;
ch_05    : h-corrector, ax = 12.50, ay = 12.50;
ch_06    : h-corrector, ax = 12.50, ay = 12.50;
ch_07    : h-corrector, ax = 12.50, ay = 12.50;
ch_08    : h-corrector, ax = 12.50, ay = 12.50;
ch_09    : h-corrector, ax = 12.50, ay = 12.50;
ch_10    : h-corrector, ax = 12.50, ay = 12.50;
ch_11    : h-corrector, ax = 12.50, ay = 12.50;

cv_01    : v-corrector, ax = 12.50, ay = 12.50;
cv_02    : v-corrector, ax = 12.50, ay = 12.50;
cv_03    : v-corrector, ax = 12.50, ay = 12.50;
cv_04    : v-corrector, ax = 12.50, ay = 12.50;
cv_05    : v-corrector, ax = 12.50, ay = 12.50;
cv_06    : v-corrector, ax = 12.50, ay = 12.50;
cv_07    : v-corrector, ax = 12.50, ay = 12.50;
cv_08    : v-corrector, ax = 12.50, ay = 12.50;
cv_09    : v-corrector, ax = 12.50, ay = 12.50;
cv_10    : v-corrector, ax = 12.50, ay = 12.50;
cv_11    : v-corrector, ax = 12.50, ay = 12.50;


{----- table of segments ----------------------------------------------------}

dl1a      : dl1a_5, dl1a_4, dl1a_3, dl1a_2, dl1a_1;
dl2a      : dl2a_5, dl2a_4, dl2a_3, dl2a_2, dl2a_1;
of1s      : dmult, of1, dmult;
arca_c2r  : dr_24, bpm_05, dr_23, ch_05, cv_05, sh2, dr_22, qf8,
            dr_21, dq1, dr_20, qf6, dr_19, bpm_04, dr_18, dcor, s, ch_04,
            cv_04, dcor, dr_17, dl2a, dr_16, qd5, dr_15, sd2, ch_03, cv_03,
            dr_14, bpm_03, dr_13, of1s, dr_12, qf4, dr_11, sf1, dr_10, qf4,
            dr_091, of1s, dr_092, bpm_02, dr_08, ch_02, cv_02, sd1, dr_07, qd3,
            dr_06, -dl1a, dr_05, qd2, dr_04, ch_01, cv_01, sh1, dr_03, qf1,
            dr_02, bpm_01, dr_01;
arca_c1r  : dl_01, bpm_06, dl_02, ch_06, cv_06, sh2, dl_03, qf8,
            dl_04, dq1, dl_05, qf6, dl_06, bpm_07, dl_07, dcor, s, ch_07,
            cv_07, dcor, dl_08, dl2a, dl_09, qd5, dl_10, sd2, ch_08, cv_08,
            dl_11, bpm_08, dl_12, of1s, dl_13, qf4l, dl_14, sf1, dl_15, qf4_c1,
            dl_161, of1s, dl_162, bpm_09, dl_17, ch_09, cv_09, sd1, dl_18,
            qd3_c1, dl_19, -dl1a, dl_20, qd2_c1, dl_21, bpm_10, dl_22, ch_10,
            cv_10, sh1, dl_23, qf1_c1, dl_24, dcor1, ch_11, cv_11, dcor1,
            dl_25, quad_add, dl_26, bpm_11, dl_27;
arca_c2   : ms, arca_c2r, ss, -arca_c2r;
arca_c21  : arca_c2r, -arca_c2r;
arca_c21m : -arca_c2r, arca_c2r;
arca_c1   : arca_c2r, ms, -arca_c1r;
arca_c11  : -arca_c1r, arca_c1r;
sp        : ls, -arca_c1r, 3*arca_c2, arca_c1r, ls;
ring      : sp, nper=6;

{s:\technical\accelerator physics\ghasem-beni\m-h6ba-20-5-1.opa}
