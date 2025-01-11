{..dt_energy_period_bendlengthangle_sextchromharm_newwpx921142y315311.opa}
{com bii standard user com}


energy = 1.200000;

    betax   = 19.5406881; alphax  = 0.0000000;
    etax    = -0.0000008; etaxp   = 0.0000000;
    betay   = 4.7290248; alphay  = 0.0000000;
    etay    = 0.0000000; etayp   = 0.0000000;

{----- variables ----------------------------------------------------}


{----- table of elements ----------------------------------------------------}

dq1     : drift, l = 0.288000, ax = 50.00, ay = 50.00;
ds1     : drift, l = 0.160000, ax = 50.00, ay = 50.00;
dq2     : drift, l = 0.260000, ax = 50.00, ay = 50.00;
ds2     : drift, l = 0.307000, ax = 50.00, ay = 50.00;
db      : drift, l = 0.260000, ax = 50.00, ay = 50.00;
ds3     : drift, l = 0.153000, ax = 50.00, ay = 50.00;
dq      : drift, l = 0.153000, ax = 50.00, ay = 50.00;
dl      : drift, l = 2.806000, ax = 50.00, ay = 50.00;
dk      : drift, l = 2.453000, ax = 50.00, ay = 50.00;

q1      : quadrupole, l = 0.250000, k = 2.425246, ax = 50.00, ay = 50.00;
q2      : quadrupole, l = 0.200000, k = -1.652303, ax = 50.00, ay = 50.00;
q3d     : quadrupole, l = 0.250000, k = -1.761750, ax = 50.00, ay = 50.00;
q4d     : quadrupole, l = 0.500000, k = 1.404000, ax = 50.00, ay = 50.00;
q3t     : quadrupole, l = 0.250000, k = -2.273750, ax = 50.00, ay = 50.00;
q4t     : quadrupole, l = 0.500000, k = 2.636000, ax = 50.00, ay = 50.00;
q5t     : quadrupole, l = 0.200000, k = -2.625000, ax = 50.00, ay = 50.00;

b       : bending, l = 1.175000, t = 22.500000, k = 0.000000,
          t1 = 11.250000, t2 = 11.250000, gap = 20.0000, ax = 50.00,
          ay = 50.00;

s1      : sextupole, l = 0.105000, k = 12.528145, n = 1, ax = 50.00,
          ay = 50.00;
s2      : sextupole, l = 0.160000, k = -10.516756, n = 1, ax = 50.00,
          ay = 50.00;
s3d     : sextupole, l = 0.160000, k = -10.184173, n = 1, ax = 50.00,
          ay = 50.00;
s4d     : sextupole, l = 0.160000, k = 6.971128, n = 1, ax = 50.00,
          ay = 50.00;
s3t     : sextupole, l = 0.160000, k = -13.795027, n = 1, ax = 50.00,
          ay = 50.00;
s4t     : sextupole, l = 0.160000, k = 20.834699, n = 1, ax = 50.00,
          ay = 50.00;

t2match : opticsmarker, betax = 1.000000, alphax = 0.000000,
          betay = 1.000000, alphay = 0.000000, etax  = 0.000000,
          etaxp  = 0.000000, etay  = 0.000000, etayp  = 0.000000, ax = 50.00,
          ay = 50.00;


{----- table of segments ----------------------------------------------------}

drh   : dl, s4d, dq, q4d, ds3, s3d, dq, q3d, db;
achlh : b, dq2, q2, ds2, s2, dq1, q1, ds1, s1;
tlh   : db, q3t, dq, s3t, ds3, q4t, dq, s4t, dq, q5t, dk;
cell  : drh, achlh, -achlh, tlh, -tlh, achlh, -achlh, -drh;
ring  : 4*cell;

{..dt_energy_period_bendlengthangle_sextchromharm_newwpx921142y315311.opa}
