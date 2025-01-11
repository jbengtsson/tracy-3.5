{ f:\work\opa\examples\mlssu.opa }
 
{----- global parameters (units: gev, m, rad) -------------------------------}
 
 
energy = 0.629000;
 
    betax   = 7.5998049; alphax  = 0.0000000;
    etax    = -0.0311487; etaxp   = 0.0000000;
    betay   = 2.1403700; alphay  = 0.0000000;
    etay    = 0.0000000; etayp   = 0.0000000;
    orbitx  = 0.0000000000; orbitxp = 0.0000000000;
    orbity  = 0.0000000000; orbityp = 0.0000000000;
    orbitdpp= 0.0000000000;
 
{----- table of elements (units: m, m^-2, deg, t; mm, mrad) ---------------- }
{      conventions: quadrupole: k>0 horizontally focusing                    }
{                   sextupole : k=m*l, m:=bpoletip/r^2/(b*rho)               }
 
d3101    : drift, l = 0.300000, ax = 50.00, ay = 50.00;
d4       : drift, l = 0.425000, ax = 50.00, ay = 50.00;
d2201    : drift, l = 0.300000, ax = 50.00, ay = 50.00;
d401     : drift, l = 0.300000, ax = 50.00, ay = 50.00;
d2202    : drift, l = 0.150000, ax = 50.00, ay = 50.00;
d3102    : drift, l = 0.150000, ax = 50.00, ay = 50.00;
d51      : drift, l = 0.125000, ax = 50.00, ay = 50.00;
d231     : drift, l = 0.125000, ax = 50.00, ay = 50.00;
d24      : drift, l = 0.150000, ax = 50.00, ay = 50.00;
d322     : drift, l = 0.125000, ax = 50.00, ay = 50.00;
d601     : drift, l = 0.300000, ax = 50.00, ay = 50.00;
d232     : drift, l = 0.125000, ax = 50.00, ay = 50.00;
d321     : drift, l = 0.125000, ax = 50.00, ay = 50.00;
d402     : drift, l = 0.150000, ax = 50.00, ay = 50.00;
d3302    : drift, l = 0.425000, ax = 50.00, ay = 50.00;
d3301    : drift, l = 0.300000, ax = 50.00, ay = 50.00;
d33      : drift, l = 0.150000, ax = 50.00, ay = 50.00;
d2401    : drift, l = 0.300000, ax = 50.00, ay = 50.00;
d52      : drift, l = 0.125000, ax = 50.00, ay = 50.00;
d8       : drift, l = 0.150000, ax = 50.00, ay = 50.00;
d25      : drift, l = 0.425000, ax = 50.00, ay = 50.00;
d6       : drift, l = 0.150000, ax = 50.00, ay = 50.00;
d2402    : drift, l = 0.425000, ax = 50.00, ay = 50.00;
d36      : drift, l = 0.150000, ax = 50.00, ay = 50.00;
d602     : drift, l = 0.425000, ax = 50.00, ay = 50.00;
d26      : drift, l = 0.150000, ax = 50.00, ay = 50.00;
d34      : drift, l = 0.425000, ax = 50.00, ay = 50.00;
d9       : drift, l = 0.150000, ax = 50.00, ay = 50.00;
d31      : drift, l = 0.425000, ax = 50.00, ay = 50.00;
d7       : drift, l = 0.425000, ax = 50.00, ay = 50.00;
d29      : drift, l = 0.150000, ax = 50.00, ay = 50.00;
d35      : drift, l = 0.150000, ax = 50.00, ay = 50.00;
d27      : drift, l = 0.150000, ax = 50.00, ay = 50.00;
end      : drift, l = 0.000000, ax = 50.00, ay = 50.00;
d101     : drift, l = 0.150000, ax = 50.00, ay = 50.00;
d12      : drift, l = 0.425000, ax = 50.00, ay = 50.00;
d2801    : drift, l = 0.150000, ax = 50.00, ay = 50.00;
d3601    : drift, l = 1.250000, ax = 50.00, ay = 50.00;
d2701    : drift, l = 3.000000, ax = 50.00, ay = 50.00;
d10      : drift, l = 3.000000, ax = 50.00, ay = 50.00;
d28      : drift, l = 3.000000, ax = 50.00, ay = 50.00;
d1       : drift, l = 1.250000, ax = 50.00, ay = 50.00;
d1001    : drift, l = 0.150000, ax = 50.00, ay = 50.00;
d901     : drift, l = 3.000000, ax = 50.00, ay = 50.00;
d2       : drift, l = 0.150000, ax = 50.00, ay = 50.00;
d11      : drift, l = 0.150000, ax = 50.00, ay = 50.00;
d30      : drift, l = 0.425000, ax = 50.00, ay = 50.00;
d3       : drift, l = 0.425000, ax = 50.00, ay = 50.00;
d13      : drift, l = 0.425000, ax = 50.00, ay = 50.00;
d1301    : drift, l = 0.300000, ax = 50.00, ay = 50.00;
d1302    : drift, l = 0.150000, ax = 50.00, ay = 50.00;
d141     : drift, l = 0.125000, ax = 50.00, ay = 50.00;
d15      : drift, l = 0.150000, ax = 50.00, ay = 50.00;
d1501    : drift, l = 0.300000, ax = 50.00, ay = 50.00;
d142     : drift, l = 0.125000, ax = 50.00, ay = 50.00;
d16      : drift, l = 0.425000, ax = 50.00, ay = 50.00;
d1502    : drift, l = 0.425000, ax = 50.00, ay = 50.00;
d19      : drift, l = 1.250000, ax = 50.00, ay = 50.00;
d1801    : drift, l = 1.250000, ax = 50.00, ay = 50.00;
d17      : drift, l = 0.150000, ax = 50.00, ay = 50.00;
d18      : drift, l = 0.150000, ax = 50.00, ay = 50.00;
d1901    : drift, l = 0.150000, ax = 50.00, ay = 50.00;
d21      : drift, l = 0.425000, ax = 50.00, ay = 50.00;
d20      : drift, l = 0.150000, ax = 50.00, ay = 50.00;
d22      : drift, l = 0.425000, ax = 50.00, ay = 50.00;
lo       : drift, l = 0.050000, ax = 50.00, ay = 50.00;
q1p2k3   : quadrupole, l = 0.200000, k = 2.962400, ax = 50.00, ay = 50.00;
q1p1k1   : quadrupole, l = 0.200000, k = 2.474600, ax = 50.00, ay = 50.00;
q1p2k1   : quadrupole, l = 0.200000, k = 2.474600, ax = 50.00, ay = 50.00;
q1p2l4   : quadrupole, l = 0.200000, k = 2.962400, ax = 50.00, ay = 50.00;
q1p1l2   : quadrupole, l = 0.200000, k = 2.962400, ax = 50.00, ay = 50.00;
q1p1l4   : quadrupole, l = 0.200000, k = 2.962400, ax = 50.00, ay = 50.00;
q2p1l4   : quadrupole, l = 0.200000, k = -4.174990, ax = 50.00, ay = 50.00;
q2p1k1   : quadrupole, l = 0.200000, k = -4.506670, ax = 50.00, ay = 50.00;
q3p1k1   : quadrupole, l = 0.200000, k = 5.005740, ax = 50.00, ay = 50.00;
q2p1l2   : quadrupole, l = 0.200000, k = -4.174990, ax = 50.00, ay = 50.00;
q3p1l2   : quadrupole, l = 0.200000, k = 5.191920, ax = 50.00, ay = 50.00;
q2p2l4   : quadrupole, l = 0.200000, k = -4.174990, ax = 50.00, ay = 50.00;
q3p1l4   : quadrupole, l = 0.200000, k = 5.191920, ax = 50.00, ay = 50.00;
q3p2k1   : quadrupole, l = 0.200000, k = 5.005740, ax = 50.00, ay = 50.00;
q3p2l4   : quadrupole, l = 0.200000, k = 5.191920, ax = 50.00, ay = 50.00;
q2p2l2   : quadrupole, l = 0.200000, k = -4.174990, ax = 50.00, ay = 50.00;
q3p2l2   : quadrupole, l = 0.200000, k = 5.191920, ax = 50.00, ay = 50.00;
q2p2k1   : quadrupole, l = 0.200000, k = -4.506670, ax = 50.00, ay = 50.00;
q1p2l2   : quadrupole, l = 0.200000, k = 2.962400, ax = 50.00, ay = 50.00;
q1p1k3   : quadrupole, l = 0.200000, k = 2.962400, ax = 50.00, ay = 50.00;
q2p1k3   : quadrupole, l = 0.200000, k = -4.506670, ax = 50.00, ay = 50.00;
q3p1k3   : quadrupole, l = 0.200000, k = 5.005740, ax = 50.00, ay = 50.00;
q2p2k3   : quadrupole, l = 0.200000, k = -4.506670, ax = 50.00, ay = 50.00;
q3p2k3   : quadrupole, l = 0.200000, k = 5.005740, ax = 50.00, ay = 50.00;
bend1    : bending, l = 1.200000, t = 45.000000, k = 0.000000, 
           t1 = 22.500000, t2 = 22.500000, gap = 50.00, 
           k1in = 0.5000, k1ex = 0.5000, k2in = 0.0000, 
           k2ex = 0.0000, ax = 50.00, ay = 50.00;
bend2    : bending, l = 1.200000, t = 45.000000, k = 0.000000, 
           t1 = 22.500000, t2 = 22.500000, gap = 50.00, 
           k1in = 0.5000, k1ex = 0.5000, k2in = 0.0000, 
           k2ex = 0.0000, ax = 50.00, ay = 50.00;
bend3    : bending, l = 1.200000, t = 45.000000, k = 0.000000, 
           t1 = 22.500000, t2 = 22.500000, gap = 50.00, 
           k1in = 0.5000, k1ex = 0.5000, k2in = 0.0000, 
           k2ex = 0.0000, ax = 50.00, ay = 50.00;
bend4    : bending, l = 1.200000, t = 45.000000, k = 0.000000, 
           t1 = 22.500000, t2 = 22.500000, gap = 50.00, 
           k1in = 0.5000, k1ex = 0.5000, k2in = 0.0000, 
           k2ex = 0.0000, ax = 50.00, ay = 50.00;
bend5    : bending, l = 1.200000, t = 45.000000, k = 0.000000, 
           t1 = 22.500000, t2 = 22.500000, gap = 50.00, 
           k1in = 0.5000, k1ex = 0.5000, k2in = 0.0000, 
           k2ex = 0.0000, ax = 50.00, ay = 50.00;
bend6    : bending, l = 1.200000, t = 45.000000, k = 0.000000, 
           t1 = 22.500000, t2 = 22.500000, gap = 50.00, 
           k1in = 0.5000, k1ex = 0.5000, k2in = 0.0000, 
           k2ex = 0.0000, ax = 50.00, ay = 50.00;
bend7    : bending, l = 1.200000, t = 45.000000, k = 0.000000, 
           t1 = 22.500000, t2 = 22.500000, gap = 50.00, 
           k1in = 0.5000, k1ex = 0.5000, k2in = 0.0000, 
           k2ex = 0.0000, ax = 50.00, ay = 50.00;
bend8    : bending, l = 1.200000, t = 45.000000, k = 0.000000, 
           t1 = 22.500000, t2 = 22.500000, gap = 50.00, 
           k1in = 0.5000, k1ex = 0.5000, k2in = 0.0000, 
           k2ex = 0.0000, ax = 50.00, ay = 50.00;
s2p2k1   : sextupole, l = 0.100000, k = -1.000000, n =1, 
           ax = 50.00, ay = 50.00;
s1p2k3   : sextupole, l = 0.100000, k = 1.000000, n =1, 
           ax = 50.00, ay = 50.00;
s1p2l4   : sextupole, l = 0.100000, k = 1.000000, n =1, 
           ax = 50.00, ay = 50.00;
s1p2k1   : sextupole, l = 0.100000, k = 1.000000, n =1, 
           ax = 50.00, ay = 50.00;
s1p1l4   : sextupole, l = 0.100000, k = 1.000000, n =1, 
           ax = 50.00, ay = 50.00;
s2p1k1   : sextupole, l = 0.100000, k = -1.000000, n =1, 
           ax = 50.00, ay = 50.00;
s1p1k1   : sextupole, l = 0.100000, k = 1.000000, n =1, 
           ax = 50.00, ay = 50.00;
s1p1l2   : sextupole, l = 0.100000, k = 1.000000, n =1, 
           ax = 50.00, ay = 50.00;
s2p1l4   : sextupole, l = 0.100000, k = -1.000000, n =1, 
           ax = 50.00, ay = 50.00;
s3p1k1   : sextupole, l = 0.100000, k = 1.000000, n =1, 
           ax = 50.00, ay = 50.00;
s2p1l2   : sextupole, l = 0.100000, k = -1.000000, n =1, 
           ax = 50.00, ay = 50.00;
s2p2k3   : sextupole, l = 0.100000, k = -1.000000, n =1, 
           ax = 50.00, ay = 50.00;
s3p1l4   : sextupole, l = 0.100000, k = 1.000000, n =1, 
           ax = 50.00, ay = 50.00;
s3p2l2   : sextupole, l = 0.100000, k = 1.000000, n =1, 
           ax = 50.00, ay = 50.00;
s3p1l2   : sextupole, l = 0.100000, k = 1.000000, n =1, 
           ax = 50.00, ay = 50.00;
s3p2l4   : sextupole, l = 0.100000, k = 1.000000, n =1, 
           ax = 50.00, ay = 50.00;
s3p2k1   : sextupole, l = 0.100000, k = 1.000000, n =1, 
           ax = 50.00, ay = 50.00;
s2p2l4   : sextupole, l = 0.100000, k = -1.000000, n =1, 
           ax = 50.00, ay = 50.00;
s1m2l2   : sextupole, l = 0.100000, k = 1.000000, n =1, 
           ax = 50.00, ay = 50.00;
s1p1k3   : sextupole, l = 0.100000, k = 1.000000, n =1, 
           ax = 50.00, ay = 50.00;
s2p2l2   : sextupole, l = 0.100000, k = -1.000000, n =1, 
           ax = 50.00, ay = 50.00;
s2p1k3   : sextupole, l = 0.100000, k = -1.000000, n =1, 
           ax = 50.00, ay = 50.00;
s3p1k3   : sextupole, l = 0.100000, k = 1.000000, n =1, 
           ax = 50.00, ay = 50.00;
s3p2k3   : sextupole, l = 0.100000, k = 1.000000, n =1, 
           ax = 50.00, ay = 50.00;
o1       : multipole, n = 4, k = 0.000, ax = 15.00, ay = 15.00;
o2       : multipole, n = 4, k = 0.000, ax = 15.00, ay = 15.00;
o3       : multipole, n = 4, k = 0.000, ax = 15.00, ay = 15.00;
o4       : multipole, n = 4, k = 0.000, ax = 15.00, ay = 15.00;
 
{----- table of segments ----------------------------------------------------}
 
oct1 : lo, o1, lo;
oct2 : lo, o2, lo;
oct3 : lo, o3, lo;
oct4 : lo, o4, lo;
ring : d1, s3p2k1, d101, q3p2k1, d2, q2p2k1, d3, bend1, d4, s2p2k1,
       d401, s1p2k1, d402, q1p2k1, d51, oct1, d52, q1p1l2, d6, s1p1l2, d601,
       s2p1l2, d602, bend2, d7, q2p1l2, d8, q3p1l2, d9, s3p1l2, d901, d10,
       s3p2l2, d1001, q3p2l2, d11, q2p2l2, d12, bend1, d13, s2p2l2, d1301,
       s1m2l2, d1302, q1p2l2, d141, oct2, d142, q1p1k3, d15, s1p1k3, d1501,
       s2p1k3, d1502, bend1, d16, q2p1k3, d17, q3p1k3, d18, s3p1k3, d1801,
       d19, s3p2k3, d1901, q3p2k3, d20, q2p2k3, d21, bend6, d22, s2p2k3, d2201,
       s1p2k3, d2202, q1p2k3, d231, oct3, d232, q1p1l4, d24, s1p1l4, d2401,
       s2p1l4, d2402, bend6, d25, q2p1l4, d26, q3p1l4, d27, s3p1l4, d2701,
       d28, s3p2l4, d2801, q3p2l4, d29, q2p2l4, d30, bend7, d31, s2p2l4, d3101,
       s1p2l4, d3102, q1p2l4, d321, oct4, d322, q1p1k1, d33, s1p1k1, d3301,
       s2p1k1, d3302, bend8, d34, q2p1k1, d35, q3p1k1, d36, s3p1k1, d3601;
 
{ f:\work\opa\examples\mlssu.opa }
