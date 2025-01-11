function dls_bts
% 
% DIAMOND BTS
% for the AT; RB 18/04/2005
%
%
global FAMLIST THERING GLOBVAL

GLOBVAL.E0 = 0.1e9;
GLOBVAL.LatticeFile = 'BTS';
FAMLIST = cell(0);

disp(' ');
disp('** Loading DIAMOND BTS **');

AP  =  aperture('AP', [-0.05, 0.05, -0.05, 0.05],'AperturePass');
% drifts
BPM     = drift('BBPM', .05E+00,'DriftPass');
BPMSEP  = drift('BPMSEP', .27E+00,'DriftPass');
BPMSEP1  = drift('BPMSEP1', .2035E+00,'DriftPass');
BSEPTTOQA   = drift('BSEPTTOQA', .53E+00,'DriftPass');
BSEPTTOQBa   = drift('BSEPTTOQBa',0.530E+00,'DriftPass');
BSEPTTOQBb   = drift('BSEPTTOQBb',1.575E+00,'DriftPass');
BR03CDIOTR03 = marker('BR03CDIOTR03', 'IdentityPass');
BSEPTTOQB = [BSEPTTOQBa BR03CDIOTR03 BSEPTTOQBb];
BSEPTTOQC   = drift('BSEPTTOQC', .1755E+00,'DriftPass');
CORRSEP = drift('CORRSEP', .3E+00,'DriftPass');
DFT10   = drift('DFT10', .01E+00,'DriftPass');
DFT120  = drift('DFT120', .12E+00,'DriftPass');
DFT1530 = drift('DFT1520', 1.53E+00,'DriftPass');
DFT180  = drift('DFT180', .18E+00,'DriftPass');
DFT200  = drift('DFT200', .2E+00,'DriftPass');
DFT2170 = drift('DFT2170', 2.17E+00,'DriftPass');
DFT250  = drift('DFT250', .25E+00,'DriftPass');
DFT300  = drift('DFT300', .3E+00,'DriftPass');
DFT520  = drift('DFT520', .52E+00,'DriftPass');
DFT580  = drift('DFT580', .58E+00,'DriftPass');
DFT60   = drift('DFT60', .06E+00,'DriftPass');
MMS   = drift('MMS', .0E+00,'DriftPass');
MPS   = drift('MPS', .0E+00,'DriftPass');
PSTRAIGHT1  = drift('PSTRAIGHT1',0.5,'DriftPass');
PSTRAIGHT10 = drift('PSTRAIGHT10',1.14,'DriftPass');
PSTRAIGHT11 = drift('PSTRAIGHT11',0.5,'DriftPass');
PSTRAIGHT12 = drift('PSTRAIGHT12',1.404,'DriftPass');
PSTRAIGHT13 = drift('PSTRAIGHT13',2.89,'DriftPass');
PSTRAIGHT14 = drift('PSTRAIGHT14',0.5,'DriftPass');
PSTRAIGHT15A= drift('PSTRAIGHT15A',0.671,'DriftPass');
PSTRAIGHT15B= drift('PSTRAIGHT15B',4.488,'DriftPass');
PSTRAIGHT15C= drift('PSTRAIGHT15C',0.285,'DriftPass');
PSTRAIGHT2  = drift('PSTRAIGHT2',1.864,'DriftPass');
PSTRAIGHT3  = drift('PSTRAIGHT3',0.5,'DriftPass');
PSTRAIGHT4Aa = drift('PSTRAIGHT4Aa',1.804,'DriftPass');
PSTRAIGHT4Ab = drift('PSTRAIGHT4Ab',2.434,'DriftPass');
BSDIOTR01 = marker('BSDIOTR01', 'IdentityPass');
PSTRAIGHT4A = [PSTRAIGHT4Aa BSDIOTR01 PSTRAIGHT4Ab];
PSTRAIGHT4B = drift('PSTRAIGHT4B',0.176,'DriftPass');
PSTRAIGHT5  = drift('PSTRAIGHT5',0.5,'DriftPass');
PSTRAIGHT6A = drift('PSTRAIGHT6A',6.2185,'DriftPass');
PSTRAIGHT6B = drift('PSTRAIGHT6B',0.4655,'DriftPass');
PSTRAIGHT7  = drift('PSTRAIGHT7',1.28,'DriftPass');
PSTRAIGHT8  = drift('PSTRAIGHT8',0.5,'DriftPass');
PSTRAIGHT9A = drift('PSTRAIGHT9A',5.539,'DriftPass');
PSTRAIGHT9B = drift('PSTRAIGHT9B',0.465,'DriftPass');
cdrift = drift('cdrift',0.206,'DriftPass');

% quads
BDQUD =  quadrupole('BQD', 0.34, -1.016147,'QuadLinearPass');
BFQUD =  quadrupole('BQF', 0.34,  1.402775,'QuadLinearPass');
DQUD1 =  quadrupole('QUAD', 0.4,  -1.3618983,'QuadLinearPass');
DQUD2 =  quadrupole('QUAD', 0.4,  -0.46133062,'QuadLinearPass');
DQUD3 =  quadrupole('QUAD', 0.4,  -0.54609978,'QuadLinearPass');
DQUD4 =  quadrupole('QUAD', 0.4,  -0.89580103,'QuadLinearPass');
DQUD5 =  quadrupole('QUAD', 0.4,  -1.0407761,'QuadLinearPass');
DQUD6 =  quadrupole('QUAD', 0.4,  -0.88235615,'QuadLinearPass');
FQUD1 =  quadrupole('QUAD', 0.4,   1.4195575,'QuadLinearPass');
FQUD2 =  quadrupole('QUAD', 0.4,   0.63588467,'QuadLinearPass');
FQUD3 =  quadrupole('QUAD', 0.4,   0.34859322,'QuadLinearPass');
FQUD4 =  quadrupole('QUAD', 0.4,   1.0404348,'QuadLinearPass');
FQUD5 =  quadrupole('QUAD', 0.4,   1.1194183,'QuadLinearPass');
FQUD6 =  quadrupole('QUAD', 0.4,   1.0832694,'QuadLinearPass');

% correctors
HBTSC1 = corrector('HSTR',1e-6,[ 0 0 ],'CorrectorPass');
VBTSC1 = corrector('VSTR',1e-6,[ 0 0 ],'CorrectorPass');
BTSC1 = [HBTSC1 cdrift VBTSC1];
HBTSC2 = corrector('HSTR',1e-6,[ 0 0 ],'CorrectorPass');
VBTSC2 = corrector('VSTR',1e-6,[ 0 0 ],'CorrectorPass');
BTSC2 = [HBTSC1 cdrift VBTSC1];
HBTSC3 = corrector('HSTR',1e-6,[ 0 0 ],'CorrectorPass');
VBTSC3 = corrector('VSTR',1e-6,[ 0 0 ],'CorrectorPass');
BTSC3 = [HBTSC1 cdrift VBTSC1];
HBTSC4 = corrector('HSTR',1e-6,[ 0 0 ],'CorrectorPass');
VBTSC4 = corrector('VSTR',1e-6,[ 0 0 ],'CorrectorPass');
BTSC4 = [HBTSC1 cdrift VBTSC1];
HBTSC5 = corrector('HSTR',1e-6,[ 0 0 ],'CorrectorPass');
VBTSC5 = corrector('VSTR',1e-6,[ 0 0 ],'CorrectorPass');
BTSC5 = [HBTSC1 cdrift VBTSC1];
HBTSC6 = corrector('HSTR',1e-6,[ 0 0 ],'CorrectorPass');
VBTSC6 = corrector('VSTR',1e-6,[ 0 0 ],'CorrectorPass');
BTSC6 = [HBTSC1 cdrift VBTSC1];
HBTSC7 = corrector('HSTR',1e-6,[ 0 0 ],'CorrectorPass');
VBTSC7 = corrector('VSTR',1e-6,[ 0 0 ],'CorrectorPass');
BTSC7 = [HBTSC1 cdrift VBTSC1];

hk = corrector('HSTRK',1e-6,[ 0 0 ],'CorrectorPass');
hkdrift = drift('hkdrift', .08E+00,'DriftPass');
HKICKER = [hkdrift hk hkdrift];
vk = corrector('VSTRK',1e-6,[ 0 0 ],'CorrectorPass');
vkdrift = drift('vkdrift', .08E+00,'DriftPass');
VKICKER = [vkdrift vk vkdrift];

% Bending
BBANGLE=0.174533;
BEND1  =   sbend('BBTS', 2.16, ...
            BBANGLE, BBANGLE/2, BBANGLE/2, 0,'BendLinearPass');
        
BANGLEFK=-0.003; % AP-BST-REP-0049: 3 mrad
FK      =  sbend('FK', 1.0, BANGLEFK, BANGLEFK/2, BANGLEFK/2, 0, 'BendLinearPass');
% FK      =  drift('FK', 1.0, 'DriftPass');

% BAPS = -0.00515; % AP-BST-REP-0049: 4.66 mrad
BAPS = -0.00466; % AP-BST-REP-0049: 4.66 mrad
PS  =  sbend('PS', 0.36, BAPS, BAPS/2, BAPS/2, 0, 'BendLinearPass');
        
% BANGLEMS=-0.108366; % AP-BST-REP-0049: 108.03 mrad
BANGLEMS=-0.10803; % AP-BST-REP-0049: 108.03 mrad
MS      =  sbend('MS', 1.2, BANGLEMS, BANGLEMS/2, BANGLEMS/2,0, 'BendLinearPass');

% the booster quads act to bend the reference trajectory
% see booster_extraction.m and AP-BST-REP-0038
thDQ1S = -0.001541/2;
thFQ1S =  0.009597/2;
thDQ2S = -0.005793/2;
% thDQ1S = -0.001541/2;
% thFQ1S =  0.009597/2;
% thDQ2S = -0.006311/2;
DQ1S      =  sbend('DQ1S', 0, thDQ1S, 0, 0 ,0, 'BendLinearPass');
FQ1S      =  sbend('FQ1S', 0, thFQ1S, 0, 0 ,0, 'BendLinearPass');
DQ2S      =  sbend('DQ2S', 0, thDQ2S, 0, 0 ,0, 'BendLinearPass');

BABEND=0.314159;
PBEND1  =  sbend('BB', 2.17, ...
            BABEND, BABEND/2, BABEND/2,0, 'BendLinearPass');
PBEND2  =  sbend('BB', 2.17, ...
            BABEND, BABEND/2, BABEND/2,0, 'BendLinearPass');
PBEND3  =  sbend('BB', 2.17, ...
            BABEND, BABEND/2, BABEND/2,0, 'BendLinearPass');

BASEPTUM = 0.15;
SRSEPTUM  =  sbend('SRSEPTUM', 1.9, ...
            BASEPTUM, BASEPTUM/2, BASEPTUM/2, 0, 'BendLinearPass');

BTSBPM1 = marker('BPM', 'IdentityPass');
BTSBPM2 = marker('BPM', 'IdentityPass');
BTSBPM3 = marker('BPM', 'IdentityPass');
BTSBPM4 = marker('BPM', 'IdentityPass');
BTSBPM5 = marker('BPM', 'IdentityPass');
BTSBPM6 = marker('BPM', 'IdentityPass');
BTSBPM7 = marker('BPM', 'IdentityPass');
ECOLL = marker('ECOLL', 'IdentityPass');
EP = marker('EP', 'IdentityPass');
IP      =  marker('IP', 'IdentityPass');
HCOLL1 = marker('HCOLL1', 'IdentityPass');
HCOLL2 = marker('HCOLL2', 'IdentityPass');
VCOLL1 = marker('VCOLL1', 'IdentityPass');
VCOLL2 = marker('VCOLL2', 'IdentityPass');

% Begin Lattice
BTS = [BFQUD, DFT1530, FK, DFT250, DFT200, VKICKER, DFT10, BPM,...
DFT60, DQ1S, BDQUD, DQ1S, DFT520, BEND1, DFT300, HKICKER, DFT120, ...
FQ1S, BFQUD, FQ1S, DFT250, MPS, PS, DFT2170, DFT200, VKICKER, DFT120, ...
DQ2S, BDQUD, DQ2S, DFT580, MMS, MS, DFT180, EP,...
BSEPTTOQA, BTSC1, BSEPTTOQB, BTSBPM1, BSEPTTOQC, HCOLL1, BPMSEP1,...
FQUD1, PSTRAIGHT1, DQUD1, CORRSEP, BTSC2, PSTRAIGHT2, BTSBPM2, BPMSEP,...
DQUD2, PSTRAIGHT3, FQUD2, CORRSEP, BTSC3, PSTRAIGHT4A, VCOLL1,...
PSTRAIGHT4B, BTSBPM3, BPMSEP, DQUD3, PSTRAIGHT5, FQUD3, CORRSEP,...
BTSC4, PSTRAIGHT6A, HCOLL2, PSTRAIGHT6B, PBEND1, PSTRAIGHT7, BTSBPM4,...
BPMSEP, DQUD4, PSTRAIGHT8, FQUD4, CORRSEP, BTSC5, PSTRAIGHT9A, ECOLL,...
PSTRAIGHT9B, PBEND2, PSTRAIGHT10, BTSBPM5, BPMSEP, FQUD5, PSTRAIGHT11,...
DQUD5, CORRSEP, BTSC6, PSTRAIGHT12, PBEND3, PSTRAIGHT13, BTSBPM6,...
BPMSEP, FQUD6, PSTRAIGHT14, DQUD6, CORRSEP, BTSC7, PSTRAIGHT15A,...
VCOLL2, PSTRAIGHT15B, BTSBPM7, PSTRAIGHT15C, SRSEPTUM, IP];
            
buildlat(BTS);
evalin('caller','global THERING FAMLIST GLOBVAL');
disp('** Done **');
