/* Tracy-2

   J. Bengtsson, CBP, LBL      1990 - 1994   Pascal version
                 SLS, PSI      1995 - 1997
   M. Boege      SLS, PSI      1998          C translation
   L. Nadolski   SOLEIL        2002          Link to NAFF, Radia field maps
   J. Bengtsson  NSLS-II, BNL  2004 -        

*/

//#include "complexeheader_naff.h"

double  pi = M_PI;

/****************************************************************************/
/* void Trac_Simple(double x, double px, double y, double py, double dp, long nmax,
                 double Tx[][NTURN], bool *status)

   Purpose:
       Single particle tracking around the closed orbit for NTURN turns
       The 6D phase trajectory is saved in a array

   Input:
       x, px, y, py 4 transverses coordinates
       dp           energy offset
       nmax         number of turns
       pos          starting position for tracking
       aperture     global physical aperture

   Output:
      lastn         last n (should be nmax if  not lost)
      lastpos       last position in the ring
      Tx            6xNTURN matrix of phase trajectory

   Return:
       none

   Global variables:
       NTURN number of turn for tracking
       Lattice.param

   Specific functions:
       Cell_Pass

   Comments:
       useful for connection with NAFF
       19/01/03 tracking around the closed orbit
       20/03/04 closed orbit searched w/ findcod (Newton-Raphson numerical method)
****************************************************************************/
void Trac_Simple(double x, double px, double y, double py, double dp, 
                 double ctau, long nmax, double Tx[][NTURN], bool *status2)
{
  bool             lostF = false; /* Lost particle Flag */
  ss_vect<double>  x1;  /* Tracking coordinates */
  long             lastpos = Lattice.param.Cell_nLoc;
  long             lastn = 0;
  Vector2          aperture = {1.0, 1.0};

  bool cod = false;

  *status2 = true; /* stable */

  if (cod) {
    // printf("Trac_Simple: relative (to COD)\n");
    /* Get closed orbit */
    //~ getcod(dp, lastpos);
    findcod(dp);
  
    if (trace && status.codflag) 
      printf("dp= % .5e %% xcod= % .5e mm zcod= % .5e mm \n",
             dp*1e2, Lattice.param.CODvect[0]*1e3, Lattice.param.CODvect[2]*1e3);

    /* Tracking coordinates around the closed orbit */
    x1[0] =  x + Lattice.param.CODvect[0]; x1[1] = px + Lattice.param.CODvect[1];
    x1[2] =  y + Lattice.param.CODvect[2]; x1[3] = py + Lattice.param.CODvect[3];
  } else {
    // printf("Trac_Simple: absolute\n");
    x1[0] = x; x1[1] = px; x1[2] = y; x1[3] = py;
  }

  x1[4] = dp; x1[5] = ctau;

  Tx[0][lastn] = x1[0]; Tx[1][lastn] = x1[1];
  Tx[2][lastn] = x1[2]; Tx[3][lastn] = x1[3];
  Tx[4][lastn] = x1[4]; Tx[5][lastn] = x1[5];
  lastn++;

  do
  { /* tracking through the ring */
    if ((lastpos == Lattice.param.Cell_nLoc) &&
        (fabs(x1[0]) < aperture[0]) && (fabs(x1[2]) < aperture[1]) && status.codflag)
    {
      Cell_Pass(0, Lattice.param.Cell_nLoc, x1, lastpos);
      Tx[0][lastn] = x1[0]; Tx[1][lastn] = x1[1];
      Tx[2][lastn] = x1[2]; Tx[3][lastn] = x1[3];
      Tx[4][lastn] = x1[4]; Tx[5][lastn] = x1[5];
    }
    else
    {
      // printf("Trac_Simple: Particle lost \n");
      // printf("%6ld plane:"
      // 	     " %1d %+10.5g %+10.5g %+10.5g %+10.5g %+10.5g %+10.5g \n", 
      // 	     lastn, status.lossplane,
      // 	     x1[0], x1[1], x1[2], x1[3], x1[4], x1[5]);
      lostF = true;
      *status2 = false;
    }
    lastn++;
  } while ((lastn < nmax) && (lastpos == Lattice.param.Cell_nLoc)
	   && (lostF == false));

  if (lastpos != Lattice.param.Cell_nLoc)
  { /* Particle lost: Error message section */
    *status2 = false;
    // printf("Trac_Simple: Particle lost \n");
    // printf("turn=%5ld plane= %1d"
    // 	   " %+10.5g %+10.5g %+10.5g %+10.5g %+10.5g %+10.5g \n",
    // 	   lastn-1, status.lossplane,
    // 	   x1[0], x1[1], x1[2], x1[3], x1[4], x1[5]);
  }
}


/****************************************************************************/
/* void Get_NAFF(int nterm, long ndata, double T[DIM][NTURN],
              double *fx, double *fz, int nb_freq[2])

   Purpose:
       Compute quasiperiodic approximation of a phase space trajectory
       using NAFF Algorithm ((c) Laskar, IMCCE)

   Input:
       nterm number of frequencies to look for
             if not multiple of 6, truncated to lower value
       ndata size of the data to analyse
       T     6D vector to analyse

   Output:
       fx frequencies found in the H-plane
       fz frequencies found in the V-plane
       nb_freq number of frequencies found out in each plane

   Return:
       none

   Global variables:
       g_NAFVariable  see modnaff.c
       M_2_PI defined in math.h
       trace ON or TRUE  for debugging

   Specific functions:
       naf_initnaf, naf_cleannaf
       naf_mftnaf

   Comments:
       none

****************************************************************************/
/* Frequency Map Analysis */
/* Analyse en Frequence */
void Get_NAFF(int nterm, long ndata, double Tab[DIM][NTURN],
              double *fx, double *fz, int nb_freq[2])
{
  /* Test whether ndata is divisible by 6 -- for NAFF -- */
  /* Otherwise truncate ndata to lower value */
  long r = 0; /* remainder of the euclidian division of ndata by 6 */
  int i;

  if ((r = ndata % 6) != 0) {
    printf("Get_NAFF: Warning ndata = %ld, \n", ndata);
    ndata -= r;
    printf("New value for NAFF ndata = %ld \n", ndata);
  }

  g_NAFVariable.DTOUR      = M_2_PI;    /* size of one "cadran" */
  g_NAFVariable.XH         = M_2_PI;    /* step */
  g_NAFVariable.T0         = 0.0;       /* time t0 */
  g_NAFVariable.NTERM      = nterm;     /* max term to find */
  g_NAFVariable.KTABS      = ndata;     /* number of data: must be a multiple of 6 */
  g_NAFVariable.m_pListFen = NULL;      /* no window */
  g_NAFVariable.TFS        = NULL;      /* will contain frequency */
  g_NAFVariable.ZAMP       = NULL;      /* will contain amplitude */
  g_NAFVariable.ZTABS      = NULL;      /* will contain data to analyze */

  /****************************************************/
  /*               internal use in naf                */
  g_NAFVariable.NERROR            = 0;
  g_NAFVariable.ICPLX             = 1;
  g_NAFVariable.IPRT              = -1;     /* 1 for diagnostics */
  g_NAFVariable.NFPRT             = stdout; /* NULL   */
  g_NAFVariable.NFS               = 0;
  g_NAFVariable.IW                = 1;
  g_NAFVariable.ISEC              = 1;
  g_NAFVariable.EPSM              = 0;
  g_NAFVariable.UNIANG            = 0;
  g_NAFVariable.FREFON            = 0;
  g_NAFVariable.ZALP              = NULL;
  g_NAFVariable.m_iNbLineToIgnore = 1;      /* unused */
  g_NAFVariable.m_dneps           = 1.e10;
  g_NAFVariable.m_bFSTAB          = FALSE;  /* unused */
  /*             end of interl use in naf             */
  /****************************************************/

  /* NAFF initialization */
  naf_initnaf();

  /**********************/
  /* Analyse in H-plane */
  /**********************/

  /* fills up complexe vector for NAFF analysis */
  for(i = 0; i < ndata; i++) {
    g_NAFVariable.ZTABS[i].reel = Tab[0][i]; /* x  */
    g_NAFVariable.ZTABS[i].imag = Tab[1][i]; /* xp */
  }

  /* Get out the mean value */
  naf_smoy(g_NAFVariable.ZTABS);

  naf_prtabs(g_NAFVariable.KTABS,g_NAFVariable.ZTABS, 20);
//  trace=on;
  naf_mftnaf(nterm,fabs(g_NAFVariable.FREFON)/g_NAFVariable.m_dneps);

  /* fill up H-frequency vector */
  for (i = 1; i <= g_NAFVariable.NFS; i++)  {
    fx[i-1] = g_NAFVariable.TFS[i];
  }

  nb_freq[0] = g_NAFVariable.NFS; /* nb of frequencies found out by NAFF */

  if (trace)   /* print out results */
  {
    printf("(x,x') phase space: NFS=%d\n",g_NAFVariable.NFS);
    for (i = 1; i <= g_NAFVariable.NFS; i++) {
      printf("AMPL=%15.8E+i*%15.8E abs(AMPL)=%15.8E arg(AMPL)=%15.8E FREQ=%15.8E\n",
              g_NAFVariable.ZAMP[i].reel,g_NAFVariable.ZAMP[i].imag,
              i_compl_module(g_NAFVariable.ZAMP[i]),
              i_compl_angle(g_NAFVariable.ZAMP[i]),
              g_NAFVariable.TFS[i]);
    }
  }

  /**********************/
  /* Analyse in V-plane */
  /**********************/

  /* fill up complexe vector for NAFF analysis */
  for (i = 0; i < ndata; i++) {
    g_NAFVariable.ZTABS[i].reel = Tab[2][i];  /* z */
    g_NAFVariable.ZTABS[i].imag = Tab[3][i];  /*zp */
  }

  naf_mftnaf(nterm,fabs(g_NAFVariable.FREFON)/g_NAFVariable.m_dneps);

  /* fills up V-frequency vector */
  for (i = 1; i <= g_NAFVariable.NFS; i++) {
    fz[i-1] =  g_NAFVariable.TFS[i];
  }

  nb_freq[1] = g_NAFVariable.NFS; /* nb of frequencies found out by NAFF */

  if (trace)    /* print out results */
  {
    printf("(z,z') phase space: NFS=%d\n",g_NAFVariable.NFS);
    for (i = 1; i <= g_NAFVariable.NFS; i++) {
      printf("AMPL=%15.8E+i*%15.8E abs(AMPL)=%15.8E arg(AMPL)=%15.8E FREQ=%15.8E\n",
              g_NAFVariable.ZAMP[i].reel,g_NAFVariable.ZAMP[i].imag,
              i_compl_module(g_NAFVariable.ZAMP[i]),
              i_compl_angle(g_NAFVariable.ZAMP[i]),
              g_NAFVariable.TFS[i]);
    }
  }

  /* free out memory used by NAFF */
  naf_cleannaf();
}

/****************************************************************************/
/* void Get_Tabshift(double Tab[DIM][NTURN], double Tab0[DIM][NTURN],
   long nbturn, long nshift)

   Purpose:   used by fmap
       Store in Tab0 values of Tab shifted by nshift turn

   Input:
       Tab    tracking tabular
       nshift nb of turns to shift
       nbturn nb of turns

   Output:
       Tab0  tracking tabular

   Return:
       none

   Global variables:
       none

   Specific functions:
       none

   Comments:
       

****************************************************************************/
void Get_Tabshift(double Tab[DIM][NTURN], double Tab0[DIM][NTURN], long nbturn, long nshift)
{
  long k = 0L, k1 = 0L;

  for (k = 0L; k < nbturn ; k++){
    k1 = k + nshift;
    Tab0[0][k] = Tab[0][k1];
    Tab0[1][k] = Tab[1][k1];
    Tab0[2][k] = Tab[2][k1];
    Tab0[3][k] = Tab[3][k1];
  }

}

/****************************************************************************/
/* void Get_freq(double *fx, double *fz, double *nux, double *nuz)

   Purpose:   used by fmap
       Looks for tunes from frequency vectors output from NAFF

   Input:
       fx vector of horizontal frequencies
       fz vector of verticcal frequencies

   Output:
       nux H tune
       nuz V tune

   Return:
       none

   Global variables:
       none

   Specific functions:
       none

   Comments:


****************************************************************************/
void Get_freq(double *fx, double *fz, double *nux, double *nuz)
{
  const double eps0 = 1e-4;
  const double eps1 = 1e-6;
  
  // case of nux
  if (fabs(fx[0]) < eps0){
    *nux = fx[1];
  }
  else {
    *nux = fx[0];
  }

  // case of nuz
  if (fabs(fz[0]) < eps0) {
     if (fabs(fabs(fz[1]) - fabs(*nux)) < eps1) {
       if (fabs(fabs(fz[2]) - fabs(*nux)) < eps1) {
          *nuz = fz[3];
       }
       else *nuz = fz[2];
     }
     else *nuz = fz[1];
  }
  else{
    if (fabs(fabs(fz[0]) - fabs(*nux)) < eps1) {
      if (fabs(fabs(fz[1]) - fabs(*nux)) < eps1) {
         *nuz = fz[2];
      }
      else *nuz = fz[1];
    }
    else *nuz = fz[0];
  }
  *nuz = fabs(*nuz);  *nux = fabs(*nux);
}
