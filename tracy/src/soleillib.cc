/* Tracy-2

   J. Bengtsson, CBP, LBL      1990 - 1994   Pascal version
                 SLS, PSI      1995 - 1997
   M. Boege      SLS, PSI      1998          C translation
   L. Nadolski   SOLEIL        2002          Link to NAFF, Radia field maps
   J. Bengtsson  NSLS-II, BNL  2004 -        

*/


/****************************************************************************/
/* void Get_Disp_dp(void)

   Purpose:
       Get dispersion w/ energy offset

   Input:
       none

   Output:
       none

   Return:
       none

   Global variables:
       trace

   specific functions:
       getcod, Ring_GetTwiss, getelem

   Comments:
       none

****************************************************************************/
void Get_Disp_dp(void)
{
  long i;
//  long lastpos = 0;
  const char nomfic[] = "dispersion.out";
  FILE *outf;
  double dP = 0e0;
  CellType Cell;

  if (trace) fprintf(stdout,"Entering Get_Disp_dp function ...\n");

  if ((outf = fopen(nomfic, "w")) == NULL) {
    fprintf(stdout, "Get_Disp_dp: Error while opening file %s\n",nomfic);
    exit_(1);
  }

  for (i = 1; i <= 20; i++) {
    dP = -0.003 + 1e-6 + i*0.0006;
    //~ getcod(dP, lastpos);
    findcod(dP);
    Lattice.Ring_GetTwiss(true, dP);  /* Compute and get Twiss parameters */
    getelem(0, &Cell);
    fprintf(outf,"%+e %+e %+e\n", dP, Cell.BeamPos[0], Cell.Eta[0]);
  }

  fclose(outf);
}

/****************************************************************************/
/* void InducedAmplitude(long spos)

   Purpose:
      Compute the induced amplitude for a particle getting for a energy offset dP
        process similar to a Touschek scattering
        The induced maplitude is trnasported to the first element of the lattice
        by scaling the maplitude with energy dependent betafunctions        

   Input:
       spos : position where Touschek scattering occurs

   Output:
       amp_ind.out

   Return:
       none

   Global variables:
       none

   specific functions:
       none

   Comments:
       none

****************************************************************************/
void InducedAmplitude(long spos)
{
  psVector        x1;     /* tracking coordinates */
  long          i = 0L, k = 0L, imax = 50;
  FILE *        outf;
  double        dP = 0.0, dP20 = 0.0, dpmax = 0.06;
  Vector2       H = {0.0, 0.0};
  const char    nomfic[] = "amp_ind.out";
  CellType      Celldebut, Cell;
  psVector        codvector[Cell_nLocMax];

  Lattice.param.Cavity_on  = false;    /* Cavity on/off */
  Lattice.param.radiation  = false;    /* radiation on/off */

  /* Ouverture fichier moustache */
  if ((outf = fopen(nomfic, "w")) == NULL) {
    fprintf(stdout, "Erreur à l'ouverture de %s\n",nomfic);
    exit_(1);
  }

  fprintf(outf, "#    dp           xind         zind       "
                " Betax(0)     Betaz(0)       Betax         betaz"
	        "       Hx          Hz            etax        etaxp\n#\n");

 for (k = 0; k <= imax ; k++)  {
    dP = -dpmax + 2*dpmax*k/imax;
    /* Coordonnees initiales */
    x1[0] = 0.0; x1[1] = 0.0;
    x1[2] = 0.0; x1[3] = 0.0;
    x1[4] = dP ; x1[5] = 0.0;

    /* Computes closed orbit and store it in a vector */
    set_vectorcod(codvector, dP) ;
    Lattice.Ring_GetTwiss(false, dP);  /* Compute and get Twiss parameters */
    getelem(1L, &Celldebut);
    getelem(spos, &Cell);

    /* compute H at s =spos */
    dP20 = ((dP == 0) ? 1.0 : dP*dP);
    i = 0; /* Horizontal */
    H[i] = ((1.0+Cell.Alpha[i]*Cell.Alpha[i])
	    /Cell.Beta[i]*codvector[spos][0]*codvector[spos][0]
	    +2.0*Cell.Alpha[i]*codvector[spos][0]*codvector[spos][1]
	    +Cell.Beta[i]*codvector[spos][1]*codvector[spos][1])/dP20;
    i = 1; /* Vertical */
    H[i] =
      ((1.0+Cell.Alpha[i]*Cell.Alpha[i])/Cell.Beta[i]*codvector[spos][2]
       *codvector[spos][2]
       +2.0*Cell.Alpha[i]*codvector[spos][2]*codvector[spos][3]
       +Cell.Beta[i]*codvector[spos][3]*codvector[spos][3])/dP20;

    fprintf(outf, "%+10.5e %+10.5e %+10.5e %+10.5e %+10.5e %+10.5e %+10.5e "
                  "%+10.5e %+10.5e %+10.5e %+10.5e \n",
	    dP, codvector[spos][0], codvector[spos][1], Celldebut.Beta[0],
	    Celldebut.Beta[1], Cell.Beta[0], Cell.Beta[1], H[0], H[1],
	    Cell.Eta[0], Cell.Etap[0]);
  }
  fclose(outf);
}

/****************************************************************************/
/* void Hfonction(long pos, double dP,Vector2 H)

   Purpose:
     Compute the Hfunction at position pos for the energy offset dP
     H is wrong at large dp since eta and eta' are computed
       by numerical differentiation, which means that
       eta(dp) = eta0 + eta2*dp*dp + O(4) instead of
       eta(dp) = eta0 + eta1*dp + eta2*dp*dp + O(3)

     A solution is to compute eta from the closed orbit by:
       xco(dp) = eta(dp)*dp => eta(dp) = xco(dp)/dp
       WARNING: this definition is true only if the lattice
       is perfect.
       Indeed in general : xco = eta(dp)*dp + x0(defaults)

   Input:
       none

   Output:
       none

   Return:
       none

   Global variables:
       none

   specific functions:
       Lattice.Ring_GetTwiss
       getelem

   Comments:
       none

****************************************************************************/

void Hfonction(long pos, double dP,Vector2 H)
{
  CellType Cell;
  long i;

  Lattice.Ring_GetTwiss(pos, dP); /* Compute and get Twiss parameters */
  getelem(pos, &Cell);    /* Position sur l'element pos */

  i = 0; /* Horizontal */
  H[i] = (1+Cell.Alpha[i]*Cell.Alpha[i])/Cell.Beta[i]*Cell.Eta[i]*Cell.Eta[i]+
          2*Cell.Alpha[i]*Cell.Eta[i]*Cell.Etap[i]+
          Cell.Beta[i]*Cell.Etap[i]*Cell.Etap[i];
  i = 1; /* Vertical */
  H[i] = (1+Cell.Alpha[i]*Cell.Alpha[i])/Cell.Beta[i]*Cell.Eta[i]*Cell.Eta[i]+
          2*Cell.Alpha[i]*Cell.BeamPos[i]*Cell.Etap[i]+
    Cell.Beta[i]*Cell.Etap[i]*Cell.Etap[i];
}

/****************************************************************************/
/* void Hcofonction(long pos, double dP,Vector2 H)

   Purpose:
       Compute the true Hfunction defined by the chromatic closed orbit
       at position pos and for a energy offset dP

       For a givien delta
       H = gamma xcod² + 2*alpha*xcod*xcod' + beta*xcod'*xcod'
       
   Input:
       none

   Output:
       none

   Return:
       none

   Global variables:
       none

   specific functions:
       getcod
       Lattice.Ring_GetTwiss
       getelem

   Comments:
       Bug: Cell.BeamPos does not give closed orbit !!!

****************************************************************************/
void Hcofonction(long pos, double dP,Vector2 H)
{
  CellType Cell;
  long i;
  long lastpos = 1;

  //~ getcod(dP, lastpos);   /* determine closed orbit */
  findcod(dP);
  if (lastpos != Lattice.param.Cell_nLoc) printf("Ring unstable for dp=%+e @ pos=%ld\n", dP, lastpos);

  Lattice.Ring_GetTwiss(pos, dP); /* Compute and get Twiss parameters */
  getelem(pos, &Cell);    /* Position sur l'element pos */

  i = 0; /* Horizontal */
  H[i] = (1+Cell.Alpha[i]*Cell.Alpha[i])/Cell.Beta[i]*Cell.BeamPos[i]*Cell.BeamPos[i]+
          2*Cell.Alpha[i]*Cell.BeamPos[i]*Cell.BeamPos[i+1]+
          Cell.Beta[i]*Cell.BeamPos[i+1]*Cell.BeamPos[i+1];
  i = 1; /* Vertical */
  H[i] = (1+Cell.Alpha[i]*Cell.Alpha[i])/Cell.Beta[i]*Cell.BeamPos[i+1]*Cell.BeamPos[i+1]+
          2*Cell.Alpha[i]*Cell.BeamPos[i+1]*Cell.BeamPos[i+2]+
          Cell.Beta[i]*Cell.BeamPos[i+2]*Cell.BeamPos[i+2];
}

  
/****************************************************************************/
/* void SetErr(void)

   Purpose:
       Set error
       Definir une distribution aleatoire de quadripoles tournes associee a chaque
       quadripole de la machine
       Distribution gaussienne d'ecart type fac et coupe a normcut*sigma

   Input:
       none

   Output:
       none

   Return:
       none

   Global variables:
       Lattice.param
       HOMmax

   specific functions:
       setrandcut, initranf
       getelem, putelem
       Mpole_SetB

   Comments:
       Only valid if quad split into two part (cf pair variable)
       Rotation inversion to do as in BETA code
       Test if normal quad sin(theta) = 0. Do not work if tilt error

****************************************************************************/
void SetErr(void)
{
  double  fac = 0.0, normcut = 0.0;
//  long    seed = 0L;
  long    i = 0L;
  CellType Cell;
  double theta = 0.0;
  int pair = 0;

  /* coupling */
//  fac = 0.00064; //normal faux
//  fac = 0.00075; //normal faux
//  fac = 0.00155; //normal faux
//  fac = 0.00076; // mais coupe en deux
//  fac = 0.001335/2.0;

// Normal
//  fac = 0.00155/2.0;
//  setrancut(normcut=2L);
//  iniranf(seed);

    fac = 0.0014/2.0;
//  fac = 0.00115/2.0;
  setrancut(normcut=2L);
//  iniranf(seed=0L);

  for (i = 1L; i <= Lattice.param.Cell_nLoc; i++)
  {
    getelem(i, &Cell);
    if (Cell.Elem.Kind == 2L)
    {
      if (Cell.Elem.M->order == 2L && Cell.dT[1] == 0)
      {
        if ((pair%2)==0) theta = fac*normranf(); /* random error every 2 elements (quad split into 2) */
        pair++;
        Cell.Elem.M->Bpar[HOMmax-2L] =
	  -Cell.Elem.M->Bpar[HOMmax+2L]*sin(2.0*theta);
        Cell.Elem.M->Bpar[HOMmax+2L] =
	  Cell.Elem.M->Bpar[HOMmax+2L]*cos(2.0*theta);
        if (trace)
	  printf("%6s % .5e % .5e % .5e\n",
		 Cell.Elem.Name, Cell.Elem.M->Bpar[HOMmax-2L],
		 Cell.Elem.M->Bpar[HOMmax+2L], theta);

        putelem(i, &Cell);
        Mpole_SetB(Cell.Fnum, Cell.Knum, -2L);
        Mpole_SetB(Cell.Fnum, Cell.Knum, 2L);
      }
    }
  }
}


/****************************************************************************/
/* void DefineCh(void)

   Purpose:  called by read_Lattice
       Defines the vacuum chamber around the ring

   Input:
       none

   Output:
       none

   Return:
       none

   Global variables:
       Lattice.param

   specific functions:
       none

   Comments:
       This function works providing that the makers have unique names into
       the lattice

****************************************************************************/
void DefineCh(void)
{
  long i;
  long isep1 = 0, isep2 = 0, hu600 = 0,
       isdm1 = 0, isdm2 = 0, isdac1 = 0, isdac2 = 0;

/*  const long isep1 = 3, isep2 = 6, hu600=8,
             isdm1 = 69, isdm2 = 98,
             isdac1 = 131, isdac2 = 140;
*/

  /* Look for indices for defining the vaccum pipe*/
  for (i = 0; i <= Lattice.param.Cell_nLoc; i++) {
    if (Lattice.Cell[i].Elem.Kind == marker){
      if (strncmp(Lattice.Cell[i].Elem.Name,"ssep",4) == 0){
        if (trace) fprintf(stdout,"trouve %s Element numero %ld \n",
			   Lattice.Cell[i].Elem.Name,i);
        isep1 = i;
      }
      if (strncmp(Lattice.Cell[i].Elem.Name,"esep",4) == 0) {
        if (trace) fprintf(stdout,"trouve %s Element numero %ld \n",
			   Lattice.Cell[i].Elem.Name,i-1);
        isep2 = i-1;
      }
      if (strncmp(Lattice.Cell[i].Elem.Name,"ehu600",6) == 0) {
        if (trace) fprintf(stdout,"trouve %s Element numero %ld \n",
			   Lattice.Cell[i].Elem.Name,i-1);
        hu600 = i-1;
      }
      if (strncmp(Lattice.Cell[i].Elem.Name,"ssdm",4) == 0) {
         if (trace) fprintf(stdout,"trouve %s Element numero %ld \n",
			    Lattice.Cell[i].Elem.Name,i);
         isdm1 = i;
      }
      if (strncmp(Lattice.Cell[i].Elem.Name,"esdm",4) == 0) {
        if (trace) fprintf(stdout,"trouve %s Element numero %ld \n",
			   Lattice.Cell[i].Elem.Name,i);
        isdm2 = i;
      }
      if (strncmp(Lattice.Cell[i].Elem.Name,"ssdac",5) == 0) {
        if (trace) fprintf(stdout,"trouve %s Element numero %ld \n",
			   Lattice.Cell[i].Elem.Name,i);
        isdac1 = i;
      }
      if (strncmp(Lattice.Cell[i].Elem.Name,"esdac",5) == 0) {
        if (trace) fprintf(stdout,"trouve %s Element numero %ld \n",
			   Lattice.Cell[i].Elem.Name,i-1);
        isdac2 = i-1;
      }
    }
  }
  trace=0;

  /* Set the vacuum chamber */

  for (i = 0; i <= Lattice.param.Cell_nLoc; i++) {
    if  ((i < isep1) || ((i > isep2) && (i < isdm1)) ||
    ((i > isdm2) && (i < isdac1)) || (i > isdac2)) {
      /* ch normale */
      Lattice.Cell[i].maxampl[X_][0] = -35.e-3;
      Lattice.Cell[i].maxampl[X_][1] =  35.e-3;
      Lattice.Cell[i].maxampl[Y_][1] =  12.5e-3;
    } else if ((i >= isdm1) && (i <= isdm2)) {
      /* SD13 */
      Lattice.Cell[i].maxampl[X_][0] = -21e-3;
      Lattice.Cell[i].maxampl[X_][1] =  21e-3;
      Lattice.Cell[i].maxampl[Y_][1] =   5e-3;
    } else if ((i >= isep1) && (i <= isep2)) {
      /* septum */
      Lattice.Cell[i].maxampl[X_][0] = -25e-3;
      Lattice.Cell[i].maxampl[X_][1] =  25e-3;
      Lattice.Cell[i].maxampl[Y_][1] =  12.5e-3;
    } else if ((i >= isdac1) && (i <= isdac2)) {
      /*  minigap */
      Lattice.Cell[i].maxampl[X_][0] = -35e-3;
      Lattice.Cell[i].maxampl[X_][1] =  35e-3;
      Lattice.Cell[i].maxampl[Y_][1] =  2.5e-3;
    }
    if  (i <= hu600) {
      /* HU640 */
      Lattice.Cell[i].maxampl[Y_][1] =  7.0e-3;
    }
  }

}

/****************************************************************************/
/* void ChamberOn(void)

   Purpose:
     Switch on the vacuum chamber

   Input:
       none

   Output:
       none

   Return:
       none

   Global variables:
       none

   specific functions:
       DefineCh

   Comments:
       none

****************************************************************************/
void ChamberOn(void)
{
  DefineCh();
  status.chambre = true;
}

/****************************************************************************/
/* void Trac_Tab(double x, double px, double y, double py, double dp,
 long nmax, long pos, long *lastn, long *lastpos, FILE *outf1, double Tx[][NTURN])

   Purpose:
       Single particle tracking over NTURN turns
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

   specific functions:
       Cell_Pass

   Comments:
       useful for connection with NAFF

****************************************************************************/
void Trac_Tab(double x, double px, double y, double py, double dp,
	      long nmax, long pos, long &lastn, long &lastpos, FILE *outf1,
	      double Tx[][NTURN])
{
  bool lostF = true; /* Lost particle Flag */
  psVector x1;            /* tracking coordinates */
  long i;
  Vector2  aperture;
  aperture[0] = 1e0; aperture[1] = 1e0;

  x1[0] =  x; x1[1] = px;
  x1[2] =  y; x1[3] = py;
  x1[4] = dp; x1[5] = 0.0;

  lastn = 0;
  (lastpos)=pos;

  Cell_Pass(pos -1, Lattice.param.Cell_nLoc, x1, lastpos);

  if(trace) fprintf(outf1, "\n");

  do {
    (lastn)++;
    if ((lastpos == Lattice.param.Cell_nLoc) &&
        (fabs(x1[0]) < aperture[0]) && (fabs(x1[2]) < aperture[1]))
     /* tracking entre debut anneau et element */
    {
     Cell_Pass(0,Lattice.param.Cell_nLoc, x1, lastpos);
     if(trace) fprintf(outf1, "%6ld %+10.5e %+10.5e %+10.5e %+10.5e"
		       " %+10.5e %+10.5e \n",
		       lastn, x1[0], x1[1], x1[2], x1[3], x1[4], x1[5]);
     i = (lastn)-1;
     Tx[0][i] = x1[0]; Tx[1][i] = x1[1];
     Tx[2][i] = x1[2]; Tx[3][i] = x1[3];
     Tx[4][i] = x1[4]; Tx[5][i] = x1[5];

    }
    else  {
      printf("Trac_Tab: Particle lost \n");
      fprintf(stdout, "%6ld %+10.5g %+10.5g %+10.5g"
	      " %+10.5g %+10.5g %+10.5g \n",
	      lastn, x1[0], x1[1], x1[2], x1[3], x1[4], x1[5]);
      lostF = false;
    }
   }
   while (((lastn) < nmax) && ((lastpos) == Lattice.param.Cell_nLoc) && (lostF == true));


   for (i = 1; i < nmax; i++) {
     fprintf(outf1, "%6ld %+10.5e %+10.5e %+10.5e %+10.5e %+10.5e %+10.5e \n", i,
                     Tx[0][i], Tx[1][i], Tx[2][i], Tx[3][i], Tx[4][i], Tx[5][i]);
   }
}


/****************************************************************************/
/* void NuDx(long Nbx, long Nbz, long Nbtour, double xmax, double zmax,
               double energy)

   Purpose:
       Compute nux, nuz with respect to x : nudx.out   if xmax!=0
                        with respect to z : nudz.out   if zmax!=0
               for an energy offset delta=energy
               over Nbtour turns of the ring
               for x varying within [-xmax, xmax] around the closed orbit
               for z varying within [-zmax, zmax] around the closed orbit

   Input:
       Nbx    horizontal point number
       Nbz    vertical point number
       Nbtour turn number
       xmax   maximum horizontal amplitude
       zmax   maximum vertical amplitude
       energy enrgy offset

   Output:
       none

   Return:
       none

   Global variables:
       none

   Specific functions:
       Trac_Simple, Get_NAFF

   Comments:
       16/01/03 add test for non zero frequency
                add variation around the closed orbit

****************************************************************************/
#define nterm  4
void NuDx(long Nbx, long Nbz, long Nbtour, double xmax, double zmax,
               double energy)
{
  FILE * outf;
  const char ficx[] = "nudx.out";
  const char ficz[] = "nudz.out";
  int i = 0;
  double Tab[6][NTURN], fx[nterm], fz[nterm];
  double x = 0.0 , xp = 0.0 , z = 0.0 , zp = 0.0;
  double x0 = 1e-6, xp0= 0.0 , z0 = 1e-6, zp0 = 0.0;
  double xstep = 0.0, zstep = 0.0;
  double nux = 0.0, nuz = 0.0;
  int nb_freq[2] = {0, 0};
  bool stable = true;
  struct tm *newtime;

  /* Get time and date */
  newtime = GetTime();

    if (trace) printf("Entering NuDx ... results in nudx.out\n\n");

    /////////////
    // H tuneshift
    /////////////

  if (fabs(xmax) > 0.0){
    
    /* Opening file */
    if ((outf = fopen(ficx, "w")) == NULL) {
      fprintf(stdout, "NuDx: error while opening file %s\n", ficx);
      exit_(1);
    }

    fprintf(outf,"# Tracy-2 v. 2.8 -- %s -- %s \n", ficx, asctime2(newtime));
    fprintf(outf,"# nu = f(x) \n");
    fprintf(outf,"#    x[m]          z[m]           fx            fz \n");

    if ((Nbx <= 1) || (Nbz <= 1))
      fprintf(stdout,"NuDx: Error Nbx=%ld Nbz=%ld\n",Nbx,Nbz);

    xstep = xmax/Nbx*2.0;
    x0 = 1e-6 - xmax;
    z0 = 1e-3;

    for (i = 0; i <= Nbx; i++) {
      x  = x0 + i*xstep ;
      xp = xp0 ;
      z  = z0  ;
      zp = zp0 ;

      Trac_Simple(x,xp,z,zp,energy,0.0,Nbtour,Tab,&stable); // tracking around closed orbit
      if (stable) {
        Get_NAFF(nterm, Nbtour, Tab, fx, fz, nb_freq); // gets frequency vectors
        Get_freq(fx,fz,&nux,&nuz);  // gets nux and nuz
      }

      else { // unstable
        nux = 0.0; nuz = 0.0;

      }
      fprintf(outf,"% 10.6e % 10.6e % 10.6e % 10.6e\n",
                    x, z, nux, nuz);
    }
    fclose(outf);
  }

    /////////////
    // V tuneshift
    /////////////

  if (fabs(zmax) > 0.0)
  {
    /* Opening file */
    if ((outf = fopen(ficz, "w")) == NULL) {
      fprintf(stdout, "NuDx: error while opening file %s\n", ficz);
      exit_(1);
    }
    
    fprintf(outf,"# tracy-3.5 -- %s -- %s \n", ficz, asctime2(newtime));
    fprintf(outf,"# nu = f(z) \n");
    fprintf(outf,"#    x[mm]         z[mm]          fx            fz \n");

    zstep = zmax/Nbz*2.0;
    x0 = 1e-3;
    z0 = 1e-6 - zmax;
    for (i = 0; i <= Nbz; i++) {
      x  = x0 ;
      xp = xp0;
      z  = z0 + i*zstep;
      zp = zp0;

      Trac_Simple(x,xp,z,zp,energy,0.0,Nbtour,Tab,&stable);
      if (stable) {
        Get_NAFF(nterm, Nbtour, Tab, fx, fz, nb_freq);
        Get_freq(fx,fz,&nux,&nuz);  // gets nux and nuz
      }
      else {
        nux = 0.0; nuz =0.0;
      }
      fprintf(outf,"% 10.6e % 10.6e % 10.6e % 10.6e\n",
                    x, z, nux, nuz);
    }

    fclose(outf);
  }
}
#undef nterm


double get_D(const double df_x, const double df_y)
{
  double  D;

  const double D_min = -2.0, D_max = -10.0;

  if ((df_x != 0.0) || (df_y != 0.0))
    D = log(sqrt(pow(df_x, 2.0)+pow(df_y, 2.0)))/log(10.0);
  else
    D = D_min;

  return max(D, D_max);
}


/****************************************************************************/
/* void fmap(long Nbx, long Nbz, long Nbtour, double xmax, double zmax,
   double energy, bool diffusion)

   Purpose:
       Compute a frequency map of Nbx x Nbz points
       For each set of initial conditions the particle is tracked over
       Nbtour for an energy offset dp

       The stepsize follows a square root law

       Results in fmap.out

   Input:
       Nbx    horizontal step number
       Nby    vertical step number
       xmax   horizontal maximum amplitude
       zmax   vertical maximum amplitude
       Nbtour number of turn for tracking
       energy particle energy offset

   Output:
       status true if stable
              false otherwise

   Return:
       none

   Global variables:
       none

   Specific functions:
       Trac_Simple, Get_NAFF

   Comments:
       15/10/03 run for the diffusion: nasty patch for retrieving the closed orbit
       16/02/03 patch removed
       
****************************************************************************/
#define NTERM2  10
void fmap(long Nbx, long Nbz, long Nbtour, double xmax, double zmax,
          double energy, bool diffusion, bool matlab)
{
 FILE * outf;
 const char fic[] = "fmap.out";
 long i = 0L, j = 0L;
 double Tab[DIM][NTURN], Tab0[DIM][NTURN];
 double fx[NTERM2], fz[NTERM2], fx2[NTERM2], fz2[NTERM2], dfx, dfz;
 double x = 0.0, xp = 0.0, z = 0.0, zp = 0.0;
 double x0 = 1e-6, xp0 = 0.0, z0 = 1e-6, zp0 = 0.0;
 const double ctau = 0.0;
 double xstep = 0.0, zstep = 0.0;
 double nux1 = 0.0, nuz1 = 0.0, nux2 = 0.0, nuz2 = 0.0;
 int nb_freq[2] = {0, 0};
 long nturn = Nbtour;
 bool status = true;
 struct tm *newtime;

 /* Get time and date */
 time_t aclock;
 time(&aclock);                 /* Get time in seconds */
 newtime = localtime(&aclock);  /* Convert time to struct */

 if (trace) printf("Entering fmap ... results in %s\n\n",fic);

 /* Opening file */
 if ((outf = fopen(fic, "w")) == NULL) {
   fprintf(stdout, "fmap: error while opening file %s\n", fic);
   exit_(1);
 }

 fprintf(outf,"# TRACY II v. 2.6 -- %s -- %s \n", fic, asctime2(newtime));
 fprintf(outf,"# nu = f(x) \n");
 fprintf(outf,"#    x[mm]          z[mm]           fx             fz"
	 "            dfx            dfz      D=log_10(sqrt(df_x^2+df_y^2))\n");

 
 if ((Nbx < 1) || (Nbz < 1))
   fprintf(stdout,"fmap: Error Nbx=%ld Nbz=%ld\n",Nbx,Nbz);
 
 // steps in both planes
 xstep = xmax/sqrt((double)Nbx);
 zstep = zmax/sqrt((double)Nbz);

 // double number of turn if diffusion to compute
 if (diffusion) nturn = 2*Nbtour;

 // px and pz zeroed
 xp = xp0;
 zp = zp0;

// Tracking part + NAFF 
 for (i = -Nbx; i <= Nbx; i++) {
   x  = x0 + sgn(i)*sqrt((double)abs(i))*xstep;
   if (!matlab) fprintf(outf,"\n");
   fprintf(stdout,"\n");
//   for (j = 0; j<= Nbz; j++) {
   for (j = -Nbz; j<= Nbz; j++) {
     z  = z0 + sgn(j)*sqrt((double)abs(j))*zstep;
     // tracking around closed orbit
     Trac_Simple(x,xp,z,zp,energy,ctau,nturn,Tab,&status);
     if (status) { // if trajectory is stable
       // gets frequency vectors
       Get_NAFF(NTERM2, Nbtour, Tab, fx, fz, nb_freq);
       Get_freq(fx,fz,&nux1,&nuz1);  // gets nux and nuz
       if (diffusion) { // diffusion
	 // shift data for second round NAFF
         Get_Tabshift(Tab,Tab0,Nbtour,Nbtour);
	 // gets frequency vectors
         Get_NAFF(NTERM2, Nbtour, Tab0, fx2, fz2, nb_freq);
         Get_freq(fx2,fz2,&nux2,&nuz2); // gets nux and nuz
       }
     } // unstable trajectory
     else { //zeroing output
      nux1 = 0.0; nuz1 = 0.0;
      nux2 = 0.0; nuz2 = 0.0;
     }
     
     // printout value
     if (!diffusion) {
       fprintf(outf,"%14.6e %14.6e %14.6e %14.6e\n",
	       1e3*x, 1e3*z, nux1, nuz1);
       fprintf(stdout,"%14.6e %14.6e %14.6e %14.6e\n",
	       1e3*x, 1e3*z, nux1, nuz1);
     } else {
       dfx = nux1 - nux2; dfz = nuz1 - nuz2;
       fprintf(outf,"%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e\n",
	       1e3*x, 1e3*z, nux1, nuz1, dfx, dfz, get_D(dfx, dfz));
       fprintf(stdout,"%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e\n",
	       1e3*x, 1e3*z, nux1, nuz1, dfx, dfz, get_D(dfx, dfz));
     }

     fflush(outf);
   }
 }

 fclose(outf);
}
#undef NTERM2

/****************************************************************************/
/* void fmapdp(long Nbx, long Nbe, long Nbtour, double xmax, double emax,
   double z, bool *status)

   Purpose:
       Compute a frequency map of Nbx x Nbz points
       For each set of initial conditions the particle is tracked over
       Nbtour for an energy offset dp

       The stepsize follows a square root law

       Results in fmapdp.out

   Input:
       Nbx    horizontal step number
       Nbdp   energy step number
       xmax   horizontal maximum amplitude
       emax   energy
       z      vertical starting amplitude
       Nbtour number of turn for tracking

   Output:
       status true if stable
              false otherwise

   Return:
       none

   Global variables:
       none

   Specific functions:
       Trac_Simple, Get_NAFF

   Comments:
       15/10/03 run for the diffusion: nasty patch for retrieving the closed orbit
       23/10/04 for 6D turn off diffusion automatically and horizontal amplitude
       is negative for negative enrgy offset since this is true for the cod

****************************************************************************/
#define NTERM2  10
void fmapdp(long Nbx, long Nbe, long Nbtour, double xmax, double emax,
              double z, bool diffusion, bool matlab)
{
 FILE * outf;
 const char fic[] = "fmapdp.out";
 long i = 0L, j = 0L;
 double Tab[DIM][NTURN], Tab0[DIM][NTURN];
 double fx[NTERM2], fz[NTERM2], fx2[NTERM2], fz2[NTERM2], dfx, dfz;
 double x = 0.0, xp = 0.0, zp = 0.0, dp = 0.0, ctau = 0.0;
 double x0 = 1e-6, xp0 = 0.0, zp0 = 0.0;
 double xstep = 0.0, estep = 0.0;
 double nux1 = 0.0, nuz1 = 0.0, nux2 = 0.0, nuz2 = 0.0;
 
 int nb_freq[2] = {0, 0};
 long nturn = Nbtour;
 bool status=true;
 struct tm *newtime;

 /* Get time and date */
 time_t aclock;
 time(&aclock);                 /* Get time in seconds */
 newtime = localtime(&aclock);  /* Convert time to struct */

 if (diffusion && Lattice.param.Cavity_on == false) nturn = 2*Nbtour;

 if (trace) printf("Entering fmap ... results in %s\n\n",fic);

 /* Opening file */
 if ((outf = fopen(fic, "w")) == NULL) {
   fprintf(stdout, "fmap: error while opening file %s\n", fic);
   exit_(1);
 }

 fprintf(outf,"# TRACY II v. 2.6 -- %s -- %s \n", fic, asctime2(newtime));
 fprintf(outf,"# nu = f(x) \n");
 fprintf(outf,"#    dp[%%]         x[mm]          fx            fz           dfx           dfz\n");

 if ((Nbx <= 1) || (Nbe <= 1))
   fprintf(stdout,"fmap: Error Nbx=%ld Nbe=%ld\n",Nbx,Nbe);

 xp = xp0;
 zp = zp0;

 xstep = xmax/sqrt((double)Nbx);
 estep = 2.0*emax/Nbe;

 for (i = 0; i <= Nbe; i++) {
   dp  = -emax + i*estep;
   if (!matlab) fprintf(outf,"\n");
   fprintf(stdout,"\n");
//   for (j = 0; j<= Nbx; j++) {
   for (j = -Nbx; j<= Nbx; j++) {

     // IF 6D Tracking diffusion turn off and x negative for dp negative
     if ((Lattice.param.Cavity_on == true) && (dp < 0.0)){
       x  = x0 - sgn(j)*sqrt((double)abs(j))*xstep;
        diffusion = false;
     }   
     else
       x  = x0 + sgn(j)*sqrt((double)abs(j))*xstep;

     Trac_Simple(x,xp,z,zp,dp,ctau,nturn,Tab,&status);
     if (status) {
       Get_NAFF(NTERM2, Nbtour, Tab, fx, fz, nb_freq);
       Get_freq(fx,fz,&nux1,&nuz1);  // gets nux and nuz
       if (diffusion) { // diffusion
         Get_Tabshift(Tab,Tab0,Nbtour,Nbtour); // shift data for second round NAFF
         Get_NAFF(NTERM2, Nbtour, Tab0, fx2, fz2, nb_freq); // gets frequency vectors
         Get_freq(fx2,fz2,&nux2,&nuz2); // gets nux and nuz
       }
     } // unstable trajectory       
     else { //zeroing output
      nux1 = 0.0; nuz1 = 0.0;
      nux2 = 0.0; nuz2 = 0.0;
     }

     // printout value
     if (!diffusion) {
       fprintf(outf,"%14.6e %14.6e %14.6e %14.6e\n",
	       1e2*dp, 1e3*x, nux1, nuz1);
       fprintf(stdout,"%14.6e %14.6e %14.6e %14.6e\n",
	       1e2*dp, 1e3*x, nux1, nuz1);
     } else {
       dfx = nux2 - nux1; dfz = nuz2 - nuz1;
       fprintf(outf,"%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e\n",
	       1e2*dp, 1e3*x, nux1, nuz2, dfx, dfz, get_D(dfx, dfz));
       fprintf(stdout,"%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e\n",
	       1e2*dp, 1e3*x, nux1, nuz2, dfx, dfz, get_D(dfx, dfz));
     }

     fflush(outf);
   }
 }

 fclose(outf);
}
#undef NTERM2

/****************************************************************************/
/* void NuDp(long Nb, long Nbtour, double emax)

   Purpose:
       Computes tunes versus energy offset by tracking
       by linear energy step between -emax and emax

   Input:
       Nb+1   numbers of points
       NbTour number of turns for tracking
       emax   maximum energy

   Output:
       none

   Return:
       none

   Global variables:
       trace

   Specific functions:
       Trac_Simple, Get_NAFF

   Comments:
       none

****************************************************************************/
#define NTERM  4
void NuDp(long Nb, long Nbtour, double emax)
{
  FILE * outf;
  const char fic[] = "nudp.out";
  long i = 0L;
//  long lastpos = 0L;
  double Tab[DIM][NTURN];
  double fx[NTERM], fz[NTERM];
  double x  = 0.0,  xp  = 0.0, z  = 0.0,  zp  = 0.0, ctau  = 0.0, dp  = 0.0;
  double x0 = 1e-6, xp0 = 0.0, z0 = 1e-6, zp0 = 0.0, ctau0 = 0.0, dp0 = 0.0;
  double nux1 = 0.0, nuz1 = 0.0;
  int nb_freq[2] = {0, 0};
  bool status = true;
  struct tm *newtime;

  /* Get time and date */
  newtime = GetTime();

  if (trace) printf("Entering NuDp ...\n\n");

  /* Opening file */
  if ((outf = fopen(fic, "w")) == NULL) {
    fprintf(stdout, "NuDp: error while opening file %s\n", fic);
    exit_(1);
  }

  fprintf(outf,"# TRACY II v. 2.6 -- %s -- %s \n", fic, asctime2(newtime));
  fprintf(outf,"#    dP/P           fx            fz          xcod         pxcod          zcod         pzcod\n");
  fprintf(stdout,"#    dP/P           fx            fz          xcod         pxcod          zcod         pzcod\n");

  if (Nb <= 1L)
    fprintf(stdout,"NuDp: Error Nb=%ld\n",Nb);

  // start loop over energy  
  dp0 = -emax;

  for (i = 0L; i < Nb; i++) {
    dp   = dp0 + i*emax/(Nb-1)*2;
    x    = x0  ;
    xp   = xp0 ;
    z    = z0  ;
    zp   = zp0 ;
    ctau = ctau0;
    
    Trac_Simple(x,xp,z,zp,dp,ctau,Nbtour,Tab,&status); // tracking around closed orbit
    if (status) {
       Get_NAFF(NTERM, Nbtour, Tab, fx, fz, nb_freq); // get frequency vectors
       Get_freq(fx,fz,&nux1,&nuz1);  // gets nux and nuz
    }
    else {
       nux1 = 0.0; nuz1 = 0.0;
       status = true;
    }
    
    //~ getcod(dp, lastpos); // get cod for printout
    findcod(dp);

    fprintf(outf,"%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e\n",
            dp, nux1,nuz1, Lattice.param.CODvect[0], Lattice.param.CODvect[1],
            Lattice.param.CODvect[2], Lattice.param.CODvect[3]);
    fprintf(stdout,"%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e\n",
            dp, nux1,nuz1, Lattice.param.CODvect[0], Lattice.param.CODvect[1],
            Lattice.param.CODvect[2], Lattice.param.CODvect[3]);
  }

  fclose(outf);
}
#undef NTERM




/****************************************************************************/
/* void Phase(double x,double xp,double y, double yp,double energy, double ctau, long Nbtour)

   Purpose:
       Compute 6D phase space
       Results in phase.out

   Input:
       x, xp, y, yp, energy, ctau starting position
       Nbtour turn number

   Output:
       none

   Return:
       none

   Global variables:
       trace

   Specific functions:
       Trac_Simple, Get_NAFF

   Comments:
       none

****************************************************************************/
void Phase(double x,double xp,double y, double yp,double energy, double ctau, long Nbtour)
{
  double Tab[6][NTURN];
  FILE *outf;
  const char fic[] = "phase.out";
  int i;
  bool status;
  struct tm *newtime;

  /* Get time and date */
  newtime = GetTime();

  if (Nbtour > NTURN) {
    fprintf(stdout, "Phase: error Nbtour=%ld > NTURN=%d\n",Nbtour,NTURN);
    exit_(1);
  }

  if ((outf = fopen(fic, "w")) == NULL)  {
    fprintf(stdout, "Phase: error while opening file %s\n", fic);
    exit_(1);
  }

  fprintf(outf,"# TRACY II v. 2.6 -- %s -- %s \n", fic, asctime2(newtime));
  fprintf(outf,"# Phase Space \n");
  fprintf(outf,
  "#    x           xp             z            zp           dp          ctau\n");

  // initialization to zero (case where unstable
  for (i = 0; i < Nbtour; i++) {
    Tab[0][i] = 0.0;
    Tab[1][i] = 0.0;
    Tab[2][i] = 0.0;
    Tab[3][i] = 0.0;
    Tab[4][i] = 0.0;
    Tab[5][i] = 0.0;
  }
  
  Trac_Simple(x,xp,y,yp,energy,ctau,Nbtour,Tab,&status);
  for (i = 0; i < Nbtour; i++) {
    fprintf(outf,"% .5e % .5e % .5e % .5e % .5e % .5e\n",
            Tab[0][i],Tab[1][i],Tab[2][i],Tab[3][i],Tab[4][i],Tab[5][i]);
  }
  fclose(outf);
}

/****************************************************************************/
/* void PhasePoly(long pos, double x0,double px0, double z0, double pz0, double delta0,
               double ctau0, long Nbtour)

   Purpose:
       Compute 6D phase space at position pos (=element number in the lattice )
       for sevelral particle: first aim was for injection study
       Results in phasepoly.out

   Input:
       x, xp, y, yp, energy, ctau starting position
       Nbtour turn number

   Output:
       none

   Return:
       none

   Global variables:
       trace

   Specific functions:
       Trac_Simple, Get_NAFF

   Comments:
       none

****************************************************************************/
void PhasePoly(long pos, double x0,double px0, double z0, double pz0, double delta0,
               double ctau0, long Nbtour)
{
  FILE *outf;
  const char  *fic="phasepoly.out";
  long        lastpos = 0,lastn = 0;
  int         i,j;
  double      x, z, px, pz, delta, ctau;
  double      ex = 1368E-9, el = 1.78E-4;
  double      betax = 9.0, /*betaz = 8.2, */betal = 45.5;
  psVector      xsynch;
  int         nx = 1, ne = 400;
  struct tm   *newtime;

  /* Get time and date */
  newtime = GetTime();

  fprintf(stdout,"Closed orbit:\n");
  fprintf(stdout,"      x            px           z           pz        delta       ctau\n");
  fprintf(stdout,"% 12.8f % 12.8f % 12.8f % 12.8f % 12.8f % 12.8f\n",
          Lattice.param.CODvect[0], Lattice.param.CODvect[1], Lattice.param.CODvect[2],
          Lattice.param.CODvect[3], Lattice.param.CODvect[4], Lattice.param.CODvect[5]);
  lastpos = pos;
  Lattice.param.CODvect = xsynch;
//  xsynch[0] = Lattice.param.CODvect[0];
//  xsynch[1] = Lattice.param.CODvect[1];
//  xsynch[2] = Lattice.param.CODvect[2];
//  xsynch[3] = Lattice.param.CODvect[3];
//  xsynch[4] = Lattice.param.CODvect[4];
//  xsynch[5] = Lattice.param.CODvect[5];
  
  if ((outf = fopen(fic, "w")) == NULL)  {
    fprintf(stdout, "Phase: error while opening file %s\n", fic);
    exit_(1);
  }

  fprintf(outf,"# TRACY II v. 2.6 -- %s -- %s \n", fic, asctime2(newtime));
  fprintf(outf,"# 6D Phase Space \n");
  fprintf(outf,
  "# num         x           xp             z            zp           dp          ctau\n");

  trace = true;
  for (j = 0; j < ne; j++){
    for (i = 0; i < nx; i++){
       x     = x0     + xsynch[0] + sqrt(ex*betax)*cos(2.0*M_PI/nx*i)*0;
       px    = px0    + xsynch[1] + sqrt(ex/betax)*sin(2.0*M_PI/nx*i)*0;
       z     = z0     + xsynch[2];
       pz    = pz0    + xsynch[3];
       delta = delta0 + xsynch[4] + sqrt(el/betal)*sin(2*M_PI/ne*j)*0 ;
       ctau  = ctau0  + xsynch[5] + sqrt(el*betal)*cos(2*M_PI/ne*j)*0 + j*0.002;
       fprintf(outf, "%6ld %+10.5e %+10.5e %+10.5e %+10.5e %+10.5e %+10.5e",
                      0L, x, px, z, pz, delta, ctau);
       Trac(x,px,z,pz,delta,ctau, Nbtour,pos, lastn, lastpos, outf);
       fprintf(outf,"\n");
    }
  }
  fclose(outf);
}

/****************************************************************************/
/* void PhasePortrait(double x0,double px0,double z0, double pz0, double delta0,
                   double end, double Nb, long Nbtour, int num)

   Purpose:
       Compute a phase portrait: Nb orbits
       Results in phaseportrait.out

   Input:
       x0, px0, z0, Pz0, delta0, starting position
       num cooordinate to vary (0 is x and 4 is delta)
       end is the last value for the varying coordinate
       Nb is the number of orbits to draw
       Nbtour turn number

   Output:
       none

   Return:
       none

   Global variables:
       none

   Specific functions:
       Trac_Simple

   Comments:
       Change of tracking routine: do not use a tabular to store data

****************************************************************************/
void PhasePortrait(double x0,double px0,double z0, double pz0, double delta0,
                   double ctau0, double end, long Nb, long Nbtour, int num)
{
  double Tab[6][NTURN];
  FILE *outf;
  const char fic[] = "phaseportrait.out";
  int i = 0, j = 0;
  double start = 0.0, step = 0.0;
  double x = 0.0, px = 0.0, z = 0.0, pz = 0.0, delta = 0.0, ctau = 0.0;
  bool status = true;
  struct tm *newtime;

  /* Get time and date */
  newtime = GetTime();

  if (Nbtour > NTURN) {
    fprintf(stdout, "Phase: error Nbtour=%ld > NTURN=%d\n",Nbtour,NTURN);
    exit_(1);
  }

  if ((outf = fopen(fic, "w")) == NULL) {
    fprintf(stdout, "Phase: error while opening file %s\n", fic);
    exit_(1);
  }

  fprintf(outf,"# TRACY II v. 2.6  -- %s \n", asctime2(newtime));
  fprintf(outf,"#  x           xp            z           zp           dp          ctau\n#\n");
  
  x = x0; px = px0;
  z = z0; pz = pz0;
  delta = delta0; 
  
  switch (num) {
    case 0:
      start = x0; break;
    case 1:
      start = px0; break;
    case 2:
      start = z0; break;
    case 3:
      start = pz0; break;
    case 4:
      start = delta0; break;
    case 5:
      start = ctau0; break;
  }

  /** Step between intila conditions **/
  step = (end - start)/Nb;

  for (j = 0; j <= Nb; j++){
    switch (num){
      case 0:
        x     = start + j*step;  break;
      case 1:
        px    = start + j*step;  break;
      case 2:
        z     = start + j*step;  break;
      case 3:
        pz    = start + j*step;  break;
      case 4:
        delta = start + j*step;  break;
      case 5:
        ctau  = start + j*step;  break;
    }

   fprintf(stdout,"% .5e % .5e % .5e % .5e % .5e % .5e\n",
            x,px,z,pz,delta,ctau);
    Trac_Simple(x,px,z,pz,delta,ctau,Nbtour,Tab,&status);
   for (i = 0; i < Nbtour; i++) {
      fprintf(outf,"% .5e % .5e % .5e % .5e % .5e % .5e\n",
            Tab[0][i],Tab[1][i],Tab[2][i],Tab[3][i],Tab[4][i],Tab[5][i]);
    }
  }
  fclose(outf);
}


/****************************************************************************/
/* void Check_Trac(double x, double px, double y, double py, double dp)

   Purpose:
       Diagnosis for tracking
       Used only for debuging
       Print particle coordinates after each element over 1 single turn

   Input:
       x, px, y, py, dp starting conditions for tracking

   Output:
       none

   Return:
       none

   Global variables:
       trace

   Specific functions:
       Trac_Simple, Get_NAFF

   Comments:
       none

****************************************************************************/
void Check_Trac(double x, double px, double y, double py, double dp)
{
  psVector x1;             /* Tracking coordinates */
  long lastpos = Lattice.param.Cell_nLoc;
  FILE *outf;
  const char fic[] = "check_ampl.out";
  int i;

  if ((outf = fopen(fic, "w")) == NULL)
  {
    fprintf(stdout, "Phase: error while opening file %s\n", fic);
    exit_(1);
  }

  x1[0] =  x; x1[1] = px;
  x1[2] =  y; x1[3] = py;
  x1[4] = dp; x1[5] = 0e0;

  fprintf(outf,"# i    x   xp  z   zp   delta cT \n");

  for (i = 1; i<= Lattice.param.Cell_nLoc; i++)
  {
    Cell_Pass(i,i+1, x1, lastpos);
    fprintf(outf,"%4d % .5e % .5e % .5e % .5e % .5e % .5e\n",
            i, x1[0],x1[1],x1[2],x1[3],x1[4],x1[5]);
  }
}

/****************************************************************************/
/* void Enveloppe(double x, double px, double y, double py, double dp, double nturn)

   Purpose:
       Diagnosis for tracking
       Used only for debuging
       Print particle coordinates after each element over 1 single turn

   Input:
       x, px, y, py, dp starting conditions for tracking

   Output:
       none

   Return:
       none

   Global variables:
       trace

   Specific functions:
       Trac_Simple, Get_NAFF

   Comments:
       none

****************************************************************************/
void Enveloppe(double x, double px, double y, double py, double dp, double nturn)
{
  psVector x1; /* Tracking coordinates */
  long lastpos = Lattice.param.Cell_nLoc;
  FILE *outf;
  const char fic[] = "enveloppe.out";
  int i,j ;
  CellType Cell;

  /* Get cod the delta = energy*/
  //~ getcod(dp, lastpos);
  findcod(dp);

  printf("xcod=%.5e mm zcod=% .5e mm \n", Lattice.param.CODvect[0]*1e3, Lattice.param.CODvect[2]*1e3);

  if ((outf = fopen(fic, "w")) == NULL)
  {
    fprintf(stdout, "Enveloppe: error while opening file %s\n", fic);
    exit_(1);
  }

  x1[0] =  x + Lattice.param.CODvect[0]; x1[1] = px + Lattice.param.CODvect[1];
  x1[2] =  y + Lattice.param.CODvect[2]; x1[3] = py + Lattice.param.CODvect[3];
  x1[4] = dp; x1[5] = 0e0;

  fprintf(outf,"# i    x   xp  z   zp   delta cT \n");

  for (j = 1; j <= nturn; j++)
  {
    for (i = 0; i< Lattice.param.Cell_nLoc; i++)
    {/* loop over full ring */

      getelem(i, &Cell);
      Cell_Pass(i,i+1, x1, lastpos);
      if (lastpos != i+1)
      {
       printf("Unstable motion ...\n"); exit_(1);
      }

      fprintf(outf,"%6.2f % .5e % .5e % .5e % .5e % .5e % .5e\n",
              Cell.S, x1[0],x1[1],x1[2],x1[3],x1[4],x1[5]);
    }
  }
}


/****************************************************************************/
/* void Multipole(void)

   Purpose:
       Set multipole in dipoles, quadrupoles, sextupoles, skew quadrupole,
           horizontal and vertical corrector.

   Input:
       none

   Output:
       none

   Return:
       none

   Global variables:
       trace

   Specific functions:
       getelem, SetKLpar, GetKpar

   Comments:
       Test for short and long quadrupole could be changed using the length
       instead of the name. Maybe more portable, in particular if periodicity
       is broken
       Should be rewritten because list already exists now ..

****************************************************************************/
void Multipole(void)
{
  int i = 0;
  int ndip  = 0,  /* Number of dipoles */
      nquad = 0,  /* Number of quadrupoles */
      nsext = 0,  /* Number of sextupoles  */
      nhcorr= 0,  /* Number of horizontal correctors */
      nvcorr= 0,  /* Number of vertical correctors */
      nqcorr= 0;  /* Number of skew quadrupoles */

  int dlist[500];     /* dipole list */
  int qlist[500];     /* Quadrupole list */
  int slist[500];     /* Sextupole list */
  int hcorrlist[120]; /* horizontal corrector list */
  int vcorrlist[120]; /* vertical corrector list */
  int qcorrlist[120]; /* skew quad list */

  CellType Cell;

  int    mOrder = 0;     /* multipole order */
  double mKL = 0.0 ;     /* multipole integrated strength */
  double corr_strength = 0.0;
  double hcorr[120], vcorr[120], qcorr[120];
  double b2 = 0.0, b3 = 0.0;
  double dBoB2 = 0.0, dBoB3 = 0.0, dBoB4 = 0.0, dBoB5 = 0.0, dBoB6 = 0.0,
         dBoB7 = 0.0, dBoB9 = 0.0, dBoB11 = 0.0, dBoB15 = 0.0, dBoB21 = 0.0,
         dBoB27;
  double dBoB6C = 0.0, dBoB6L = 0.0, dBoB10C = 0.0, dBoB10L = 0.0,
         dBoB14C = 0.0, dBoB14L = 0.0;
  double x0i = 0.0, x02i = 0.0, x03i = 0.0, x04i = 0.0, x05i = 0.0,
         x06i = 0.0,
//         x07i = 0.0,
         x08i = 0.0, x012i = 0.0, x010i = 0.0,
         x018i = 0.0, x024i = 0.0;
  double theta = 0.0, brho = 0.0, dummyf = 0.0 ;
  char *dummy = NULL;
  
  FILE *fi;
  const char fic_hcorr[] = "hcor.dat";
  const char fic_vcorr[] = "vcor.dat";
  const char fic_skew[]  = "vcor.dat";
/*********************************************************/

  printf("Enter multipole ... \n");

/* Make lists of dipoles, quadrupoles and  sextupoles */
  for (i = 0; i <= Lattice.param.Cell_nLoc; i++)
  {
    getelem(i, &Cell); /* get element */

    if (Cell.Elem.Kind == Mpole)
    {
      if (fabs(Cell.Elem.M->irho) > 0.0)
      {
        dlist[ndip] = i;
        ndip++;
        if (trace) printf("%s % f\n",Cell.Elem.Name, Cell.Elem.M->B[0 + HOMmax]);
      }
      else if (fabs(Cell.Elem.M->Bpar[2L + HOMmax]) > 0.0)
      {
        qlist[nquad] = i;
        nquad++;
        if (trace) printf("%s % f\n",Cell.Elem.Name, Cell.Elem.M->Bpar[2L + HOMmax]);
      }
      else if (fabs(Cell.Elem.M->Bpar[3L + HOMmax]) > 0.0)
      {
        slist[nsext] = i;
        nsext++;
        if (trace) printf("%s % f\n",Cell.Elem.Name, Cell.Elem.M->Bpar[3L + HOMmax]);
      }
      else if ( Cell.Elem.Name[0] == 'c' && Cell.Elem.Name[1] == 'h')
      {
        hcorrlist[nhcorr] = i;
        nhcorr++;
        if (trace) printf("%s \n",Cell.Elem.Name);
      }
      else if ( Cell.Elem.Name[0] == 'c' && Cell.Elem.Name[1] == 'v')
      {
        vcorrlist[nvcorr] = i;
        nvcorr++;
        if (trace) printf("%s \n",Cell.Elem.Name);
      }
      else if ( Cell.Elem.Name[0] == 'q' && Cell.Elem.Name[1] == 't')
      {
        qcorrlist[nqcorr] = i;
        nqcorr++;
        if (trace) printf("%s \n",Cell.Elem.Name);
      }
    }
  }

 if (!trace) printf("Elements: ndip=%d nquad=%d  nsext=%d nhcorr=%d nvcorr=%d nqcorr=%d\n",
                     ndip,nquad,nsext,nhcorr,nvcorr,nqcorr);

 /***********************************************************************************/
 /*                                                                                 */
 /***********                Set multipoles for dipole               ****************/
 /*
  *                        x0ni w/ n = p-1 for a 2p-poles
  */
 /***********************************************************************************/
 
  x0i   = 1.0/20e-3;  /* 1/radius */
  x02i  = x0i*x0i;
  x03i  = x02i*x0i;
  x04i  = x02i*x02i;
  x05i  = x04i*x0i;
  x06i  = x03i*x03i;
//  x07i  = x06i*x0i;

  dBoB2 =  1.7e-4*0;  /* gradient */
  dBoB3 = -3.7e-4*0;  /* hexapole */
  dBoB4 = -4.1e-5*0;  /* octupole */
  dBoB5 = -9.6e-5*0;  /* decapole */
  dBoB6 = -5.7e-5*0;  /* 12-poles */
  dBoB7 = -4.3e-5*0;  /* 14-poles */

 for (i = 0; i < ndip; i++)
 {
   getelem(dlist[i], &Cell);
   theta = Cell.Elem.L*Cell.Elem.M->irho;

   /* gradient error */
   mKL = dBoB2*theta*x0i;
   SetKLpar(Cell.Fnum, Cell.Knum, mOrder=2L, mKL);

   if (trace) printf("num= %4d name = %s Fnum = %3d, Knum=%3d theta=% e mKl=% e\n",
		     i, Cell.Elem.Name,Cell.Fnum, Cell.Knum, theta, mKL);

   /* sextupole error */
   mKL = dBoB3*theta*x02i;
   SetKLpar(Cell.Fnum, Cell.Knum, mOrder=3L, mKL);

   /* octupole error */
   mKL = dBoB4*theta*x03i;
   SetKLpar(Cell.Fnum, Cell.Knum, mOrder=4L, mKL);

   /* decapole error */
   mKL = dBoB5*theta*x04i;
   SetKLpar(Cell.Fnum, Cell.Knum, mOrder=5L, mKL);

   /* 12-pole error */
   mKL = dBoB6*theta*x05i;
   SetKLpar(Cell.Fnum, Cell.Knum, mOrder=6L, mKL);

   /* 14-pole error */
   mKL = dBoB7*theta*x06i;
   SetKLpar(Cell.Fnum, Cell.Knum, mOrder=7L, mKL);
 }

 /***********************************************************************************/
 /*                                                                                 */
 /***********                Set multipoles for quadripole           ****************/
 /*
  *                          x0ni w/ n = p-2 for a 2p-poles
  */
 /***********************************************************************************/

 x0i  = 1.0/35e-3;       /* 1/Radius in meters */
 b2   = 0.0;             /* Quadrupole strength */
 x02i = x0i*x0i;
 x04i = x02i*x02i;       /* 12-poles */
 x08i = x04i*x04i;       /* 20-poles */
 x012i= x08i*x04i;       /* 28-poles */

 dBoB6C  =  1.6e-5*0;
 dBoB10C = -1.7e-4*0;
 dBoB6L  =  6.5e-5*0;
 dBoB10L =  1.0e-4*0;
 dBoB14L =  1.0e-4*0;
 dBoB14C =  1.0e-4*0;

 for (i = 0; i < nquad; i++)
 {
   getelem(qlist[i], &Cell);
   b2 = Cell.Elem.L*GetKpar(Cell.Fnum, Cell.Knum, 2L);

   /* 12-pole multipole error */
   if ((strncmp(Cell.Elem.Name,"qp2",3)==0) || (strncmp(Cell.Elem.Name,"qp7",3)==0))
      mKL= b2*dBoB6L*x04i;
   else
      mKL= b2*dBoB6C*x04i;
   SetKLpar(Cell.Fnum, Cell.Knum, mOrder=6L, mKL);
   if (trace) printf("num= %4d name = %s Fnum = %3d, Knum=%3d b2=% e mKl=% e\n",i,
               Cell.Elem.Name,Cell.Fnum, Cell.Knum, b2, mKL);

   /* 20-pole multipole error */
   if ((strncmp(Cell.Elem.Name,"qp2",3)==0) || (strncmp(Cell.Elem.Name,"qp7",3)==0))
     mKL= b2*dBoB10L*x08i;
   else
     mKL= b2*dBoB10C*x08i;
   SetKLpar(Cell.Fnum, Cell.Knum, mOrder=10L, mKL);
   if (trace) printf("num= %4d name = %s Fnum = %3d, Knum=%3d b2=% e mKl=% e\n",i,
               Cell.Elem.Name,Cell.Fnum, Cell.Knum, b2, mKL);

   /* 28-pole multipole error */
   if ((strncmp(Cell.Elem.Name,"qp2",3)==0) || (strncmp(Cell.Elem.Name,"qp7",3)==0))
     mKL= b2*dBoB14L*x012i;
   else
     mKL= b2*dBoB14C*x012i;
   SetKLpar(Cell.Fnum, Cell.Knum, mOrder=14L, mKL);
   if (trace) printf("num= %4d name = %s Fnum = %3d, Knum=%3d b2=% e mKl=% e\n",i,
               Cell.Elem.Name,Cell.Fnum, Cell.Knum, b2, mKL);

 }

 /***********************************************************************************/
 /*                                                                                 */
 /***********              Set multipoles for sextupole              ****************/
 /*
  *                        x0ni w/ n = p-3 for a 2p-poles 
  */
 /***********************************************************************************/
 
  b3    = 0.0;
  x0i   = 1.0/35e-3;
  x02i  = x0i*x0i;
  x04i  = x02i*x02i;
  x06i  = x04i*x02i;   /* 18-poles */
  x012i = x06i*x06i;   /* 30-poles */
  x018i = x012i*x06i;  /* 42-poles */
  x024i = x012i*x012i; /* 54-poles */
  dBoB9  =  8.1e-3*0;
  dBoB15 = -1.1e-3*0;
  dBoB21 = -1.1e-3*1;
  dBoB27 = -1.1e-3*0;

 for (i = 0; i < nsext; i++)
 {
   getelem(slist[i], &Cell);
   b3 = GetKpar(Cell.Fnum, Cell.Knum, 3L);

   /* 18-pole multipole error */
   mKL= b3*dBoB9*x06i;
   SetKLpar(Cell.Fnum, Cell.Knum, mOrder=9L, mKL);
   if (trace) printf("num= %4d name = %s Fnum = %3d, Knum=%3d b3=% e mKl=% e\n",i,
               Cell.Elem.Name,Cell.Fnum, Cell.Knum, b3, mKL);

   /* 30-pole multipole error */
   mKL= b3*dBoB15*x012i;
   SetKLpar(Cell.Fnum, Cell.Knum, mOrder=15L, mKL);
   if (trace) printf("num= %4d name = %s Fnum = %3d, Knum=%3d b3=% e mKl=% e\n",i,
               Cell.Elem.Name,Cell.Fnum, Cell.Knum, b3, mKL);

   /* 42-pole multipole error */
   mKL= b3*dBoB21*x018i;
   SetKLpar(Cell.Fnum, Cell.Knum, mOrder=21L, mKL);
   if (trace) printf("num= %4d name = %s Fnum = %3d, Knum=%3d b3=% e mKl=% e\n",i,
               Cell.Elem.Name,Cell.Fnum, Cell.Knum, b3, mKL);

   /* 54-pole multipole error */
   mKL= b3*dBoB27*x024i;
   SetKLpar(Cell.Fnum, Cell.Knum, mOrder=27L, mKL);
   if (trace) printf("num= %4d name = %s Fnum = %3d, Knum=%3d b3=% e mKl=% e\n",i,
               Cell.Elem.Name,Cell.Fnum, Cell.Knum, b3, mKL);
}

 /***********************************************************************************/
 /*                                                                                 */
 /******            Set multipoles for horizontal correctors         ****************/
 /*
  *                x0ni w/ n = p-1 for a 2p-poles
  */
 /***********************************************************************************/
  x0i   = 1.0/35e-3;  /* 1/radius */
  x02i  = x0i*x0i;
  x03i  = x02i*x0i;
  x04i  = x02i*x02i;
  x05i  = x04i*x0i;
  x06i  = x03i*x03i;
  x010i = x05i*x05i;

  dBoB5  = 0.430*0;  /* decapole */
  dBoB7  = 0.067*0;  /* 14-poles */
  dBoB11 =-0.017*0;  /* 22-poles */

  brho = 2.75/0.299792458; /* magnetic rigidity */

  /* open H corrector file */
  if ((fi = fopen(fic_hcorr,"r")) == NULL)
  {
    fprintf(stdout, "Error while opening file %s \n",fic_hcorr);
    exit_(1);
  }

  for (i = 0; i < nhcorr; i++)
  {
    fscanf(fi,"%le \n", &hcorr[i]);
  }
  fclose(fi); /* close H corrector file */

 for (i = 0; i < nhcorr; i++)
 {
   getelem(hcorrlist[i], &Cell);
   corr_strength = hcorr[i]/brho;

   /* gradient error */
   mKL = dBoB5*corr_strength*x04i;
   SetKLpar(Cell.Fnum, Cell.Knum, mOrder=5L, mKL);

   if (trace) printf("num= %4d name = %s Fnum = %3d, Knum=%3d BL/brho=% e mKl=% e\n",i,
               Cell.Elem.Name,Cell.Fnum, Cell.Knum, corr_strength, mKL);
   /* 14-pole error */
   mKL = dBoB7*corr_strength*x06i;
   SetKLpar(Cell.Fnum, Cell.Knum, mOrder=7L, mKL);

   if (trace) printf("num= %4d name = %s Fnum = %3d, Knum=%3d BL/brho=% e mKl=% e\n",i,
            Cell.Elem.Name,Cell.Fnum, Cell.Knum, corr_strength, mKL);

   /* 22-pole error */
   mKL = dBoB11*corr_strength*x010i;
   SetKLpar(Cell.Fnum, Cell.Knum, mOrder=11, mKL);

   if (trace) printf("num= %4d name = %s Fnum = %3d, Knum=%3d BL/brho=% e mKl=% e\n",i,
               Cell.Elem.Name,Cell.Fnum, Cell.Knum, corr_strength, mKL);
 }

 /***********************************************************************************/
 /*                                                                                 */
 /******            Set multipoles for vertical correctors           ****************/
 /*
  *                    x0ni w/ n = p-1 for a 2p-poles
  */
 /***********************************************************************************/

  x0i   = 1.0/35e-3;  /* 1/radius */
  x02i  = x0i*x0i;
  x03i  = x02i*x0i;
  x04i  = x02i*x02i;
  x05i  = x04i*x0i;
  x06i  = x03i*x03i;
  x010i = x05i*x05i;

  dBoB5  = -0.430*0;  /* decapole */
  dBoB7  =  0.063*0;  /* 14-poles */
  dBoB11 =  0.037*0;  /* 22-poles */

  brho = 2.75/0.299792458; /* magnetic rigidity */

  /* open V corrector file */
  if ((fi = fopen(fic_vcorr,"r")) == NULL)
  {
    fprintf(stdout, "Error while opening file %s \n",fic_vcorr);
    exit_(1);
  }

  for (i = 0; i < nvcorr; i++)
  {
    fscanf(fi,"%s %le %le %le \n", dummy,&dummyf,&dummyf,&vcorr[i]);
  }
  fclose(fi); /* close V corrector file */

 for (i = 0; i < nvcorr; i++)
 {
   getelem(vcorrlist[i], &Cell);
   corr_strength = vcorr[i]/brho;

   /* skew decapole error */
   mKL = dBoB5*corr_strength*x04i;
   SetKLpar(Cell.Fnum, Cell.Knum, mOrder=-5L, mKL);

   if (trace) printf("num= %4d name = %s Fnum = %3d, Knum=%3d BL/brho=% e mKl=% e\n",i,
               Cell.Elem.Name,Cell.Fnum, Cell.Knum, corr_strength, mKL);
   /* skew 14-pole error */
   mKL = dBoB7*corr_strength*x06i;
   SetKLpar(Cell.Fnum, Cell.Knum, mOrder=-7L, mKL);

   if (trace) printf("num= %4d name = %s Fnum = %3d, Knum=%3d BL/brho=% e mKl=% e\n",i,
            Cell.Elem.Name,Cell.Fnum, Cell.Knum, corr_strength, mKL);

   /* skew 22-pole error */
   mKL = dBoB11*corr_strength*x010i;
   SetKLpar(Cell.Fnum, Cell.Knum, mOrder=-11, mKL);

   if (trace) printf("num= %4d name = %s Fnum = %3d, Knum=%3d BL/brho=% e mKl=% e\n",i,
               Cell.Elem.Name,Cell.Fnum, Cell.Knum, corr_strength, mKL);
 }

 /***********************************************************************************/
 /*                                                                                 */
 /******                Set multipoles for skew quadripole           ****************/
 /*
  *                        x0ni w/ n = p-2 for a 2p-poles
  */
 /***********************************************************************************/

 /* Set multipoles for skew quad */
  x0i   = 1.0/35e-3;  /* 1/radius */
  x02i  = x0i*x0i;

  dBoB4  = -0.680*0;  /* Octupole */

  /* open skew quad file */
  if ((fi = fopen(fic_skew,"r")) == NULL)
  {
    fprintf(stdout, "Error while opening file %s \n",fic_skew);
    exit_(1);
  }

  for (i = 0; i < nqcorr; i++)
  {
    fscanf(fi,"%le \n", &qcorr[i]);
  }
  fclose(fi); /* close skew quad file */

 for (i = 0; i < nqcorr; i++)
 {
   getelem(qcorrlist[i], &Cell);

   /* skew octupole */
   mKL = dBoB4*qcorr[i]*x02i;
   SetKLpar(Cell.Fnum, Cell.Knum, mOrder=-4L, mKL);

   if (trace) printf("num= %4d name = %s Fnum = %3d, Knum=%3d BL/brho=% e mKl=% e\n",i,
               Cell.Elem.Name,Cell.Fnum, Cell.Knum, corr_strength, mKL);
 }
}

/****************************************************************************/
/* void SetSkewQuad(void)

   Purpose:
       Set SkewQuad in normal quadrupole
       The name of each quadrupole has to be unique

   Input:
       none

   Output:
       none

   Return:
       none

   Global variables:
       trace

   Specific functions:
       GetElem, SetKLpar, GetKpar

   Comments:
       none

****************************************************************************/
void SetSkewQuad(void)
{
  FILE *fi;
  const char fic_skew[] = "QT-solamor_2_3.dat";
  int i;
  double theta[500]; /* array for skew quad tilt*/
  double b2, mKL;
  CellType Cell;
  long mOrder;

  int nquad = 0;  /* Number of skew quadrupoles */
  int qlist[500];  /* Quadrupole list */

  /* make quadrupole list */
  for (i = 0; i <= Lattice.param.Cell_nLoc; i++)
  {
    getelem(i, &Cell); /* get element */

    if (Cell.Elem.Kind == Mpole)
    {
      if (fabs(Cell.Elem.M->Bpar[2L + HOMmax]) > 0.0)
      {
        qlist[nquad] = i;
        nquad++;
        if (trace) printf("%s % f\n",Cell.Elem.Name,
                           Cell.Elem.M->Bpar[2L + HOMmax]);
      }
    }
  }

  /* open skew quad file */
  if ((fi = fopen(fic_skew,"r")) == NULL)
  {
    fprintf(stdout, "Error while opening file %s \n",fic_skew);
    exit_(1);
  }

  /* read tilt in radians */
  for (i = 0; i < nquad; i++)
  {
    fscanf(fi,"%le \n", &theta[i]);
    theta[i+1] = theta[i];
    i++;
  }
  fclose(fi);


  for (i = 0; i < nquad; i++)
  {
    if (trace) fprintf(stdout,"%le \n", theta[i]);

    getelem(qlist[i], &Cell);

    /* Get KL for a quadrupole */
    b2 = Cell.Elem.L*GetKpar(Cell.Fnum, Cell.Knum, 2L);

    mKL = b2*sin(2*theta[i]);
    SetKLpar(Cell.Fnum, Cell.Knum, mOrder=-2L, mKL);
    mKL = b2*cos(2*theta[i]);
    SetKLpar(Cell.Fnum, Cell.Knum, mOrder=2L, mKL);

    if (trace) printf("num= %4d name = %s Fnum = %3d, Knum=%3d KL=% e, KtiltL=% e\n"
                ,i,
                Cell.Elem.Name,Cell.Fnum, Cell.Knum,
                Cell.Elem.M->Bpar[HOMmax+2],
                Cell.Elem.M->Bpar[HOMmax-2]);
 }
}

/****************************************************************************/
/* void MomentumAcceptance(long deb, long fin,
                           double ep_min, double ep_max, long nstepp,
                           double em_min, double em_max, long nstepm)
   Purpose:
        Compute momemtum acceptance along the ring

   Input:
       deb first element for momentum acceptance
       fin last element for momentum acceptance

       ep_min minimum energy deviation for positive momentum acceptance
       ep_max maximum energy deviation for positive momentum acceptance
       nstepp number of energy steps for positive momentum acceptance

       em_min minimum energy deviation for negative momentum acceptance
       em_max maximum energy deviation for negative momentum acceptance
       nstepm number of energy steps for negative momentum acceptance


       * 1 grande section droite
       * 13 entree premier bend
       * 22 sortie SX4
       * 41 section droite moyenne
       * 173 fin superperiode

   Output:
       output file soleil.out : file of results

   Return:
       none

   Global variables:
       none

   specific functions:
       set_vectorcod

   Comments:
       30/06/03 add fflush(NULL) to force writing at the end to coorect
                unexpected bug: rarely the output file is not finished
       31/07/03 add closed orbit a element: usefull for 6D tracking
                delta_closed_orbite = dp(cavite)/2
       21/10/03 add array for vertical intial conditions using tracking
                removed choice of tracking: now this should be done outside

****************************************************************************/
void MomentumAcceptance(long deb, long fin, double ep_min, double ep_max,
			long nstepp, double em_min, double em_max, long nstepm)
{
  double        dP = 0.0, dp1 = 0.0, dp2 = 0.0;
  long          lastpos = 0L,lastn = 0L;
  long          i = 0L, j = 0L, pos = 0L;
  CellType      Cell, Clost;
  double        x = 0.0, px = 0.0, z = 0.0, pz = 0.0, ctau0 = 0.0, delta = 0.0;
  psVector         x0;
  const long    nturn = 1000L;
  FILE          *outf2, *outf1;
  // Nonzero vertical amplitude
  const double  zmax = 0.3e-3; // 0.3 mm at the ring entrance (element 1)
  double        **tabz0, **tabpz0;
  struct tm     *newtime;  // for time
  psVector        codvector[Cell_nLocMax];
  bool          cavityflag, radiationflag;
  
  x0.zero();

  /* Get time and date */
  newtime = GetTime();
  
  /************************/
  /* Fin des declarations */

  /* File opening for writing */

  outf1 = fopen("phase.out", "w");
  outf2 = fopen("soleil.out", "w");

  fprintf(outf2,"# TRACY II v. 2.6  -- %s \n", asctime2(newtime));
  fprintf(outf2,"#  i        s         dp      s_lost  name_lost \n#\n");

  fprintf(outf1,"# TRACY II v. 2.6  -- %s \n", asctime2(newtime));
  fprintf(outf1,"#  i        x           xp            z           zp           dp          ctau\n#\n");
  

  pos = deb; /* starting position in the ring */

  /***************************************************************/
  fprintf(stdout,"Computing initial conditions ... \n");
  /***************************************************************/

  // cod search has to be done in 4D since in 6D it is zero
  cavityflag = Lattice.param.Cavity_on;
  radiationflag = Lattice.param.radiation;  
  Lattice.param.Cavity_on = false;  /* Cavity on/off */
  Lattice.param.radiation = false;  /* radiation on/off */  

   // Allocation of an array of pointer array
  tabz0  = (double **)malloc((nstepp)*sizeof(double*));
  tabpz0 = (double **)malloc((nstepp)*sizeof(double*));
  if (tabz0 == NULL || tabpz0 == NULL){
    fprintf(stdout,"1 out of memory \n"); return;
  }

  for (i = 1L; i <= nstepp; i++){ // loop over energy
    // Dynamical allocation 0 to nstepp -1
    tabz0[i-1L]  = (double *)malloc((fin+1L)*sizeof(double));
    tabpz0[i-1L] = (double *)malloc((fin+1L)*sizeof(double));
    if (tabz0[i-1L] == NULL || tabpz0[i-1L] == NULL){
      fprintf(stdout,"2 out of memory \n"); return;
    }

    // compute dP
    if (nstepp != 1L) {
      dP = ep_max - (nstepp - i)*(ep_max - ep_min)/(nstepp - 1L);
    }
    else {
      dP = ep_max;
    }

    // find and store closed orbit for dP energy offset
    set_vectorcod(codvector, dP);
       
   // coordinates around closed orbit specially usefull for 6D
    x0[0] = codvector[0][0];
    x0[1] = codvector[0][1];
    x0[2] = codvector[0][2] + zmax;
    x0[3] = codvector[0][3];
    x0[4] = codvector[0][4];
    x0[5] = codvector[0][5];

  if (0) fprintf(stdout,"dP=% e : %e %e %e %e %e %e\n",
          dP,x0[0],x0[1],x0[2],x0[3],x0[4],x0[5]);
    // Store vertical initial conditions
    // case where deb is not element 1
    if (deb > 1L){
       Cell_Pass(1L, deb - 1L, x0, lastpos); // track from 1 to deb-1L element
       j = deb -1L;
       if (lastpos != j){ // look if stable
         tabz0 [i- 1L][j] = 1.0;
         tabpz0[i- 1L][j] = 1.0;
       }
       else{ // stable case
         tabz0 [i - 1L][j] = x0[2] - codvector[deb-1L][2];
         tabpz0[i - 1L][j] = x0[3] - codvector[deb-1L][3];
       }
    }
    else { // case where deb is element 1
      j = deb - 1L;
      tabz0 [i - 1L][j] = x0[2] - codvector[j][2];
      tabpz0[i - 1L][j] = x0[3] - codvector[j][3];
   }

    for (j = deb; j < fin; j++){ // loop over elements
      Cell_Pass(j -1L, j, x0, lastpos);
      if (lastpos != j){ // look if stable
        tabz0 [i - 1L][j] = 1.0;
        tabpz0[i - 1L][j] = 1.0;
      }
      else{ // stable case
        tabz0 [i - 1L][j] = x0[2] - codvector[j][2];
        tabpz0[i - 1L][j] = x0[3] - codvector[j][3];
//        fprintf(stdout,"z0= % e pz0= % e\n", tabz0 [i - 1L][j], tabpz0 [i - 1L][j]);
      }
    }
  }

  Lattice.param.Cavity_on = cavityflag;
  Lattice.param.radiation = radiationflag;

  /***************************************************************/
  fprintf(stdout,"Computing positive momentum acceptance ... \n");
  /***************************************************************/

  do
  {
    //~ getcod(dP=0.0, lastpos);       /* determine closed orbit */
    findcod(dP=0.0);
  getelem(pos,&Cell);
    // coordinates around closed orbit which is non zero for 6D tracking
    x     = Cell.BeamPos[0];
    px    = Cell.BeamPos[1];
    z     = Cell.BeamPos[2];
    pz    = Cell.BeamPos[3];
    delta = Cell.BeamPos[4];
    ctau0 = Cell.BeamPos[5];
    fprintf(stdout,"%3ld %6.4g %6.4g %6.4g %6.4g %6.4g %6.4g\n",
            pos, x, px, z, pz, delta, ctau0);

    dp1 = 0.0;
    dp2 = 0.0;
    i   = 0L;
    do /* Tracking over nturn */
    {
      i++;
      dp1 = dp2;
      if (nstepp != 1L) {
        dp2= ep_max - (nstepp - i)*(ep_max - ep_min)/(nstepp - 1L);
      }
      else {
        dp2 = ep_max;
      }      
      if (trace)  printf("i=%4ld pos=%4ld dp=%6.4g\n",i,pos,dp2);
      if (0) fprintf(stdout,"pos=%4ld z0 =% 10.5f  pz0 =% 10.5f  \n", pos, tabz0[i-1L][pos-1L], tabpz0[i-1L][pos-1L]);
      Trac(x, px, tabz0[i-1L][pos], tabpz0[i-1L][pos-1L], dp2+delta , ctau0, nturn, pos, lastn, lastpos, outf1);
    }
    while (((lastn) == nturn) && (i != nstepp));
    
    if ((lastn) == nturn) dp1 = dp2;

    getelem(lastpos,&Clost);
    getelem(pos,&Cell);
    fprintf(stdout,"pos=%4ld z0 =% 10.5f  pz0 =% 10.5f  \n", pos, tabz0[i-1L][pos-1L], tabpz0[i-1L][pos-1L]);
    fprintf(stdout,"%4ld %10.5f %10.5f %10.5f %*s\n", pos,Cell.S,dp1,Clost.S,5,Clost.Elem.Name);
    fprintf(outf2,"%4ld %10.5f %10.5f %10.5f %*s\n", pos,Cell.S,dp1,Clost.S,5,Clost.Elem.Name);
    pos++;
  }
  while(pos != fin);

  // free memory
  for (i = 1L; i <= nstepp; i++){
    free(tabz0 [i - 1L]);
    free(tabpz0[i - 1L]);
  }
  free(tabz0);
  free(tabpz0);

  /***************************************************************/
  /***************************************************************/
  // NEGATIVE MOMENTUM ACCEPTANCE
  /***************************************************************/
  /***************************************************************/
  
  fprintf(outf2,"\n"); /* A void line */

  pos = deb; /* starting position in the ring */
  
  /***************************************************************/
  fprintf(stdout,"Computing initial conditions ... \n");
  /***************************************************************/

  // cod search has to be done in 4D since in 6D it is zero
  cavityflag        = Lattice.param.Cavity_on;
  radiationflag     = Lattice.param.radiation;
  Lattice.param.Cavity_on = false;  /* Cavity on/off */
  Lattice.param.radiation = false;  /* radiation on/off */  
  
   // Allocation of an array of pointer array
  tabz0  = (double **)malloc((nstepm)*sizeof(double*));
  tabpz0 = (double **)malloc((nstepm)*sizeof(double*));
  if (tabz0 == NULL || tabpz0 == NULL){
    fprintf(stdout,"1 out of memory \n"); return;
  }

  for (i = 1L; i <= nstepm; i++){ // loop over energy
    // Dynamical allocation
    tabz0[i-1L]  = (double *)malloc((fin+1L)*sizeof(double));
    tabpz0[i-1L] = (double *)malloc((fin+1L)*sizeof(double));
    if (tabz0[i-1L] == NULL || tabpz0[i-1L] == NULL){
      fprintf(stdout,"2 out of memory \n"); return;
    }

    // compute dP
    if (nstepm != 1L) {
      dP = em_max - (nstepm - i)*(em_max - em_min)/(nstepm - 1L);
    }
    else {      
      dP = em_max;
    }
    // store closed orbit
    set_vectorcod(codvector, dP);

   // coordinates around closed orbit specially usefull for 6D
    x0[0] = codvector[0][0];
    x0[1] = codvector[0][1];
    x0[2] = codvector[0][2] + zmax;
    x0[3] = codvector[0][3];
    x0[4] = codvector[0][4];
    x0[5] = codvector[0][5];

    // Store vertical initial conditions
    // case where deb is not element 1
    if (deb > 1L){
       Cell_Pass(1L, deb - 1L, x0, lastpos); // track from 1 to deb-1L element
       j = deb -1L;
       if (lastpos != j){ // look if stable
         tabz0 [i- 1L][j] = 1.0;
         tabpz0[i- 1L][j] = 1.0;
       }
       else{ // stable case
         tabz0 [i - 1L][j] = x0[2] - codvector[deb-1L][2];
         tabpz0[i - 1L][j] = x0[3] - codvector[deb-1L][3];
       }
    }
    else { // case where deb is element 1
      j = deb - 1L;
      tabz0 [i - 1L][j] = x0[2] - codvector[j][2];
      tabpz0[i - 1L][j] = x0[3] - codvector[j][3];
//      fprintf(stdout,"z0= % e pz0= % e\n", tabz0 [i - 1L][j], tabpz0 [i - 1L][j]);
   }

    for (j = deb; j < fin; j++){ // loop over elements
      Cell_Pass(j -1L, j, x0, lastpos);
      if (lastpos != j){ // look if stable
        tabz0 [i - 1L][j] = 1.0;
        tabpz0[i - 1L][j] = 1.0;
      }
      else{ // stable case
        tabz0 [i - 1L][j] = x0[2] - codvector[j][2];
        tabpz0[i - 1L][j] = x0[3] - codvector[j][3];
//        fprintf(stdout,"dP= % e pos= %ld z0= % e pz0= % e\n", dP, j, tabz0 [i - 1L][j], tabpz0 [i - 1L][j]);
      }
    }
  }

  Lattice.param.Cavity_on = cavityflag;  
  Lattice.param.radiation = radiationflag;

  /***************************************************************/
  fprintf(stdout,"Computing negative momentum acceptance ... \n");
  /***************************************************************/
    
  do {
    //~ getcod(dP=0.0, lastpos);       /* determine closed orbit */
    findcod(dP=0.0);
  getelem(pos,&Cell);
    // coordinates around closed orbit which is non zero for 6D tracking
    x     = Cell.BeamPos[0];
    px    = Cell.BeamPos[1];
    z     = Cell.BeamPos[2];
    pz    = Cell.BeamPos[3];
    delta = Cell.BeamPos[4];
    ctau0 = Cell.BeamPos[5];
    fprintf(stdout,"%3ld %6.4g %6.4g %6.4g %6.4g %6.4g %6.4g\n",
            pos, x, px, z, pz, delta, ctau0);

    dp1 = 0.0;
    dp2 = 0.0;
    i   = 0L;
    do /* Tracking over nturn */
    {
      i++;
      dp1 = dp2;
      /*
       printf("i= %d, dp1=%g, pos=%d lastn=%d\n", i, dp1,pos, lastn);
      */
      if (nstepm != 1L) {
        dp2= em_max - (nstepm - i)*(em_max - em_min)/(nstepm - 1L);
      }
      else {
        dp2 = em_max;
      }
      if (!trace) printf("i=%4ld pos=%4ld dp=%6.4g\n",i,pos,dp2);
      Trac(x, px, tabz0[i-1L][pos], tabpz0[i-1L][pos-1L], dp2+delta , ctau0, nturn, pos, lastn, lastpos, outf1);
    }
    while (((lastn) == nturn) && (i != nstepm));

    if ((lastn) == nturn) dp1 = dp2;

    getelem(lastpos,&Clost);
    getelem(pos,&Cell);
    if (!trace)  printf("i=%4ld pos=%4ld dp=%6.4g\n",i,pos,dp2);
    fprintf(stdout,"pos=%4ld z0 =% 10.5f  pz0 =% 10.5f  \n", pos, tabz0[i-1L][pos-1L], tabpz0[i-1L][pos-1L]);
    fprintf(stdout,"%4ld %10.5f %10.5f %10.5f %*s\n", pos,Cell.S,dp1,Clost.S, 5, Clost.Elem.Name);
    fprintf(outf2,"%4ld %10.5f %10.5f %10.5f %*s\n", pos,Cell.S,dp1,Clost.S, 5, Clost.Elem.Name);
    pos++;
  }
  while(pos != fin);

  // free memory
  for (i = 1L; i <= nstepp; i++){
    free(tabz0 [i - 1L]);
    free(tabpz0[i - 1L]);
  }
  free(tabz0);
  free(tabpz0);
  
  fflush(NULL); // force writing at the end (BUG??)
  fclose(outf1);
  fclose(outf2);
}
  
/****************************************************************************/
/* set_vectorcod(double codvector[Cell_nLocMax][6], double dP)

   Purpose:
      Store closed orbit computed for a Dp energy offset

   Input:
       dP  offset energy 

   Output:
       codvector : closed orbit all around the ring

   Return:
       none

   Global variables:
       status

   Specific functions:
       getcod

   Comments:
       Does not work for a transfer line

****************************************************************************/
void set_vectorcod(psVector  codvector[], double dP)
{
  long      k = 0;
  CellType  Cell;
  psVector    zerovector;

  zerovector.zero();
  
  //~ getcod(dP, lastpos);  /* determine closed orbit */
  findcod(dP);
  
  if (status.codflag == 1) { /* cod exists */
    for (k = 1L; k <= Lattice.param.Cell_nLoc; k++){
      getelem(k,&Cell);
      codvector[k] = Cell.BeamPos;
    }
    // cod at entrance of the ring is the one at the exit (1-periodicity)
    CopyVec(6L, Cell.BeamPos, codvector[0]);
  }
  else { /* nostable cod */
    for (k = 1L; k <= Lattice.param.Cell_nLoc; k++)
      codvector[k] = zerovector;
  }
}

// LAURENT
/****************************************************************************/
/* void spectrum(long Nbx, long Nbz, long Nbtour, double xmax, double zmax,
   double energy, bool *status)

   Purpose:
       Compute a frequency map of Nbx x Nbz points
       For each set of initial conditions the particle is tracked over
       Nbtour for an energy offset dp

       The stepsize follows a square root law

       Results in fmap.out

   Input:
       Nbx    horizontal step number
       Nby    vertical step number
       xmax   horizontal maximum amplitude
       zmax   vertical maximum amplitude
       Nbtour number of turn for tracking
       energy particle energy offset

   Output:
       status true if stable
              false otherwise

   Return:
       none

   Global variables:
       none

   Specific functions:
       Trac_Simple, Get_NAFF

   Comments:
       15/10/03 run for the diffusion: nasty patch for retrieving the closed orbit

****************************************************************************/
void spectrum(long Nbx, long Nbz, long Nbtour, double xmax, double zmax,
              double energy, bool diffusion)
{
 FILE *xoutf, *zoutf;
 const char xfic[] = "xspectrum.out";
 const char zfic[] = "zspectrum.out";
 long i, j, k;
 #define nterm2  20
 double Tab[6][NTURN], fx[nterm2], fz[nterm2];
 double x = 0.0, xp = 0.0, z = 0.0, zp = 0.0;
 double x0 = 1e-6, xp0 = 0.0, z0 = 1e-6, zp0 = 0.0;
 double xstep = 0.0, zstep = 0.0;
 int nb_freq[2] = {0, 0};
 long nturn = Nbtour;
 bool status=true;
 struct tm *newtime;

 /* Get time and date */
 time_t aclock;
 time(&aclock);                 /* Get time in seconds */
 newtime = localtime(&aclock);  /* Convert time to struct */

 if (diffusion) nturn = 2*Nbtour;

// if (trace) printf("Entering fmap ... results in %s\n\n",fic);

 /* Opening file */
 if ((xoutf = fopen(xfic, "w")) == NULL) {
   fprintf(stdout, "fmap: error while opening file %s\n", xfic);
   exit_(1);
 }

 if ((zoutf = fopen(zfic, "w")) == NULL) {
   fprintf(stdout, "fmap: error while opening file %s\n", zfic);
   exit_(1);
 }

 fprintf(xoutf,"# TRACY II v. 2.6 -- %s -- %s \n", xfic, asctime2(newtime));
 fprintf(zoutf,"# TRACY II v. 2.6 -- %s -- %s \n", zfic, asctime2(newtime));
// fprintf(outf,"# nu = f(x) \n");
// fprintf(outf,"#    x[m]          z[m]           fx            fz           dfx           dfz\n");

 if ((Nbx <= 1) || (Nbz <= 1))
   fprintf(stdout,"fmap: Error Nbx=%ld Nbz=%ld\n",Nbx,Nbz);

 xp = xp0;
 zp = zp0;

 xstep = xmax/sqrt((double)Nbx);
 zstep = zmax/sqrt((double)Nbz);

 for (i = 0; i <= Nbx; i++) {
   x  = x0 + sqrt((double)i)*xstep;
   for (j = 0; j<= Nbz; j++) {
     z  = z0 + sqrt((double)j)*zstep;
     Trac_Simple(x,xp,z,zp,energy,0.0,nturn,Tab,&status);
     if (status) {
      Get_NAFF(nterm2, Nbtour, Tab, fx, fz, nb_freq);
     }
     else {
      fx[0]  = 0.0; fz[0]  = 0.0;
     }

     // printout value
         if (!diffusion){

       fprintf(xoutf,"%14.6e %14.6e", x, z);
       fprintf(zoutf,"%14.6e %14.6e", x, z);
       fprintf(stdout,"%14.6e %14.6e", x, z);

       for (k = 0; k < nb_freq[0]; k++){
         fprintf(xoutf," %14.6e", fx[k]);
         fprintf(stdout," %14.6e", fx[k]);
       }

       for (k = 0; k < nb_freq[1]; k++){
         fprintf(zoutf," %14.6e", fz[k]);
         fprintf(stdout," %14.6e", fz[k]);
       }

       fprintf(stdout,"\n");
       fprintf(xoutf,"\n");
       fprintf(zoutf,"\n");       
     }
//     else {
//       fprintf(outf,"%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e\n",
//        x, z, fx[0], fz[0], fx[0]-fx2[0], fz[0]-fz2[0]);
//       fprintf(stdout,"%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e\n",
//        x, z, fx[0], fz[0], fx[0]-fx2[0], fz[0]-fz2[0]);
//     }
   }
 }

 fclose(xoutf);
 fclose(zoutf);
}

/****************************************************************************/
/* void TracCO(double x, double px, double y, double py, double dp, double ctau,
          long nmax, long pos, long *lastn, long *lastpos, FILE *outf1)

   Purpose:
      Single particle tracking
      Same as Trac but with respect to closed orbit

   Input:
      x, px, y, py 4 transverses coordinates
      dp           energy offset
      nmax         number of turns
      pos          starting position for tracking
      aperture     global physical aperture


   Output:
      lastn       last n (should be nmax if  not lost)
      lastpos     last position in the ring

   Return:
       none

   Global variables:
       Lattice.param

   specific functions:
       Cell_Pass

   Comments:
       BUG: last printout is wrong because not at pos but at the end of the ring
       26/04/03 print output for phase space is for position pos now

****************************************************************************/
void TracCO(double x, double px, double y, double py, double dp, double ctau,
	    long nmax, long pos, long &lastn, long &lastpos, FILE *outf1)
{
  psVector x1;     /* tracking coordinates */
  CellType Cell;

  /* Get closed orbit */
  Lattice.Ring_GetTwiss(true, 0.0);
  Lattice.getcod(dp, lastpos);
  getelem(pos-1,&Cell);

  if (!trace) printf("dp= % .5e %% xcod= % .5e mm zcod= % .5e mm \n",
             dp*1e2, Cell.BeamPos[0]*1e3, Cell.BeamPos[2]*1e3);

  /* Tracking coordinates around the closed orbit */
    x1[0] =  x + Cell.BeamPos[0]; x1[1] = px   + Cell.BeamPos[1];
    x1[2] =  y + Cell.BeamPos[2]; x1[3] = py   + Cell.BeamPos[3];
    x1[4] = dp; x1[5] = ctau; // line true in 4D tracking
//    x1[4] = dp + Cell.BeamPos[4]; x1[5] = ctau + Cell.BeamPos[5];

    lastn = 0;

    (lastpos) = pos;

    if (!trace) fprintf(outf1, "\n");

    do
    {
      (lastn)++;
      if (!trace) { // print initial conditions
        fprintf(outf1, "%6ld %+10.5e %+10.5e %+10.5e %+10.5e"
		" %+10.5e %+10.5e \n",
		lastn, x1[0], x1[1], x1[2], x1[3], x1[4], x1[5]);
      }

      Cell_Pass(pos-1L, Lattice.param.Cell_nLoc, x1, lastpos);
      Cell_Pass(0,pos-1L, x1, lastpos);
    }
    while (((lastn) < nmax) && ((lastpos) == pos-1L));

    if (lastpos != pos-1L)
    {
      printf("TracCO: Particle lost \n");
      fprintf(stdout, "turn=%6ld %+10.5g %+10.5g %+10.5g"
	      " %+10.5g %+10.5g %+10.5g \n",
	      lastn, x1[0], x1[1], x1[2], x1[3], x1[4], x1[5]);
    }
  }


/****************************************************************************/
/*   void getA4antidamping()

   Purpose:

   Input:
       none

   Output:
       none

   Return:
       none

   Global variables:
       none

   specific functions:
       none

   Comments:

****************************************************************************/
void getA4antidamping()
  {
  /* function to get A for anti damping condition */
  /* See publication at ALS for off momentum particle dynamics */

  CellType Cell;
  int qlist[320];
  int nquad=0, i;
  double A = 0.0;

  for (i = 0; i <= Lattice.param.Cell_nLoc; i++)
  {
    getelem(i, &Cell); /* get element */

    if (Cell.Elem.Kind == Mpole)
    {
      if (fabs(Cell.Elem.M->Bpar[2L + HOMmax]) > 0.0)
      {
        qlist[nquad] = i;
        nquad++;
        if (!trace) printf("%s % f\n",
			   Cell.Elem.Name, Cell.Elem.M->Bpar[2L + HOMmax]);
      }
    }
  }
  fprintf(stdout,"Nombre de quadrupoles %d\n", nquad);

  Lattice.Ring_GetTwiss(true, 0.0);
  for (i = 0; i < nquad; i++)
  {
    getelem(qlist[i],&Cell);
    fprintf(stdout,"%d Name = %s L=%g A= %g etax=%g \n",
	    i, Cell.Elem.Name, Cell.Elem.L, A,Cell.Eta[0]);
    A +=
      Cell.Elem.L*2.0*(Cell.Elem.M->Bpar[2L + HOMmax]*Cell.Eta[0])*
      (Cell.Elem.M->Bpar[2L + HOMmax]*Cell.Eta[0]);
    i++;
  }
  fprintf(stdout,"A= %g\n", A*1.706);
  }


/****************************************************************************/
/* void fmapfull(long Nbx, long Nbz, long Nbtour, double xmax, double zmax,
   double energy, bool *status)

   Purpose:
       Compute a frequency map of Nbx x Nbz points
       For each set of initial conditions the particle is tracked over
       Nbtour for an energy offset dp

       The stepsize follows a square root law

       Results in fmap.out

   Input:
       Nbx    horizontal step number
       Nby    vertical step number
       xmax   horizontal maximum amplitude
       zmax   vertical maximum amplitude
       Nbtour number of turn for tracking
       energy particle energy offset

   Output:
       status true if stable
              false otherwise

   Return:
       none

   Global variables:
       none

   Specific functions:
       Trac_Simple, Get_NAFF

   Comments:
       Note enough precision for diffusion

****************************************************************************/
#define NTERM  10
void fmapfull(long Nbx, long Nbz, long Nbtour, double xmax, double zmax,
              double energy, bool diffusion)
{
 FILE * outf;
 const char fic[] = "fmapfull.out";
 int i, j, k;
 double Tab[DIM][NTURN], Tab0[DIM][NTURN];
 double fx[NTERM], fz[NTERM], fx2[NTERM], fz2[NTERM];
 double x  = 0.0, xp = 0.0, z = 0.0, zp = 0.0;
 double x0 = 1e-6, xp0 = 0.0, z0 = 1e-6, zp0 = 0.0;
 double xstep = 0.0, zstep = 0.0;
 int nb_freq[2] = {0, 0};
 double nux1[NTERM], nuz1[NTERM],nux2[NTERM], nuz2[NTERM];
 long nturn = Nbtour;
 bool status=true;
 struct tm *newtime;
 char name[14];

 /* Get time and date */
 time_t aclock;
 time(&aclock);                 /* Get time in seconds */
 newtime = localtime(&aclock);  /* Convert time to struct */

 if (diffusion) nturn = 2*Nbtour;

 if (trace) printf("Entering fmap ... results in %s\n\n",fic);

 /* Opening file */
 if ((outf = fopen(fic, "w")) == NULL) {
   fprintf(stdout, "fmapfull: error while opening file %s\n", fic);
   exit_(1);
 }

 fprintf(outf,"# TRACY II v. 2.6 -- %s -- %s \n", fic, asctime2(newtime));
 fprintf(outf,"# Frequency map freq = f(x,z) \n");
 fprintf(outf,"#    x[m]          z[m]          ");

 for (k = 0; k < NTERM; k++){
   sprintf(name,"f%dx           ",k);
   fprintf(outf,"%s",name);
 }
 for (k = 0; k < NTERM; k++){
   sprintf(name,"f%dz           ",k);
   fprintf(outf,"%s",name);
 }

 if (!diffusion){
   fprintf(outf,"\n");
 }
 else{
   for (k = 0; k < NTERM; k++){
     sprintf(name,"df%dx          ",k);
     fprintf(outf,"%s",name);
   }
   for (k = 0; k < NTERM; k++){
     sprintf(name,"df%dz          ",k);
     fprintf(outf,"%s",name);
   }
   fprintf(outf,"\n");
 }

 if ((Nbx <= 1) || (Nbz <= 1))
   fprintf(stdout,"fmap: Error Nbx=%ld Nbz=%ld\n",Nbx,Nbz);

 xp = xp0;
 zp = zp0;

 xstep = xmax/sqrt((double)Nbx);
 zstep = zmax/sqrt((double)Nbz);

 for (i = 0; i <= Nbx; i++) {
   x  = x0 + sqrt((double)i)*xstep;
   for (j = 0; j<= Nbz; j++) {
     z  = z0 + sqrt((double)j)*zstep;
     Trac_Simple(x,xp,z,zp,energy,0.0,nturn,Tab,&status);

     if (status) {
       Get_NAFF(NTERM, Nbtour, Tab, fx, fz, nb_freq);

       for (k = 0; k < nb_freq[0]; k++){
         nux1[k] = fx[k];
       }
       for (k = 0; k < nb_freq[1]; k++){
         nuz1[k] = fz[k];
       }
       for (k = nb_freq[0]; k < NTERM; k++){
         nux1[k] = 0.0;
       }
       for (k = nb_freq[1]; k < NTERM; k++){
         nuz1[k] = 0.0;
       }          
       if (diffusion){
         Get_Tabshift(Tab,Tab0,Nbtour,Nbtour); // shift data for second round NAFF
         Get_NAFF(NTERM, Nbtour, Tab0, fx2, fz2, nb_freq); // gets frequency vectors

         for (k = 0; k < nb_freq[0]; k++){
           nux2[k] = fx2[k];
         }
         for (k = 0; k < nb_freq[1]; k++){
           nuz2[k] = fz2[k];
         }
         for (k = nb_freq[0]; k < NTERM; k++){
           nux2[k] = 0.0;
         }
         for (k = nb_freq[1]; k < NTERM; k++){
           nuz2[k] = 0.0;
         }
       }
     }
     else {
      for (k = 0; k < NTERM; k++){
        nux1[k] = 0.0;
        nuz1[k] = 0.0;
        nux2[k] = 0.0;
        nuz2[k] = 0.0;
      }
     }
     
     // printout value
     if (!diffusion){
       fprintf(outf,"%14.6e %14.6e ", x, z);
       for (k = 0; k < NTERM; k++){
         fprintf(outf,"%14.6e ", nux1[k]);
       }
       for (k = 0; k < NTERM; k++){
         fprintf(outf,"%14.6e ", nuz1[k]);
       }
       fprintf(outf,"\n");
//       fprintf(stdout,"%14.6e %14.6e %14.6e %14.6e\n", x, z, nux1, nuz1);
     }
     else {
       fprintf(outf,"%14.6e %14.6e ", x, z);
       for (k = 0; k < NTERM; k++){
         fprintf(outf,"%14.6e ", nux1[k]);
       }
       for (k = 0; k < NTERM; k++){
         fprintf(outf,"%14.6e ", nuz1[k]);
       }
       for (k = 0; k < NTERM; k++){
         fprintf(outf,"%14.6e ", nux2[k]);
       }
       for (k = 0; k < NTERM; k++){
         fprintf(outf,"%14.6e ", nuz2[k]);
       }
       fprintf(outf,"\n");
//       fprintf(stdout,"%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e\n",
//        x, z, nux1, nuz1, fx[0]-fx2[0], fz[0]-fz2[0]);
     }
   }
 }

 fclose(outf);
}
#undef NTERM

/****************************************************************************/
/* void Dyna(long Nbx, long Nbz, long Nbtour, double xmax, double zmax,
   double energy, bool *status)

   Purpose:
       Compute a frequency map of Nbx x Nbz points
       For each set of initial conditions the particle is tracked over
       Nbtour for an energy offset dp

       The stepsize follows a square root law

       Results in fmap.out

   Input:
       Nbx    horizontal step number
       Nby    vertical step number
       xmax   horizontal maximum amplitude
       zmax   vertical maximum maplitude
       Nbtour number of turn for tracking
       energy particle energy offset

   Output:
       status true if stable
              false otherwise

   Return:
       none

   Global variables:
       none

   Specific functions:
       Trac_Simple, Get_NAFF

   Comments:
       none

****************************************************************************/
#define NTERM2  2
void Dyna(long Nbx, long Nbz, long Nbtour, double xmax, double zmax,
               double energy, bool diffusion)
{
  FILE * outf;
  const char fic[] = "dyna.out";
  long i, j;
  double Tab[6][NTURN], fx[NTERM2], fz[NTERM2];
  double x = 0.0, xp = 0.0, z = 0.0, zp = 0.0;
  double x0 = 1e-6, xp0 = 0.0, z0 = 1e-6, zp0 = 0.0;
  double xstep = 0.0, zstep = 0.0;
  int nb_freq[2] = {0, 0};
  long nturn = Nbtour;
  bool status=true;
  struct tm *newtime;

  /* Get time and date */
  newtime = GetTime();

  if (diffusion) nturn = 2*Nbtour;

  if (trace) printf("Entering fmap ... results in %s\n\n",fic);

  /* Opening file moustache */
  if ((outf = fopen(fic, "w")) == NULL)
  {
    fprintf(stdout, "fmap: error while opening file %s\n", fic);
    exit_(1);
  }

  fprintf(outf,"# TRACY II v. 2.6 -- %s -- %s \n", fic, asctime2(newtime));
  fprintf(outf,"# nu = f(x) \n");
  fprintf(outf,"#    x[m]          z[m]           fx            fz \n");

  if ((Nbx <= 1) || (Nbz <= 1))
    fprintf(stdout,"fmap: Error Nbx=%ld Nbz=%ld\n",Nbx,Nbz);

  xp = xp0;
  zp = zp0;

  xstep = xmax/sqrt((double)Nbx);
  zstep = zmax/sqrt((double)Nbz);

  for (i = 0; i <= Nbx; i++) {
    x  = x0 + sqrt((double)i)*xstep;
    for (j = 0; j<= Nbz; j++) {
      z  = z0 + sqrt((double)j)*zstep;
      Trac_Simple(x,xp,z,zp,energy,0.0,nturn,Tab,&status);
      if (status) Get_NAFF(NTERM2, Nbtour, Tab, fx, fz, nb_freq);
      else {
       fx[0] = 0.0; fz[0] = 0.0;
      }
      fprintf(outf,"%14.6e %14.6e %14.6e %14.6e %d\n", x, z, fx[0], fz[0], status);
      fprintf(stdout,"%14.6e %14.6e %14.6e %14.6e %d\n", x, z, fx[0], fz[0], status);
      if (diffusion) {
        if (status) Get_NAFF(NTERM2, Nbtour, Tab, fx, fz, nb_freq);
        fprintf(outf,"%14.6e %14.6e %14.6e %14.6e %d\n", x, z, fx[0], fz[0], status);
        fprintf(stdout,"%14.6e %14.6e %14.6e %14.6e %d\n", x, z, fx[0], fz[0], status);
      }
    }
  }

  xp = xp0;
  zp = zp0;

  for (i = 0; i <= Nbx; i++)  {
    x  = x0 - sqrt((double)i)*xstep;
    for (j = 0; j<= Nbz; j++) {
      z  = z0 + sqrt((double)j)*zstep;
      Trac_Simple(x,xp,z,zp,energy,0.0,nturn,Tab,&status);
      if (status) Get_NAFF(NTERM2, Nbtour, Tab, fx, fz, nb_freq);
      else {
       fx[0] = 0.0; fz[0] =0.0;
      }
      fprintf(outf,"%14.6e %14.6e %14.6e %14.6e %d\n", x, z, fx[0], fz[0], status);
      fprintf(stdout,"%14.6e %14.6e %14.6e %14.6e %d\n", x, z, fx[0], fz[0], status);
      if (diffusion) {
        if (status) Get_NAFF(NTERM2, Nbtour, Tab, fx, fz, nb_freq);
        fprintf(outf,"%14.6e %14.6e %14.6e %14.6e\n", x, z, fx[0], fz[0]);
        fprintf(stdout,"%14.6e %14.6e %14.6e %14.6e\n", x, z, fx[0], fz[0]);
      }
    }
  }
  
  fclose(outf);
}

/****************************************************************************/
/* void Phase2(long pos, double x,double xp,double y, double yp,double energy, double ctau,
               long Nbtour)

   Purpose:
       Compute 6D phase space at position pos (=element number in the lattice )
       Results in phase.out

   Input:
       x, xp, y, yp, energy, ctau starting position
       Nbtour turn number

   Output:
       none

   Return:
       none

   Global variables:
       trace

   Specific functions:
       Trac_Simple, Get_NAFF

   Comments:
       none

****************************************************************************/
void Phase2(long pos, double x,double px,double y, double py,double energy,
            double ctau, long Nbtour)
{
  FILE *outf;
  const char fic[] = "phase2.out";
  long lastpos = 0,lastn = 0;
  struct tm *newtime;

  /* Get time and date */
  newtime = GetTime();

  lastpos = pos;

  if ((outf = fopen(fic, "w")) == NULL) {
    fprintf(stdout, "Phase: error while opening file %s\n", fic);
    exit_(1);
  }

  fprintf(outf,"# TRACY II v. 2.6 -- %s -- %s \n", fic, asctime2(newtime));
  fprintf(outf,"# Phase Space \n");
  fprintf(outf,
  "# num         x           xp             z            zp           dp          ctau\n");

  trace = true;
  Trac(x,px,y,py,energy,ctau, Nbtour,pos, lastn, lastpos, outf);
  fclose(outf);
}

void Phase3(long pos, double x,double px,double y, double py,double energy,
            double ctau, long Nbtour)
{
  FILE *outf;
  const char  *fic="phase3.out";
  long        lastpos = 0,lastn = 0;
  struct tm   *newtime;
  psVector      x1;
  
  /* Get time and date */
  newtime = GetTime();

  lastpos = pos;

  if ((outf = fopen(fic, "w")) == NULL) {
    fprintf(stdout, "Phase: error while opening file %s\n", fic);
    exit_(1);
  }

  fprintf(outf,"# TRACY II v. 2.6 -- %s -- %s \n", fic, asctime2(newtime));
  fprintf(outf,"# Phase Space \n");
  fprintf(outf,
  "# num         x           xp             z            zp           dp          ctau\n");

  trace = true;
  x1[0] = x;   x1[1] = px;     x1[2] = y;
  x1[3] = py;  x1[4] = energy; x1[5] = ctau;  
  Cell_Pass(0L, pos-1L, x1, lastpos);

  x  = x1[0];       px= x1[1];   y = x1[2];
  py = x1[3];  energy = x1[4]; ctau =x1[5];
  
  Trac(x,px,y,py,energy,ctau, Nbtour, pos, lastn, lastpos, outf);
  fclose(outf);
}

/****************************************************************************/
/* void PhaseLongitudinalHamiltonien(void)

   Purpose:
       Compute longitudinal phase space from analytical model
                                                         2              3
                                (                   delta          delta  )
      H(phi,delta) =    omegaRF*(dCoC delta + alpha1----- + alpha2*-----  )
                                (                     2              3    )

                       eVRF (                                               )
                     - -----( cos(phi) - cos(phis) + (phi - phis) sin(phis) )
                        ET  (                                               )
                        

      Integration method Ruth integrator H(phi, delta) = A(delta) + B(phi)
      
   Parameters:   
       omegaRF RF frequency/2pi
       eVRF    RF voltage in electron volt
       phis    synchronous phase
       alpha1  first order momentum compaction factor
       alpha2  second order momentum compaction factor
       dCoC    betatron path lengthening

   Input:
       none
               
   Output:
       longitudinale.out

   Return:
       none

   Global variables:
       trace

   Specific functions:
       PassA, PassB, Hsynchrotron

   Comments:
       none

****************************************************************************/
/* SOLEIL value for SOLAMOR2 */
#define alpha1 4.38E-4
#define alpha2 4.49E-3
#define dCoC  0E-6
#define phis  -0.238
#define E 2.75E3
#define eVRF 4
#define T 1.181E-6
#define omegaRF 352.202E6

void PhaseLongitudinalHamiltonien(void)
{
  long i,j;
  const double t = T;        // To get a one turn map
  double phi, delta, H0;
  long imax = 1000L,         // turn number 
       jmax = 25L;          // starting condition number

  /* Constant stepsize for Ruth's and Forest's Integrator */
  /* Laskar's integrator is not a good idea here, since the correction factor is
     not integrable */
  const double D1 = 0.675603595979829E0;
	const double D2 =-0.175603595979829E0;
	const double C2 = 0.135120719195966E1;
	const double C3 =-0.170241438391932E1;
  
  FILE *outf;
  const char fic[] = "longitudinal.out";
  struct tm *newtime;

  /* Get time and date */
  time_t aclock;
  time(&aclock);                 /* Get time in seconds */
  newtime = localtime(&aclock);  /* Convert time to struct */

  if ((outf = fopen(fic, "w")) == NULL)
  {
    fprintf(stdout, "PhaseLongitudinalHamiltonien: error while opening file %s\n", fic);
    exit_(1);
  }
    
  printf("Last stable orbit %f\n", acos(1.0-T*E/eVRF*Hsynchrotron(0.0,-0.098)));  

  fprintf(outf,"# TRACY II v. 2.6  -- %s \n", asctime2(newtime));
  fprintf(outf,"#  i          ctau              dp             DH/H               H \n#\n");

  for (j = 0L; j < jmax; j++)
  {  
    phi = 0.061417777*j; delta = 0.0001;
    H0 = Hsynchrotron(phi,delta);
    fprintf(outf,"%4ld % 16.8f % 16.8f % 16.8e % 16.8f\n",0L,fmod(phi,2.0*M_PI)*0.8512/2.0/M_PI,delta, 0.0, H0);

    for (i = 0L; i < imax; i++){
  // Leap Frog integrator
  //    PassA(&phi, delta, t*0.5);
  //    PassB(phi, &delta, t);
  //    PassA(&phi, delta, t*0.5);
  // 4th order symplectic integrator
      PassA(&phi, delta, t*D1);
      PassB(phi, &delta, t*C2);
      PassA(&phi, delta, t*D2);
      PassB(phi, &delta, t*C3);
      PassA(&phi, delta, t*D2);
      PassB(phi, &delta, t*C2);
      PassA(&phi, delta, t*D1);
      fprintf(outf,"%4ld % 16.8f % 16.8f % 16.8e % 16.8f\n",i,fmod(phi,2.0*M_PI)*0.8512/2.0/M_PI,
              delta,(H0-Hsynchrotron(phi,delta))/H0,Hsynchrotron(phi,delta));
    }
      fprintf(outf,"\n");
  }
  fclose(outf);      
}


/****************************************************************************/
/* void PassA(double *phi, double delta0, double step)

   Purpose:
       Integrate exp(step*liederivativeof(H(delta,phi))
                                                         2              3
                                (                   delta          delta  )
      H(phi,delta) =    omegaRF*(dCoC delta + alpha1----- + alpha2*-----  )
                                (                     2              3    )


   parameters:
       omegaRF RF frequency/2pi
       eVRF    RF voltage in electron volt
       phis    synchronous phase
       alpha1  first order momentum compaction factor
       alpha2  second order momentum compaction factor
       dCoC    betatron path lengthening

   Input:
       phi, delta coordinates
       step stepsize for integration

   Output:
       phi new phase after t=step

   Return:
       none

   Global variables:
       trace

   Specific functions:
       none

   Comments:
       none

****************************************************************************/
void PassA(double *phi, double delta0, double step)
{
  *phi -= omegaRF*2.0*M_PI*(dCoC + alpha1*delta0 + alpha2*delta0*delta0)*step;
}

/****************************************************************************/
/* void PassB(double phi0, double *delta, double step)

   Purpose:
       Integrate exp(step*liederivativeof(H(delta,phi))

                       eVRF (                                               )
      H(phi,delta) = - -----( cos(phi) - cos(phis) + (phi - phis) sin(phis) )
                        ET  (                                               )


   parameters:
       omegaRF RF frequency/2pi
       eVRF    RF voltage in electron volt
       phis    synchronous phase
       alpha1  first order momentum compaction factor
       alpha2  second order momentum compaction factor
       dCoC    betatron path lengthening

   Input:
       phi, delta coordinates
       step stepsize for integration

   Output:
       phi new phase after t=step

   Return:
       none

   Global variables:
       trace

   Specific functions:
       none

   Comments:
       none

****************************************************************************/
void PassB(double phi0, double *delta, double step)
{
  *delta += eVRF/E/T*(sin(phi0) - sin(phis))*step;
}

/****************************************************************************/
/* double Hsynchrotron(double phi, double delta)

   Purpose:
       Compute Hamiltonian
                                                         2              3
                                (                   delta          delta  )
      H(phi,delta) =    omegaRF*(dCoC delta + alpha1----- + alpha2*-----  )
                                (                     2              3    )

                       eVRF (                                               )
                     - -----( cos(phi) - cos(phis) + (phi - phis) sin(phis) )
                        ET  (                                               )


   Input:
       omegaRF RF frequency/2pi
       eVRF    RF voltage in electron volt
       phis    synchronous phase
       alpha1  first order momentum compaction factor
       alpha2  second order momentum compaction factor
       dCoC    betatron path lengthening

   Output:
       none

   Return:
       Hamiltonian computed in phi and delta

   Global variables:
       none

   Specific functions:
       none
       
   Comments:
       none

****************************************************************************/
double Hsynchrotron(double phi, double delta)
{
  double H = 0.0;
  
  H  = omegaRF*2.0*M_PI*(dCoC*delta + alpha1*delta*delta/2.0 + alpha2*delta*delta*delta/3.0);
  H -= eVRF/E/T*(cos(phi) - cos(phis) + (phi-phis)*sin(phis));
  return H;
}


double EnergySmall(double *X, double irho)
{
 double A, B;
 double h = irho;

 A = (1.0+h*X[0])*(X[1]*X[1]+X[3]*X[3])/2.0/(1.0+X[4]);
 B = -h*X[4]*X[0]+h*h*X[0]*X[0]/0.5;
 return (A+B);
}

double EnergyDrift(double *X)
{
 double A;

 A = (X[1]*X[1]+X[3]*X[3])/2.0/(1.0+X[4]);
 return (A);
}

/****************************************************************************/
/* void Enveloppe2(double x, double px, double y, double py, double dp, double nturn)

   Purpose:
       Diagnosis for tracking
       Used only for debuging
       Print particle coordinates after each element over 1 single turn

   Input:
       x, px, y, py, dp starting conditions for tracking

   Output:
       none

   Return:
       none

   Global variables:
       trace

   Specific functions:
       Trac_Simple, Get_NAFF

   Comments:
       none

****************************************************************************/
void Enveloppe2(double x, double px, double y, double py, double dp, double nturn)
{
  psVector x1; /* Tracking coordinates */
  long lastpos = Lattice.param.Cell_nLoc;
  FILE *outf;
  const char fic[] = "enveloppe2.out";
  int i,j ;
  CellType Cell;
  /* Array for Enveloppes */
  double Envxp[Cell_nLocMax], Envxm[Cell_nLocMax];
  double Envzp[Cell_nLocMax], Envzm[Cell_nLocMax];


  /* Get cod the delta = energy*/
  Lattice.getcod(dp, lastpos);
//  /* initialization to chromatic closed orbit */
//  for (i = 0; i<= Lattice.param.Cell_nLoc; i++)
//  {
//   getelem(i, &Cell);
//   Envxm[i] = Cell.BeamPos[0];   Envxp[i] = Cell.BeamPos[0];
//   Envzm[i] = Cell.BeamPos[2];   Envzp[i] = Cell.BeamPos[2];
//  }

  printf("xcod=%.5e mm zcod=% .5e mm \n",
	 Lattice.param.CODvect[0]*1e3, Lattice.param.CODvect[2]*1e3);

  if ((outf = fopen(fic, "w")) == NULL) {
    fprintf(stdout, "Enveloppe: error while opening file %s\n", fic);
    exit_(1);
  }

  x1[0] =  x + Lattice.param.CODvect[0]; x1[1] = px + Lattice.param.CODvect[1];
  x1[2] =  y + Lattice.param.CODvect[2]; x1[3] = py + Lattice.param.CODvect[3];
  x1[4] = dp; x1[5] = 0e0;

  fprintf(outf,"# s       envx(+)       envx(-)       envz(+)       envz(-)"
	  "     delta \n");

  for (i = 0; i < Lattice.param.Cell_nLoc; i++) {
    /* loop over full ring: one turn for intialization */

    getelem(i,&Cell);
    Cell_Pass(i,i+1, x1, lastpos);
    if (lastpos != i+1) {
     printf("Unstable motion ...\n"); exit_(1);
    }

    Envxp[i] = x1[0]; Envxm[i] = x1[0]; Envzp[i] = x1[2]; Envzm[i] = x1[2];
  }

  for (j = 1; j < nturn; j++) {
    /* loop over full ring */
   for (i = 0; i<= Lattice.param.Cell_nLoc; i++) {
 
      getelem(i, &Cell);
      Cell_Pass(i, i+1, x1, lastpos);
      if (lastpos != i+1) {
	printf("Unstable motion ...\n"); exit_(1);
      }
      if (x1[0] >= Envxp[i]) Envxp[i] = x1[0];
      if (x1[0] <= Envxm[i]) Envxm[i] = x1[0];
      if (x1[2] >= Envzp[i]) Envzp[i] = x1[2];
      if (x1[2] <= Envzm[i]) Envzm[i] = x1[2];
      }
  }

  for (i = 0; i <= Lattice.param.Cell_nLoc; i++) {
    getelem(i, &Cell);
    fprintf(outf,"%6.2f % .5e % .5e % .5e % .5e % .5e\n",
            Cell.S, Envxp[i],Envxm[i],Envzp[i],Envzm[i],dp);
  }

}
