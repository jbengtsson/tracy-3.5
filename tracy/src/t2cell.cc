/* Tracy-2

   J. Bengtsson, CBP, LBL      1990 - 1994   Pascal version
                 SLS, PSI      1995 - 1997
   M. Boege      SLS, PSI      1998          C translation
   L. Nadolski   SOLEIL        2002          Link to NAFF, Radia field maps
   J. Bengtsson  NSLS-II, BNL  2004 -        

*/


long    ntransfmat;
Matrix  transfmat[maxtransfmat];
long    kicks[maxtransfmat][maxkicks];


template<typename T>
inline bool CheckAmpl(const ss_vect<T> &x, const long int loc)
{
  bool  not_lost;

  if (globval.Aperture_on)
    not_lost = is_double<T>::cst(x[x_]) > Cell[loc].maxampl[X_][0] &&
               is_double<T>::cst(x[x_]) < Cell[loc].maxampl[X_][1] && 
               fabs(is_double<T>::cst(x[y_])) < Cell[loc].maxampl[Y_][1];
  else
    not_lost = is_double<T>::cst(x[x_]) > -max_ampl &&
               is_double<T>::cst(x[x_]) < max_ampl &&
               fabs(is_double<T>::cst(x[y_])) < max_ampl;

  if (!not_lost) {
    if (is_double<T>::cst(x[x_]) < Cell[loc].maxampl[X_][0] ||
        is_double<T>::cst(x[x_]) > Cell[loc].maxampl[X_][1])
      status.lossplane = 1;
    else if (fabs(is_double<T>::cst(x[y_])) > Cell[loc].maxampl[Y_][1])
      status.lossplane = 2;
	    
    if (trace)
      printf("CheckAmpl: Particle lost in plane %d at element:"
	     " %5ld s = %10.5f, x = %12.5e, z= %12.5e\n",
	     status.lossplane, loc, Cell[loc].S,
	     is_double<T>::cst(x[x_]), is_double<T>::cst(x[y_]));
  }

  return not_lost;
}


template<typename T>
void Elem_Pass(const long i, ss_vect<T> &x)
{

  switch (Cell[i].Elem.Pkind) {
    case drift:
      Drift_Pass(Cell[i], x);
      break;
    case Mpole:
      Mpole_Pass(Cell[i], x);
      break;
    case Wigl:
      Wiggler_Pass(Cell[i], x);
      break;
    case FieldMap:
      FieldMap_Pass(Cell[i], x);
      break;
    case Insertion:
      Insertion_Pass(Cell[i], x);
      break;
    case Cavity:
      Cav_Pass(Cell[i], x);
      break;
    case marker:
      Marker_Pass(Cell[i], x);
      break;
    case Spreader:
      break;
    case Recombiner:
      break;
    case Solenoid:
      Solenoid_Pass(Cell[i], x);
      break;
    default:
      printf("Elem_Pass ** undefined type\n");
      break;
  }

  is_tps<T>::get_ps(x, Cell[i]);
}


void Elem_Pass_M(const long i, Vector &xref, Matrix &x)
{
  /* Purpose:
       Transport vector xref through matrix x
       xref = Mi(xref)
       x    = M*x   M: transport matrix of element i                         */

  switch (Cell[i].Elem.Pkind) {
    case drift:
      Drift_Pass_M(Cell[i], xref, x);
      break;
    case Mpole:
      Mpole_Pass_M(Cell[i], xref, x);
      break;
    case Wigl:
      Wiggler_Pass_M(Cell[i], xref, x);
      break;
    case Insertion:
      Insertion_Pass_M(Cell[i], xref, x);
      break;
    case Cavity:   /* nothing */
      break;
    case marker:   /* nothing */
      break;
    default:
      fprintf(stdout,"Elem_Pass_M: ** undefined type\n");
      break;
  }

  Cell[i].BeamPos = xref;
}


void Cell_SetdP(const double dP)
{
  int          i, j;
  ElemFamType  *elemfamp;
  elemtype     *elemp;

  globval.dPparticle = dP;
  
  for (i = 1; i <= globval.Elem_nFam; i++) {
    elemfamp  = &ElemFam[i-1]; elemp = &elemfamp->ElemF;
    switch (elemp->Pkind) {
    case drift:
      for (j = 1; j <= elemfamp->nKid; j++)
        Drift_SetMatrix(i, j);
      break;
    case Mpole:
      for (j = 1; j <= elemfamp->nKid; j++)
        Mpole_SetPB(i, j, 2);
      break;
    case Insertion:
      for (j = 1; j <= elemfamp->nKid; j++)
        Insertion_SetMatrix(i, j);
      break;
    case Wigl:
      for (j = 1; j <= elemfamp->nKid; j++)
        Wiggler_SetPB(i, j, 2);
      break;
    case FieldMap:
      break;
    case Cavity:
      break;
    case marker:
      break;
    case Spreader:
      break;
    case Recombiner:
      break;
    case Solenoid:
      break;
    default:
      printf("** Cell_SetdP: undefined type\n");
      break;
    }
  }
}


template<typename T>
void Cell_Pass(const long i0, const long i1, ss_vect<T> &x, long &lastpos)
{
  long int  i = 0;

  if (globval.MatMeth && (x[delta_] != globval.dPparticle))
    Cell_SetdP(is_double<T>::cst(x[delta_]));
    
  if (globval.radiation) globval.dE = 0e0;

  if (globval.emittance) {
    I2 = 0e0; I4 = 0e0; I5 = 0e0;

    for (i = 0; i < DOF; i++)
      globval.D_rad[i] = 0e0;
  }

  if (!CheckAmpl(x, i0))
    lastpos = i0;
  else {
    lastpos = i1;
    for (i = i0; i <= i1; i++) {
      Elem_Pass(i, x);
      if (!CheckAmpl(x, i)) { lastpos = i; break; }
    }
  }
}


void Cell_Pass(const long i0, const long i1, tps &sigma, long &lastpos)
{
  // Note: Sigma_k+1 = M_k Sigma_k M_k^T = (M_k (M_k Sigma_k)^T)^T
  const int  n = 9;

  int           i, j, jj[n][nv_tps];
  ss_vect<tps>  Id, A;

  const double  deps = 1e-20;

  Id.identity();

  map = Id + globval.CODvect; Cell_Pass(0, i0, map, lastpos);

  if (lastpos == i0) {
    map = Id + map.cst(); Cell_Pass(i0, i1, map, lastpos);

    if (lastpos == i1) {
      // x_1 = zeta(x_0) => f_1(x) = f_0(zeta^-1(x))

      // deterministic part
      sigma = sigma*Inv(map-map.cst());

      if (globval.emittance) {
	// stochastic part

	for (i = 0; i < n; i++)
	  for (j = 0; j < nv_tps; j++)
	    jj[i][j] = 0;

	jj[0][x_]  = 2; jj[1][x_]  = 1; jj[1][px_]    = 1; jj[2][px_]    = 2;
	jj[3][y_]  = 2; jj[4][y_]  = 1; jj[4][py_]    = 1; jj[5][py_]    = 2;
	jj[6][ct_] = 2; jj[7][ct_] = 1; jj[7][delta_] = 1; jj[8][delta_] = 2;

	putlinmat(6, globval.Ascr, A); sigma = sigma*A;

	for (i = 0; i < 3; i++) {
	  if (globval.eps[i] > deps) {
	    sigma.pook(jj[3*i], sigma[jj[3*i]]-globval.D_rad[i]/2.0);
	    sigma.pook(jj[3*i+2], sigma[jj[3*i+2]]-globval.D_rad[i]/2.0);
	  }
	}

	sigma = sigma*Inv(A);
      }
    }
  }
}


void Cell_Pass_M(long i0, long i1, Vector &xref, Matrix &mat, long &lastpos)
{
  long i;

  if (xref[4] != globval.dPparticle) Cell_SetdP(xref[4]);

  if (!CheckAmpl(xref, i0))
    lastpos = i0;
  else {
    Cell[i0].BeamPos = xref;
    lastpos = i1;
    for (i = i0+1; i <= i1; i++) {
      Elem_Pass_M(i, xref, mat);
      if (!CheckAmpl(xref, i)) { lastpos = i; break; }
//      CopyVec(6, xref, Cell[i0].BeamPos);
    }
  }
}


#define n 4

bool linearel(long i)
{
  /* Purpose: called by Cell_Concat
       bool function
         true if linear element 
           ie element is:
             straight section
             dipole w/s index (w/o skew quad)
             normal quadrupole (w/o skew quad)
             insertion (focalisation)
             wiggler (focalisation)
             RF cavity
             marker
         false otherwise                                                     */

  bool     status = false;
  CellType *cellp;
  elemtype *elemp;

  cellp = &Cell[i]; elemp = &cellp->Elem;

  switch (elemp->Pkind) {
  case drift: /* straight section */
    status = true;
    break;
  case Mpole:
    if (elemp->M->Pthick == thick && elemp->M->Porder <= Quad &&
	elemp->M->PB[HOMmax-Quad] == 0e0)
      status = true; /* normal quad */
    else
      status = false;
    break;
  case Wigl:
    status = true;
    break;
  case FieldMap:
    status = true;
    break;
  case Insertion:
    status = true;
    break;
  case Cavity:
    status = true;
    break;
  case marker:
    status = true;
    break;
  case Spreader:
    status = false;
    break;
  case Recombiner:
    status = false;
    break;
  case Solenoid:
    status = false;
    break;
  case undef:
    break;
  }
  
  return status;
}

void GtoL_dP(Matrix &mat, Vector2 &dT)
{
  long     k = 0;
  Vector2  dS0;
  Vector   x;

  dS0[0] = 0e0; dS0[1] = 0e0;

  for (k = 0; k < n; k++)
    x[k] = mat[k][n];

  GtoL(x, dS0, dT, 0e0, 0e0, 0e0);

  for (k = 0; k < n; k++)
    mat[k][n] = x[k];
}


static void LtoG_dP(Matrix &mat, Vector2 &dT)
{
  long     k;
  Vector2  dS0;
  Vector   x;

  dS0[0] = 0e0; dS0[1] = 0e0;

  for (k = 0; k < n; k++)
    x[k] = mat[k][n];

  LtoG(x, dS0, dT, 0e0, 0e0, 0e0);

  for (k = 0; k < n; k++)
    mat[k][n] = x[k];
}


void Cell_Concat(double dP)
{
  long      j = 0, i1 = 0;
  double    PB2 = 0e0;
  CellType  *cellp;
  elemtype  *elemp;
  MpoleType *M;

  if (dP != globval.dPparticle) {
    Cell_SetdP(dP); cellconcat = false;
  }
  
  if (cellconcat) return;

  if (trace) printf("concatenating\n");

  cellconcat = true; i1 = 0; ntransfmat = 1;
  UnitMat(n+1, transfmat[ntransfmat-1]);
  kicks[ntransfmat-1][0] = 0;

  do {
    while ((linearel(i1+1) && (i1+1) <= globval.Cell_nLoc)) {
      i1++;
      cellp  = &Cell[i1]; elemp = &cellp->Elem;
      switch (elemp->Pkind) {
        case drift:
          MulLsMat(elemp->D->D55, transfmat[ntransfmat-1]);
          break;

        case Mpole:
          M = elemp->M;
          GtoL_M(transfmat[ntransfmat-1], cellp->dT);
          GtoL_dP(transfmat[ntransfmat-1], cellp->dT);
          GtoL(transfmat[ntransfmat-1][n], cellp->dS, cellp->dT,
          M->Pc0, M->Pc1, M->Ps1);
          if (M->Pthick == thick)
          { /* Drift Kick Drift */
            MulLsMat(M->AU55, transfmat[ntransfmat-1]);
            /* Assuming there is no quadrupole kick */
            PB2 = M->PB[Quad+HOMmax];
            M->PB[Quad+HOMmax] = 0e0;
            thin_kick(M->Porder, M->PB, elemp->PL, 0e0, 0e0,
		      transfmat[ntransfmat-1][n]);
            M->PB[Quad+HOMmax] = PB2;
            MulLsMat(M->AD55, transfmat[ntransfmat-1]);
          } else {
            /* Dipole kick */
            thin_kick(M->Porder, M->PB, 1e0, 0e0, 0e0,
		      transfmat[ntransfmat-1][n]);
          }
          LtoG_M(transfmat[ntransfmat-1], cellp->dT);
          LtoG_dP(transfmat[ntransfmat-1], cellp->dT);
          LtoG(transfmat[ntransfmat-1][n], cellp->dS, cellp->dT,
               M->Pc0, M->Pc1, M->Ps1);
          break;

        case Wigl:
          MulLsMat(elemp->W->W55, transfmat[ntransfmat-1]);
          break;

        case Insertion:
           MulLsMat(elemp->ID->KD55, transfmat[ntransfmat-1]);
          break;

        case Cavity:   /* nothing */
          break;

        case marker:   /* nothing */
          break;

        default:
          printf("**Cell_Concat: undefined type\n");
	  break;
      }
    }
    j = 0;
    while (!linearel(i1+j+1) && (i1+j+1) <= globval.Cell_nLoc) {
      j++;
      if (j >= maxkicks) {
	printf("Cell_Concat maxkicks exceeded: %ld (%d)\n", j, maxkicks);
	exit_(1);
      }
      cellp  = &Cell[i1+j]; elemp = &cellp->Elem;

      if (elemp->Pkind != Mpole) continue;

      M = elemp->M;

      GtoL_M(transfmat[ntransfmat-1], cellp->dT);
      GtoL_dP(transfmat[ntransfmat-1], cellp->dT);
      GtoL(transfmat[ntransfmat-1][n], cellp->dS, cellp->dT, M->Pc0,
           M->Pc1, M->Ps1);

      if (M->Pthick == thick)
        MulLsMat(M->AU55, transfmat[ntransfmat-1]);

      kicks[ntransfmat-1][j-1] = i1 + j; kicks[ntransfmat-1][j] = 0;

      ntransfmat++;
      if (ntransfmat >= maxtransfmat) {
	printf("Cell_Concat maxtransfmat exceeded: %ld (%d)\n",
	       ntransfmat, maxtransfmat);
	exit_(1);
      }
      UnitMat(n+1, transfmat[ntransfmat-1]);
      kicks[ntransfmat-1][0] = 0;

      if (M->Pthick == thick)
        MulLsMat(M->AD55, transfmat[ntransfmat-1]);

      LtoG_M(transfmat[ntransfmat-1], cellp->dT);
      LtoG_dP(transfmat[ntransfmat-1], cellp->dT);
      LtoG(transfmat[ntransfmat-1][n], cellp->dS, cellp->dT, M->Pc0,
	   M->Pc1, M->Ps1);
    }
    i1 += j;
  } while (i1 != globval.Cell_nLoc);
}

#undef n


#define n 4
void Cell_fPass(ss_vect<double> &x, long &lastpos)
{
  /* Purpose:
       Compute the oneturn matrix and propagates xref through it
       Nota: f means full                                                    */

  long       i, j;
  double     PB2;
  CellType   *cellp;
  elemtype   *elemp;
  MpoleType  *M;


  if (!CheckAmpl(x, 1))
    lastpos = 1;
  else { 
    lastpos = globval.Cell_nLoc;
    for (i = 0; i < ntransfmat; i++) {
      LinsTrans(transfmat[i], x);
      j = 0;
      while (kicks[i][j] != 0) {
        j++;
        cellp = &Cell[kicks[i][j-1]]; elemp = &cellp->Elem; M = elemp->M;
        CopyVec(n, x, Cell[kicks[i][j-1]].BeamPos);
        if (M->Pthick == thick) {
          PB2 = M->PB[Quad+HOMmax];
          M->PB[Quad+HOMmax] = 0e0;
          thin_kick(M->Porder, M->PB, elemp->PL, 0e0, 0e0, x);
          M->PB[Quad+HOMmax] = PB2;
        } else
          thin_kick(M->Porder, M->PB, 1e0, 0e0, 0e0, x);

        if (!CheckAmpl(x, kicks[i][j-1])) {
	  lastpos = kicks[i][j-1]; return;
        }
      }
    }
  }
}
#undef n


#define n 4
void Cell_fPass_M(ss_vect<double> &xref, Matrix &mat, long &lastpos)
{
  /* Purpose: called by Cell_MatGetCOD
       Compute the oneturn matrix and propagates xref through it using
       matrix formalism
       Nota: f means full 
       
   Input:
       lastpos last position if unstable
       x  starting vector

   Output:
       mat oneturnmatrix

   Return:
       none

   Global variables:
       transfmat contains transfert matrix for each linear element
       ntransfmat number of transfer matrices
       Cell contains all elements
       
   Specific functions:
        CheckAmpl, MulLsMat, LinsTrans
        thin_kick_M, thin_kick                                               */

  long       i= 0, j = 0;
  double     PB2 = 0e0;
  CellType   *cellp;
  elemtype   *elemp;
  MpoleType  *M;

  if (!CheckAmpl(xref, 1))
    lastpos = 1;
  else {
    lastpos = globval.Cell_nLoc;
    for (i = 0; i < ntransfmat; i++) {
      MulLsMat(transfmat[i], mat); LinsTrans(transfmat[i], xref);
      j = 0;
      while (kicks[i][j] != 0) {
        j++;
        cellp  = &Cell[kicks[i][j-1]]; elemp = &cellp->Elem; M = elemp->M;

        if (M->Pthick == thick) {
          PB2 = M->PB[Quad+HOMmax];
          M->PB[Quad+HOMmax] = 0e0;
          thin_kick_M(M->Porder, M->PB, elemp->PL, 0e0, xref, mat);
          thin_kick(M->Porder, M->PB, elemp->PL, 0e0, 0e0, xref);
          M->PB[Quad+HOMmax] = PB2;
        } else {
          thin_kick_M(M->Porder, M->PB, 1e0, 0e0, xref, mat);
          thin_kick(M->Porder, M->PB, 1e0, 0e0, 0e0, xref);
        }

        if (!CheckAmpl(xref, kicks[i][j-1])) {
	  lastpos = kicks[i][j-1]; return;
        }
      }
    }
  }
}
#undef n


#define n 4
bool Cell_GetCOD_M(long imax, double eps, double dP, long &lastpos)
{
  long    i = 0, j = 0;
  double  dxabs = 0e0;
  Vector  x0, x1, dx;
  Matrix  A;

  if (globval.MatMeth) Cell_Concat(dP);

  CopyVec(n, globval.CODvect, x0);
  x0[delta_] = dP; x0[ct_] = 0e0; i = 0;

  do {
    i++;
    UnitMat(n+2, globval.OneTurnMat); x1 = x0;

    Cell_fPass_M(x1, globval.OneTurnMat, lastpos); /* compute oneturn matrix */
//     Cell_Pass_M(0, globval.Cell_nLoc, x1, globval.OneTurnMat, lastpos);

    if (lastpos == globval.Cell_nLoc) globval.CODvect = x0;

    CopyVec(n, x0, dx); SubVec(n, x1, dx);
    CopyMat(n, globval.OneTurnMat, A);

    /* A = A - Id */
    for (j = 0; j < n; j++)
      A[j][j]--;

    if (InvMat(n, A)) {
      if (lastpos == globval.Cell_nLoc) {
        LinTrans(n, A, dx);
        for (j = 0; j < n; j++)
          x0[j] += dx[j];
      }
    } else
      printf("  *** A is singular\n");

    dxabs = xabs(4, dx);

    if (trace) {
      printf("--- CODLOOP%3ld, Err=% .3E/% .3E\n", i, dxabs, eps);
      printf("% .5E % .5E\n", x0[0], x0[1]);
      printf("% .5E % .5E\n", x0[2], x0[3]);
      printf("% .5E % .5E\n", x0[4], x0[5]);
    }
  } while (i < imax && dxabs > eps && lastpos == globval.Cell_nLoc);

  x0 = globval.CODvect; Cell_Pass(0, globval.Cell_nLoc, x0, lastpos);

  if (dxabs <= eps && lastpos == globval.Cell_nLoc)
    status.codflag = true;
  else {
    printf("Cell_GetCOD_M: GetCOD failed\n");
    status.codflag = false;
  }

  return status.codflag;
}
#undef n


bool Cell_getCOD(long imax, double eps, double dP, long &lastpos)
{
  long             j, n, n_iter;
  int              no;
  double           dxabs;
  iVector          jj;
  ss_vect<double>  x0, x1, dx;
  ss_vect<tps>     I, dx0, map;

  no = no_tps; danot_(1);
  
  n = (globval.Cavity_on)? 6 : 4; globval.dPparticle = dP;

  if (n == 6) {
    // initial guess is zero for 3 D.O.F.
    x0[x_] = 0e0; x0[px_] = 0e0; x0[y_] = 0e0; x0[py_] = 0e0;
    x0[delta_] = 0e0; x0[ct_] = 0e0;
  } else {
    // or eta*dP for 2 1/2 D.O.F.
    x0[x_] = Cell[0].Eta[X_]*dP; x0[px_] = Cell[0].Etap[X_]*dP;
    x0[y_] = Cell[0].Eta[Y_]*dP; x0[py_] = Cell[0].Etap[Y_]*dP;
    x0[delta_] = dP; x0[ct_] = 0e0;
  }

  for (j = 0; j < ss_dim; j++)
    jj[j] = (j < n)? 1 : 0;

  if (trace) {
    cout << endl;
    cout << "Cell_getCOD:" << endl;
  }
  n_iter = 0; I.identity();
  do {
    n_iter++; map.identity(); map += x0;

    Cell_Pass(0, globval.Cell_nLoc, map, lastpos); 

    if (lastpos == globval.Cell_nLoc) {
      x1 = map.cst(); dx = x0 - x1; dx0 = PInv(map-I-x1, jj)*dx;
      dxabs = xabs(n, dx); x0 += dx0.cst();
    } else {
      prt_beampos("beampos.dat");
      prt_cod("cod.out", globval.bpm, true);
      prtmfile("flat_file_dbg.dat");

//      prt_trace();

      dxabs = NAN; break;
    }

    if (trace) {
      cout << scientific << setprecision(1)
	   << setw(3) << n_iter
	   << " err = " << setw(7) << dxabs << "/" << setw(7) << eps
	   << setprecision(5)
	   << "  x0 =" << setw(13) << x0 << endl;
    }
  } while ((dxabs >= eps) && (n_iter <= imax));

  status.codflag = dxabs < eps;

  if (status.codflag) {
    globval.CODvect = x0; getlinmat(6, map, globval.OneTurnMat);
    Cell_Pass(0, globval.Cell_nLoc, x0, lastpos);
  } else {
    cout << "Cell_getCOD: failed to converge after " << n_iter << " iterations"
	  << ", dP=" << setw(13) << dP
	 << ", particle lost at element " << lastpos << endl;
    cout << scientific << setprecision(5)
	 << " x0=" << setw(13) << x0 << endl;
    cout << scientific << setprecision(5)
	 << "  x=" << setw(13) << map.cst() << endl;
  }

  danot_(no);
  
  return status.codflag;
}


bool GetCOD(long imax, double eps, double dP, long &lastpos)
{
  bool  cod;

  if (globval.MatMeth)
    cod = Cell_GetCOD_M(imax, eps, dP, lastpos);
  else
    cod = Cell_getCOD(imax, eps, dP, lastpos);

  return cod;
}


void Cell_Init(void)
{
  long         i;
  double       Stotal;
  ElemFamType  *elemfamp;
  elemtype     *elemp;

  char  first_name[] = "begin          ";

  if (debug)
    printf("**  Cell_Init\n");

  SI_init();  /* Initializes the constants for symplectic integrator */

  memcpy(Cell[0].Elem.PName, first_name, sizeof(first_name));

  for (i = 1; i <= globval.Elem_nFam; i++) {
    elemfamp  = &ElemFam[i-1]; /* Get 1 of all elements stored in ElemFam
				  array */
    elemp = &elemfamp->ElemF; // For switch structure: choice on element type
    if (debug)
      printf("Cell_Init, i:=%3ld: %*s\n", i, SymbolLength, elemp->PName);
    switch (elemp->Pkind) {
    case drift:
      Drift_Init(i);
      break;
      
    case Mpole:
      Mpole_Init(i);
      break;
      
    case Wigl:
      Wiggler_Init(i);
      break;
      
    case FieldMap:
      FieldMap_Init(i);
      break;
      
    case Insertion:
      Insertion_Init(i);
      break;

    case Cavity:
      Cav_Init(i);
      break;

    case marker:
      Marker_Init(i);
      break;

    case Spreader:
      Spreader_Init(i);
      break;

    case Recombiner:
      Recombiner_Init(i);
      break;

    case Solenoid:
      Solenoid_Init(i);
      break;

    default:
      printf("Cell_Init: undefined type\n");
      break;
    }
  }

  /* Computes s-location of each element in the structure */
  Stotal = 0e0;
  for (i = 0; i <= globval.Cell_nLoc; i++) {
    Stotal += Cell[i].Elem.PL; Cell[i].S = Stotal;
  }
}
