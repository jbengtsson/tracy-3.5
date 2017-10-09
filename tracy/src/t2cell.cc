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
    not_lost =
      is_double<T>::cst(x[x_]) > Lattice.Cell[loc].maxampl[X_][0] &&
      is_double<T>::cst(x[x_]) < Lattice.Cell[loc].maxampl[X_][1] && 
      fabs(is_double<T>::cst(x[y_])) < Lattice.Cell[loc].maxampl[Y_][1];
  else
    not_lost =
      is_double<T>::cst(x[x_]) > -max_ampl &&
      is_double<T>::cst(x[x_]) < max_ampl &&
      fabs(is_double<T>::cst(x[y_])) < max_ampl;

  if (!not_lost) {
    if (is_double<T>::cst(x[x_]) < Lattice.Cell[loc].maxampl[X_][0] ||
        is_double<T>::cst(x[x_]) > Lattice.Cell[loc].maxampl[X_][1])
      status.lossplane = 1;
    else if (fabs(is_double<T>::cst(x[y_])) > Lattice.Cell[loc].maxampl[Y_][1])
      status.lossplane = 2;
	    
    if (trace)
      printf("CheckAmpl: Particle lost in plane %d at element:"
	     " %5ld s = %10.5f, x = %12.5e, z= %12.5e\n",
	     status.lossplane, loc, Lattice.Cell[loc].S,
	     is_double<T>::cst(x[x_]), is_double<T>::cst(x[y_]));
  }

  return not_lost;
}


template<typename T>
void Elem_Pass(const long i, ss_vect<T> &x)
{

  switch (Lattice.Cell[i].Elem.Kind) {
    case drift:
      Drift_Pass(Lattice.Cell[i], x);
      break;
    case Mpole:
      Mpole_Pass(Lattice.Cell[i], x);
      break;
    case Wigl:
      Wiggler_Pass(Lattice.Cell[i], x);
      break;
    case FieldMap:
      FieldMap_Pass(Lattice.Cell[i], x);
      break;
    case Insertion:
      Insertion_Pass(Lattice.Cell[i], x);
      break;
    case Cavity:
      Cav_Pass(Lattice.Cell[i], x);
      break;
    case marker:
      Marker_Pass(Lattice.Cell[i], x);
      break;
    case Spreader:
      break;
    case Recombiner:
      break;
    case Solenoid:
      Solenoid_Pass(Lattice.Cell[i], x);
      break;
    default:
      printf("Elem_Pass ** undefined type\n");
      break;
  }

  is_tps<T>::get_ps(x, Lattice.Cell[i]);
}


template<typename T>
void Cell_Pass(const long i0, const long i1, ss_vect<T> &x, long &lastpos)
{
  long int  i = 0;

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



bool Cell_getCOD(long imax, double eps, double dP, long &lastpos)
{
  long            j, n, n_iter;
  int             no;
  double          dxabs;
  iVector         jj;
  ss_vect<double> x0, x1, dx;
  ss_vect<tps>    I, dx0, map;

  no = no_tps; danot_(1);
  
  n = (globval.Cavity_on)? 6 : 4; globval.dPparticle = dP;

  if (n == 6) {
    // initial guess is zero for 3 D.O.F.
    x0[x_] = 0e0; x0[px_] = 0e0; x0[y_] = 0e0; x0[py_] = 0e0;
    x0[delta_] = 0e0; x0[ct_] = 0e0;
  } else {
    // or eta*dP for 2 1/2 D.O.F.
    x0[x_] = Lattice.Cell[0].Eta[X_]*dP; x0[px_] = Lattice.Cell[0].Etap[X_]*dP;
    x0[y_] = Lattice.Cell[0].Eta[Y_]*dP; x0[py_] = Lattice.Cell[0].Etap[Y_]*dP;
    x0[delta_] = dP; x0[ct_] = 0e0;
  }

  for (j = 0; j < ss_dim; j++)
    jj[j] = (j < n)? 1 : 0;

  if (trace) {
    std::cout << std::endl;
    std::cout << "Cell_getCOD:" << std::endl;
  }
  n_iter = 0; I.identity();
  do {
    n_iter++; map.identity(); map += x0;

    Cell_Pass(0, globval.Cell_nLoc, map, lastpos); 

    if (lastpos == globval.Cell_nLoc) {
      x1 = map.cst(); dx = x0 - x1; dx0 = PInv(map-I-x1, jj)*dx;
      dxabs = xabs(n, dx); x0 += dx0.cst();
    } else {
      dxabs = NAN; break;
    }

    if (trace) {
      std::cout << std::scientific << std::setprecision(1)
	   << std::setw(3) << n_iter
	   << " err = " << std::setw(7) << dxabs << "/" << std::setw(7) << eps
	   << std::setprecision(5)
	   << "  x0 =" << std::setw(13) << x0 << std::endl;
    }
  } while ((dxabs >= eps) && (n_iter <= imax));

  status.codflag = dxabs < eps;

  if (status.codflag) {
    globval.CODvect = x0; getlinmat(6, map, globval.OneTurnMat);
    Cell_Pass(0, globval.Cell_nLoc, x0, lastpos);
  } else {
    std::cout << std::scientific << std::setprecision(5)
	      << "\nCell_getCOD: failed to converge after " << n_iter
	      << " iterations:\n"
	      << "  dP =" << std::setw(12) << dP
	      << ", particle lost at element " << lastpos << "\n"
	      << std::scientific << std::setprecision(5)
	      << "  x0 =" << std::setw(13) << x0 << "\n"
	      << std::scientific << std::setprecision(5)
	      << "   x =" << std::setw(13) << map.cst() << "\n";
  }

  danot_(no);
  
  return status.codflag;
}


bool GetCOD(long imax, double eps, double dP, long &lastpos)
{
  bool  cod;

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

  memcpy(Lattice.Cell[0].Elem.Name, first_name, sizeof(first_name));

  for (i = 1; i <= globval.Elem_nFam; i++) {
    elemfamp  = &Lattice.ElemFam[i-1]; /* Get 1 of all elements stored in
					  ElemFam array */
    elemp = &elemfamp->ElemF; // For switch structure: choice on element type
    if (debug)
      printf("Cell_Init, i:=%3ld: %*s\n", i, SymbolLength, elemp->Name);
    switch (elemp->Kind) {
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
    Stotal += Lattice.Cell[i].Elem.L; Lattice.Cell[i].S = Stotal;
  }
}
