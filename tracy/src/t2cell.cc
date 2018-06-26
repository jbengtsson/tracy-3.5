/* Tracy-2

   J. Bengtsson, CBP, LBL      1990 - 1994   Pascal version
                 SLS, PSI      1995 - 1997
   M. Boege      SLS, PSI      1998          C translation
   L. Nadolski   SOLEIL        2002          Link to NAFF, Radia field maps
   J. Bengtsson  NSLS-II, BNL  2004 -        

*/


template<typename T>
inline bool CheckAmpl(const ss_vect<T> &x, const long int loc)
{
  bool not_lost;

  if (Lattice.param.Aperture_on)
    not_lost =
      is_double<T>::cst(x[x_]) > Lattice.Cell[loc]->maxampl[X_][0] &&
      is_double<T>::cst(x[x_]) < Lattice.Cell[loc]->maxampl[X_][1] && 
      fabs(is_double<T>::cst(x[y_])) < Lattice.Cell[loc]->maxampl[Y_][1];
  else
    not_lost =
      is_double<T>::cst(x[x_]) > -max_ampl &&
      is_double<T>::cst(x[x_]) < max_ampl &&
      fabs(is_double<T>::cst(x[y_])) < max_ampl;

  if (!not_lost) {
    if (is_double<T>::cst(x[x_]) < Lattice.Cell[loc]->maxampl[X_][0] ||
        is_double<T>::cst(x[x_]) > Lattice.Cell[loc]->maxampl[X_][1])
      status.lossplane = 1;
    else if (fabs(is_double<T>::cst(x[y_])) > Lattice.Cell[loc]->maxampl[Y_][1])
      status.lossplane = 2;
	    
    if (trace)
      printf("CheckAmpl: Particle lost in plane %d at element:"
	     " %5ld s = %10.5f, x = %12.5e, z= %12.5e\n",
	     status.lossplane, loc, Lattice.Cell[loc]->S,
	     is_double<T>::cst(x[x_]), is_double<T>::cst(x[y_]));
  }

  return not_lost;
}


template<typename T>
void LatticeType::Cell_Pass(const long i0, const long i1, ss_vect<T> &x,
			    long &lastpos)
{
  long int i = 0;

  if (Lattice.param.radiation) Lattice.param.dE = 0e0;

  if (Lattice.param.emittance) {
    I2 = 0e0; I4 = 0e0; I5 = 0e0;

    for (i = 0; i < DOF; i++)
      Lattice.param.D_rad[i] = 0e0;
  }

  if (!CheckAmpl(x, i0))
    lastpos = i0;
  else {
    lastpos = i1;
    for (i = i0; i <= i1; i++) {
      Cell[i]->Propagate(x);
      is_tps<T>::get_ps(x, *Cell[i]);
      if (!CheckAmpl(x, i)) { lastpos = i; break; }
    }
  }
}


void LatticeType::Cell_Pass(const long i0, const long i1, tps &sigma,
			    long &lastpos)
{
  // Note: Sigma_k+1 = M_k Sigma_k M_k^T = (M_k (M_k Sigma_k)^T)^T
  const int n = 9;

  int          i, j, jj[n][nv_tps];
  ss_vect<tps> Id, A;

  const double deps = 1e-20;

  Id.identity();

  map = Id + Lattice.param.CODvect; Lattice.Cell_Pass(0, i0, map, lastpos);

  if (lastpos == i0) {
    map = Id + map.cst(); Lattice.Cell_Pass(i0, i1, map, lastpos);

    if (lastpos == i1) {
      // x_1 = zeta(x_0) => f_1(x) = f_0(zeta^-1(x))

      // deterministic part
      sigma = sigma*Inv(map-map.cst());

      if (Lattice.param.emittance) {
	// stochastic part

	for (i = 0; i < n; i++)
	  for (j = 0; j < nv_tps; j++)
	    jj[i][j] = 0;

	jj[0][x_]  = 2; jj[1][x_]  = 1; jj[1][px_]    = 1; jj[2][px_]    = 2;
	jj[3][y_]  = 2; jj[4][y_]  = 1; jj[4][py_]    = 1; jj[5][py_]    = 2;
	jj[6][ct_] = 2; jj[7][ct_] = 1; jj[7][delta_] = 1; jj[8][delta_] = 2;

	putlinmat(6, Lattice.param.Ascr, A); sigma = sigma*A;

	for (i = 0; i < 3; i++) {
	  if (Lattice.param.eps[i] > deps) {
	    sigma.pook(jj[3*i], sigma[jj[3*i]]-Lattice.param.D_rad[i]/2.0);
	    sigma.pook(jj[3*i+2], sigma[jj[3*i+2]]-Lattice.param.D_rad[i]/2.0);
	  }
	}

	sigma = sigma*Inv(A);
      }
    }
  }
}



bool LatticeType::Cell_getCOD(const long imax, const double eps,
			      const double dP, long &lastpos)
{
  long            j, n, n_iter;
  int             no;
  double          dxabs;
  iVector         jj;
  ss_vect<double> x0, x1, dx;
  ss_vect<tps>    I, dx0, map;

  no = no_tps; danot_(1);
  
  n = (Lattice.param.Cavity_on)? 6 : 4; Lattice.param.dPparticle = dP;

  if (n == 6) {
    // initial guess is zero for 3 D.O.F.
    x0[x_] = 0e0; x0[px_] = 0e0; x0[y_] = 0e0; x0[py_] = 0e0;
    x0[delta_] = 0e0; x0[ct_] = 0e0;
  } else {
    // or eta*dP for 2 1/2 D.O.F.
    // x0[x_] = Lattice.Cell[0]->Eta[X_]*dP;
    // x0[px_] = Lattice.Cell[0]->Etap[X_]*dP;
    // x0[y_] = Lattice.Cell[0]->Eta[Y_]*dP;
    // x0[py_] = Lattice.Cell[0]->Etap[Y_]*dP;
    x0[x_] = 0e0; x0[px_] = 0e0; x0[y_] = 0e0; x0[py_] = 0e0;
    x0[delta_] = dP; x0[ct_] = 0e0;
  }

  for (j = 0; j < ss_dim; j++)
    jj[j] = (j < n)? 1 : 0;

  if (trace) {
    std::cout << std::endl;
    std::cout << "Cell_getCOD:" << std::endl;
    std::cout << std::scientific << std::setprecision(5)
	      << "  0                        x0 ="
	      << std::setw(13) << x0 << std::endl;
  }
  n_iter = 0; I.identity();
  do {
    n_iter++; map.identity(); map += x0;

    Lattice.Cell_Pass(0, Lattice.param.Cell_nLoc, map, lastpos); 

    if (lastpos == Lattice.param.Cell_nLoc) {
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
    Lattice.param.CODvect = x0; getlinmat(6, map, Lattice.param.OneTurnMat);
    Lattice.Cell_Pass(0, Lattice.param.Cell_nLoc, x0, lastpos);
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


bool LatticeType::GetCOD(const long imax, const double eps, const double dP,
			 long &lastpos)
{
  bool cod;

  cod = Cell_getCOD(imax, eps, dP, lastpos);

  return cod;
}


void LatticeType::Cell_Init(void)
{
  bool        Reverse;
  long        i, j;
  int         Fnum, Knum;
  double      Stotal, phi;
  ElemFamType *elemfamp;
  CellType    *cellp;
  MpoleType   *M, *M2;

  // SymbolLength = 15.
  const char first_name[] = "begin          ";
  const int  n_prt        = 10;

  if (debug) printf("\nCell_Init:\n");

  SI_init();

  this->Cell[0] = new MarkerType();
  this->Cell[0]->Elem.Kind = PartsKind(undef);
  this->Cell[0]->Elem.Reverse = false;
  strcpy(this->Cell[0]->Name, first_name);

  for (i = 1; i <= this->param.Elem_nFam; i++) {
    elemfamp = &this->ElemFam[i-1];

    if (elemfamp->CellF->Elem.Kind == PartsKind(Mpole)) {
      M = static_cast<MpoleType*>(elemfamp->CellF);
      memcpy(M->B, M->Bpar, sizeof(mpolArray));
      M->Updateorder();
    }
    
    if (debug)
      printf("  %3ld %1d |%s|",
	     i, elemfamp->CellF->Elem.Kind, elemfamp->CellF->Name);

    for (j = 1; j <= elemfamp->nKid; j++) {
      cellp = this->Cell[elemfamp->KidList[j-1]];
      Fnum = cellp->Fnum; Knum = cellp->Knum; Reverse = cellp->Elem.Reverse;

      if (cellp->Elem.Kind == PartsKind(marker))
	delete cellp;
      else
	printf("\nCell_Init: %d %d is not MarkerType %d\n",
	       cellp->Fnum, cellp->Knum, cellp->Elem.Kind);

      this->Cell[elemfamp->KidList[j-1]] = elemfamp->CellF->clone();
      cellp = this->Cell[elemfamp->KidList[j-1]];
      cellp->Fnum = Fnum; cellp->Knum = Knum; cellp->Elem.Reverse = Reverse;

      if (cellp->Elem.Kind == PartsKind(Mpole)) {
	M2 = static_cast<MpoleType*>(cellp);

	if (Lattice.param.reverse_elem && (cellp->Elem.Reverse == true)) {
	  // Swap entrance and exit angles.
	  printf("Swapping entrance and exit angles for %8s %2ld\n",
		 cellp->Name, i);
	  phi = M->Tx1; M2->Tx1 = M->Tx2; M2->Tx2 = phi; 
	}

	// Set entrance and exit angles.
	cellp->dR[0] = cos(dtor(M2->dRpar));
	cellp->dR[1] = sin(dtor(M2->dRpar));

	// Set displacement to zero.
	cellp->dS[0] = 0e0; cellp->dS[1] = 0e0;

	if (cellp->L != 0e0 || M2->irho != 0e0) {
	  // Thick element or non zero bend radius.
	  M2->thick = thicktype(thick_);
	  // sin(L*irho/2) = sin(theta/2) half the angle.
	  M2->c0 = sin(cellp->L*M2->irho/2e0);
	  // cos roll: sin(theta/2)*cos(dR).
	  M2->c1 = cellp->dR[0]*M2->c0;
	  // sin roll: sin(theta/2)*sin(dR).
	  M2->s1 = cellp->dR[1]*M2->c0;
	} else
	  // Thin lens.
	  M2->thick = thicktype(thin_);
      }

      if (debug) {
	printf(" %4d", elemfamp->KidList[j-1]);
	if (j % n_prt == 0) printf("\n                         ");
      }
    }
    if (debug) printf("\n");
  }

  Stotal = 0e0;
  for (i = 0; i <= this->param.Cell_nLoc; i++) {
    Stotal += this->Cell[i]->L; this->Cell[i]->S = Stotal;
  }
}
