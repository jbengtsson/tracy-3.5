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


// template<typename T>
// void LatticeType::Elem_Pass(const long i, ss_vect<T> &x)
// {

//   switch (Lattice.Cell[i]->Kind) {
//     case drift:
//       Lattice.Cell[i]->Pass(x);
//       break;
//     case Mpole:
//       Lattice.Cell[i]->Pass(x);
//       break;
//     case Wigl:
//       Lattice.Cell[i]->Pass(x);
//       break;
//     case FieldMap:
//       Lattice.Cell[i]->Pass(x);
//       break;
//     case Insertion:
//       Lattice.Cell[i]->Pass(x);
//       break;
//     case Cavity:
//       Lattice.Cell[i]->Pass(x);
//       break;
//     case marker:
//       Lattice.Cell[i]->Pass(x);
//       break;
//     case Spreader:
//       break;
//     case Recombiner:
//       break;
//     case Solenoid:
//       Lattice.Cell[i]->Pass(x);
//       break;
//     default:
//       printf("Elem_Pass ** undefined type\n");
//       break;
//   }

//   is_tps<T>::get_ps(x, Lattice.Cell[i]);
// }


template<typename T>
void LatticeType::Elem_Pass(const long i, ss_vect<T> &x)
{
  Cell[i]->Pass(x);
}


template<typename T>
void LatticeType::Cell_Pass(const long i0, const long i1, ss_vect<T> &x,
			    long &lastpos)
{
  long int  i = 0;

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
      Elem_Pass(i, x);
      if (!CheckAmpl(x, i)) { lastpos = i; break; }
    }
  }
}


void LatticeType::Cell_Pass(const long i0, const long i1, tps &sigma,
			    long &lastpos)
{
  // Note: Sigma_k+1 = M_k Sigma_k M_k^T = (M_k (M_k Sigma_k)^T)^T
  const int  n = 9;

  int           i, j, jj[n][nv_tps];
  ss_vect<tps>  Id, A;

  const double  deps = 1e-20;

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
    x0[x_] = Lattice.Cell[0]->Eta[X_]*dP; x0[px_] = Lattice.Cell[0]->Etap[X_]*dP;
    x0[y_] = Lattice.Cell[0]->Eta[Y_]*dP; x0[py_] = Lattice.Cell[0]->Etap[Y_]*dP;
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
  bool  cod;

  cod = Cell_getCOD(imax, eps, dP, lastpos);

  return cod;
}


void LatticeType::Cell_Init(void)
{
  long        i, j;
  double      Stotal;
  ElemFamType *elemfamp;

  // SymbolLength = 15.
  const char first_name[] = "begin          ";
  const int  n_prt        = 10;

  if (debug) printf("\nCell_Init:\n");

  SI_init();

  strcpy(this->Cell[0]->Name, first_name);

  for (i = 1; i <= this->param.Elem_nFam; i++) {
    elemfamp = &this->ElemFam[i-1];
    if (debug)
      printf("  %2ld |%*s|", i, SymbolLength, elemfamp->Name);
    // this->ElemFam[i-1].CellF.Init(i);
    for (j = 1; j <= elemfamp->nKid; j++) {
      if (Cell[elemfamp->KidList[j-1]]->Kind == PartsKind(marker))
	this->Cell[elemfamp->KidList[j-1]] = new MarkerType();
      else if (Cell[elemfamp->KidList[j-1]]->Kind == PartsKind(drift))
	this->Cell[elemfamp->KidList[j-1]] = new DriftType();
      else if (Cell[elemfamp->KidList[j-1]]->Kind == PartsKind(Mpole))
	this->Cell[elemfamp->KidList[j-1]] = new MpoleType();
      else if (Cell[elemfamp->KidList[j-1]]->Kind == PartsKind(Wigl))
	this->Cell[elemfamp->KidList[j-1]] = new WigglerType();
      else if (Cell[elemfamp->KidList[j-1]]->Kind == PartsKind(Insertion))
	this->Cell[elemfamp->KidList[j-1]] = new InsertionType();
      else if (Cell[elemfamp->KidList[j-1]]->Kind == PartsKind(FieldMap))
	this->Cell[elemfamp->KidList[j-1]] = new FieldMapType();
      else if (Cell[elemfamp->KidList[j-1]]->Kind == PartsKind(Cavity))
	this->Cell[elemfamp->KidList[j-1]] = new CavityType();
      else if (Cell[elemfamp->KidList[j-1]]->Kind == PartsKind(Spreader))
	this->Cell[elemfamp->KidList[j-1]] = new SpreaderType();
      else if (Cell[elemfamp->KidList[j-1]]->Kind == PartsKind(Recombiner))
	this->Cell[elemfamp->KidList[j-1]] = new RecombinerType();
      else if (Cell[elemfamp->KidList[j-1]]->Kind == PartsKind(Solenoid))
	this->Cell[elemfamp->KidList[j-1]] = new SolenoidType();
      else {
	printf("\nCell_Init: unknown Kind %d\n",
	       Cell[elemfamp->KidList[j-1]]->Kind);
	exit(1);
      }
      *this->Cell[elemfamp->KidList[j-1]] = *elemfamp->CellF;
      if (debug) {
	printf(" %3d", elemfamp->KidList[j-1]);
	if (j % n_prt == 0) printf("\n                      ");
      }
    }
    if (debug) printf("\n");
  }

  Stotal = 0e0;
  for (i = 0; i <= this->param.Cell_nLoc; i++) {
    Stotal += this->Cell[i]->L; this->Cell[i]->S = Stotal;
  }
}


// Instantiate
template void LatticeType::Cell_Pass(const long i0, const long i1,
				     ss_vect<double> &x, long &lastpos);
template void LatticeType::Cell_Pass(const long i0, const long i1,
				     ss_vect<tps> &x, long &lastpos);
