/* Tracy-2

   J. Bengtsson, CBP, LBL      1990 - 1994   Pascal version
                 SLS, PSI      1995 - 1997
   M. Boege      SLS, PSI      1998          C translation
   L. Nadolski   SOLEIL        2002          Link to NAFF, Radia field maps
   J. Bengtsson  NSLS-II, BNL  2004 -        

*/


template<typename T>
inline bool ElemType::CheckAmpl(ConfigType &conf, const ss_vect<T> &ps)
{
  bool not_lost;

  if (conf.Aperture_on)
    not_lost = is_double<T>::cst(ps[x_]) > maxampl[X_][0] &&
               is_double<T>::cst(ps[x_]) < maxampl[X_][1] && 
               fabs(is_double<T>::cst(ps[y_])) < maxampl[Y_][1];
  else
    not_lost = is_double<T>::cst(ps[x_]) > -max_ampl &&
               is_double<T>::cst(ps[x_]) < max_ampl &&
               fabs(is_double<T>::cst(ps[y_])) < max_ampl;

  if (!not_lost) {
    if (is_double<T>::cst(ps[x_]) < maxampl[X_][0] ||
        is_double<T>::cst(ps[x_]) > maxampl[X_][1])
      conf.lossplane = 1;
    else if (fabs(is_double<T>::cst(ps[y_])) > maxampl[Y_][1])
      conf.lossplane = 2;
  }

  return not_lost;
}


template<typename T>
void ElemType::Cell_Pass(ConfigType &conf, ss_vect<T> &ps)
{
  Elem_Pass(conf, ps);

  is_tps<T>::get_ps(ps, this);
}


template<typename T>
void LatticeType::Cell_Pass(const long i0, const long i1, ss_vect<T> &ps,
			    long &lastpos)
{
  long int i = 0;

  if (conf.radiation) conf.dE = 0e0;

  if (conf.emittance)
    for (i = 0; i < DOF; i++)
      conf.D_rad[i] = 0e0;

  if (!elems[i0]->CheckAmpl(conf, ps))
    lastpos = i0;
  else {
    lastpos = i1;
    for (i = i0; i <= i1; i++) {
      elems[i]->Cell_Pass(conf, ps);
      if (!elems[i]->CheckAmpl(conf, ps)) {
	if (trace)
	  printf("CheckAmpl: Particle lost in plane %d at element:"
		 " %5ld s = %10.5f, x = %12.5e, z= %12.5e\n",
		 conf.lossplane, i, elems[i]->S,
		 is_double<T>::cst(ps[x_]), is_double<T>::cst(ps[y_]));
	lastpos = i;
	break;
      }
    }
  }
}


void LatticeType::Cell_Pass(const long i0, const long i1, tps &sigma,
			    long &lastpos)
{
  // Note: Sigma_k+1 = M_k Sigma_k M_k^T = (M_k (M_k Sigma_k)^T)^T
  const int  n = 9;

  int           i, j;
  long int      jj[n][nv_tps];
  ss_vect<tps>  Id, A;

  const double  deps = 1e-20;

  Id.identity();

  map = Id + conf.CODvect; Cell_Pass(0, i0, map, lastpos);

  if (lastpos == i0) {
    map = Id + map.cst(); Cell_Pass(i0, i1, map, lastpos);

    if (lastpos == i1) {
      // x_1 = zeta(x_0) => f_1(x) = f_0(zeta^-1(x))

      // deterministic part
      sigma = sigma*Inv(map-map.cst());

      if (conf.emittance) {
	// stochastic part

	for (i = 0; i < n; i++)
	  for (j = 0; j < nv_tps; j++)
	    jj[i][j] = 0;

	jj[0][x_]  = 2; jj[1][x_]  = 1; jj[1][px_]    = 1; jj[2][px_]    = 2;
	jj[3][y_]  = 2; jj[4][y_]  = 1; jj[4][py_]    = 1; jj[5][py_]    = 2;
	jj[6][ct_] = 2; jj[7][ct_] = 1; jj[7][delta_] = 1; jj[8][delta_] = 2;

	A = putlinmat(6, conf.Ascr); sigma = sigma*A;

	for (i = 0; i < 3; i++) {
	  if (conf.eps[i] > deps) {
	    sigma.pook(jj[3*i], sigma[jj[3*i]]-conf.D_rad[i]/2.0);
	    sigma.pook(jj[3*i+2], sigma[jj[3*i+2]]-conf.D_rad[i]/2.0);
	  }
	}

	sigma = sigma*Inv(A);
      }
    }
  }
}


bool LatticeType::Cell_getCOD(long imax, double eps, double dP, long &lastpos)
{
  long            j, n, n_iter;
  int             no;
  long int        jj[ss_dim];
  double          dxabs;
  ss_vect<double> x0, x1, dx;
  ss_vect<tps>    I, dx0, map;

  no = no_tps; danot_(1);
  
  if (conf.mat_meth && (dP != conf.dPparticle))
    // Recompute transport matrices.
    get_lin_maps(dP);

  conf.dPparticle = dP;

  n = (conf.Cavity_on)? 6 : 4;

  x0.zero(); x0[delta_] = dP;

  // if (n == 4) {
  //   // For 2 1/2 D.O.F.: eta*dP. 
  //   x0[x_] = elems[0]->Eta[X_]*dP; x0[px_] = elems[0]->Etap[X_]*dP;
  //   x0[y_] = elems[0]->Eta[Y_]*dP; x0[py_] = elems[0]->Etap[Y_]*dP;
  // }

  for (j = 0; j < ss_dim; j++)
    jj[j] = (j < n)? 1 : 0;

  if (trace) {
    std::cout << "\nCell_getCOD:" << "\n";
    std::cout << std::scientific << std::setprecision(5)
	      << "  0                        x0 ="
	      << std::setw(13) << x0 << "\n";
  }
  n_iter = 0; I.identity();
  do {
    n_iter++; map.identity(); map += x0;

    Cell_Pass(0, conf.Cell_nLoc, map, lastpos); 

    if (lastpos == conf.Cell_nLoc) {
      x1 = map.cst(); dx = x0 - x1; dx0 = PInv(map-I-x1, jj)*dx;
      dxabs = xabs(n, dx); x0 += dx0.cst();
    } else {
      dxabs = NAN; break;
    }

    if (trace)
      std::cout
	<< std::scientific << std::setprecision(1)
	<< std::setw(3) << n_iter
	<< " err = " << std::setw(7) << dxabs << "/" << std::setw(7) << eps
	<< std::setprecision(5)	<< "  x0 =" << std::setw(13) << x0 << "\n";
  } while ((dxabs >= eps) && (n_iter <= imax));

  conf.codflag = dxabs < eps;

  if (conf.codflag) {
    conf.CODvect = x0; getlinmat(6, map, conf.OneTurnMat);
    Cell_Pass(0, conf.Cell_nLoc, x0, lastpos);
  } else
    std::cout << std::scientific << std::setprecision(5)
	      << "\nCell_getCOD: failed to converge after " << n_iter
	      << " iterations:\n"
	      << "  dP =" << std::setw(12) << dP
	      << ", particle lost at element " << lastpos << "\n"
	      << std::scientific << std::setprecision(5)
	      << "  x_0   =" << std::setw(13) << x0 << "\n"
	      << std::scientific << std::setprecision(5)
	      << "  x_k-1 =" << std::setw(13) << elems[lastpos-1]->BeamPos
	      << "\n"
	      << std::scientific << std::setprecision(5)
	      << "  x_k   =" << std::setw(13) << map.cst() << "\n";

  danot_(no);
  
  return conf.codflag;
}


bool LatticeType::GetCOD(long imax, double eps, double dP, long &lastpos)
{
  return Cell_getCOD(imax, eps, dP, lastpos);
}


bool LatticeType::getcod(double dP, long &lastpos)
{
  return GetCOD(conf.CODimax, conf.CODeps, dP, lastpos);
}


void LatticeType::Lat_Init(void)
{
  long int i;
  double   Stotal;

  char first_name[] = "begin          ";

  if (debug) printf("Lat_Init: Cell_nLoc = %ld\n", conf.Cell_nLoc);

  SI_init();  /* Initializes the constants for symplectic integrator */

  // Allocate space for lattice.
  elems.resize(conf.Cell_nLoc+1);

  // Assign element 0 ("begin").
  elems[0] = Marker_Alloc();
  memcpy(elems[0]->PName, first_name, sizeof(first_name));
  elems[0]->PL = 0e0; elems[0]->Fnum = 0; elems[0]->Knum = 0;
  elems[0]->Pkind = marker;
  elems[0]->dT[X_] = 1e0; elems[0]->dT[Y_] = 0e0;
  elems[0]->dS[X_] = 0e0; elems[0]->dS[Y_] = 0e0;

  for (i = 1; i <= conf.Elem_nFam; i++) {
    // Allocate element.
    elemf[i-1].ElemF->Elem_Init(*this, i);
    if (debug)
      printf("\nLat_Init i = %3ld %*s\n",
	     i, SymbolLength, elemf[i-1].ElemF->PName);
  }

  /* Computes s-location of each element in the structure */
  Stotal = 0e0;
  for (i = 0; i <= conf.Cell_nLoc; i++) {
    Stotal += elems[i]->PL; elems[i]->S = Stotal;
  }

  if (debug) prt_elem();
}
