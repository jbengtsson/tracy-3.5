// Insertion-device linear-optics correction — see correction/id_corr.h.

// Defined in nsls-ii_lib, which is not included here.
void set_ID_scl(const int Fnum, const double scl);

namespace corr {

double Bet(double bq, double nus, double nuq, double NuQ)
{
  return bq*cos(2.0*M_PI*(2.0*fabs(nus-nuq)-NuQ))/(2.0*sin(2.0*M_PI*NuQ));
}


double Nus(double bq, double nus, double nuq, double NuQ)
{
  double Nu, sgn;

  sgn = ((nus-nuq) <= 0)? -1: 1;

  Nu = -bq*sgn*(sin(2.0*M_PI*NuQ)+sin(2.0*M_PI*(2.0*fabs(nus-nuq)-NuQ)))
       /(8.0*M_PI*sin(2.0*M_PI*NuQ));

  return Nu;
}


id_corr::~id_corr()
{
  // Allocation happens once in ini_ID_corr; a non-null A1 marks it. Freeing at
  // teardown rather than at the end of ID_corr is what allows ID_corr to be
  // called repeatedly against one allocation.
  if (A1 != 0) {
    free_dvector(Xsext, 1, Nconstr); free_dvector(Xsext0, 1, Nconstr);
    free_dvector(b2Ls_, 1, Nquad); free_dmatrix(A1, 1, Nconstr, 1, Nquad);
    free_dmatrix(U1, 1, Nconstr, 1, Nquad); free_dvector(w1, 1, Nquad);
    free_dmatrix(V1, 1, Nquad, 1, Nquad);
  }
}


void id_corr::get_IDs(void)
{
  int k;

  printf("\n");
  n_ID_Fams = 0;
  for (k = 0; k < globval.Elem_nFam; k++)
    switch (ElemFam[k].ElemF.Pkind) {
    case Wigl:
      printf("found ID family:   %s %12.5e\n",
	     ElemFam[k].ElemF.PName, ElemFam[k].ElemF.W->BoBrhoV[0]);
      n_ID_Fams++; ID_Fams[n_ID_Fams-1] = k + 1;
      break;
    case Insertion:
      printf("found ID family:   %s %12.5e",
	     ElemFam[k].ElemF.PName, ElemFam[k].ElemF.ID->scaling);
      if (ElemFam[k].ElemF.ID->scaling != 0e0) {
	printf("\n");
	n_ID_Fams++; ID_Fams[n_ID_Fams-1] = k + 1;
      } else
	printf("  not included\n");
      break;
    case FieldMap:
      printf("found ID family:   %s %12.5e\n",
	     ElemFam[k].ElemF.PName, ElemFam[k].ElemF.FM->scl);
      n_ID_Fams++; ID_Fams[n_ID_Fams-1] = k + 1;
      break;
    default:
      break;
    }
}


void id_corr::set_IDs(const double scl)
{
  int k;

  printf("\n");
  for (k = 0; k < n_ID_Fams; k++) {
    switch (ElemFam[ID_Fams[k]-1].ElemF.Pkind) {
    case Wigl:
      printf("setting ID family: %s %12.5e\n",
	     ElemFam[ID_Fams[k]-1].ElemF.PName,
	     scl*ElemFam[ID_Fams[k]-1].ElemF.W->BoBrhoV[0]);

      set_Wiggler_BoBrho(ID_Fams[k],
			 scl*ElemFam[ID_Fams[k]-1].ElemF.W->BoBrhoV[0]);
      break;
    case Insertion:
      printf("setting ID family: %s %12.5e\n",
	     ElemFam[ID_Fams[k]-1].ElemF.PName, scl);

      set_ID_scl(ID_Fams[k], scl);
      break;
    case FieldMap:
      printf("setting ID family: %s %12.5e\n",
	     ElemFam[ID_Fams[k]-1].ElemF.PName, scl);

      set_ID_scl(ID_Fams[k], scl);
      break;
    default:
      std::cout << "set_IDs: unknown element type" << std::endl;
      exit_(1);
      break;
    }
  }
}


// Restore each quad family's b_2 to the design snapshot quad_config captured.
// TODO: rename — this is not a general quad reset, and the name reads like one.
void id_corr::reset_quads(const int N_Fam, const int Q_Fam[])
{
  int k;

  if (N_Fam > N_Fam_max) {
    printf("reset_quads: N_Fam > N_Fam_max: %d (%d)\n", N_Fam, N_Fam_max);
    exit_(0);
  }

  for (k = 0; k < N_Fam; k++) {
    // Note, actual values can differ from the original values
/*    printf("setting quad family: %s %12.5e\n",
	   ElemFam[Q_Fam[k]-1].ElemF.PName,
	   ElemFam[Q_Fam[k]-1].ElemF.M->PBpar[HOMmax+Quad]);

    set_bn_design_fam(Q_Fam[k], Quad,
		       ElemFam[Q_Fam[k]-1].ElemF.M->PBpar[HOMmax+Quad], 0.0);*/

    printf("setting quad family: %s %12.5e\n",
	   ElemFam[Q_Fam[k]-1].ElemF.PName, b2[k]);

    set_bn_design_fam(Q_Fam[k], Quad, b2[k], 0.0);
  }
}


void id_corr::SVD(const int m, const int n, double **M,
		  double beta_nu[], double b2Ls_[], const bool first,
		  const double ID_s_cut)
{
  if (trace) {
    printf("\n");
    printf("SVD: first = %1d, m = %1d n = %1d\n", first, m, n);
  }

  if (first) corr::svd_decomp_cut(M, m, n, U1, w1, V1, ID_s_cut, true);

  corr::svd_backsub(U1, w1, V1, m, n, beta_nu, b2Ls_);
}


void id_corr::quad_config(const int N_Fam, const int Q_Fam[])
// Collect quadrupole trims. For the block diagonal. Linear optics response matrix (b_2).
// Quadrupole families get picked up from "ID_quads" in the parameter file.
{
  int    i, j;
  double an;

  if (N_Fam > N_Fam_max) {
    printf("quad_config: N_Fam > N_Fam_max: %d (%d)\n", N_Fam, N_Fam_max);
    exit_(0);
  }

  Nquad = 0;
  for (i = 0; i < N_Fam; i++) {
    for (j = 1; j <= GetnKid(Q_Fam[i]); j++) {
      Nquad++;

      if (Nquad > n_b2_max) {
        printf("quad_config: max no of quadrupoles exceeded %d (%d)\n",
               Nquad, n_b2_max);
        exit_(1);
      }

      quad_prms[Nquad-1] = Elem_GetPos(Q_Fam[i], j);

      if (j == 1) get_bn_design_elem(Q_Fam[i], j, Quad, b2[i], an);
    }
  }

  printf("\n");
  printf("quad_config: Nquad = %d\n", Nquad);
}


bool id_corr::get_SQ(void)
{
  int  j, k;
//  Vector2  alpha3[3], beta3[3], nu3[3], eta3[3], etap3[3];
  FILE *outf = NULL;

  /* Note, IDs are split for evaluation of the driving terms at the center:
       id1  1, 2
       id2  1, 2
       ...                                                                  */

  // Get Twiss params, no dispersion
  Ring_GetTwiss(false, 0e0);

  if (!status.codflag || !globval.stable) return false;

  // Get global tunes
  Nu_X = globval.TotalTune[X_]; Nu_Y = globval.TotalTune[Y_];

  if (trace) {
    printf("\n");
    printf("nu_x = %8.12f, nu_y = %8.12f\n", Nu_X, Nu_Y);

    // Get Twiss params in sext
    printf("\n");
    printf("Lattice functions at sextupoles:\n");

    outf = file_write("latfunS.out");

    fprintf(outf, "s betax nux betay nuy\n");
  }

  Nsext = 0;
  for (k = 0; k < globval.Cell_nLoc; k++) {
    if ((Cell[k].Elem.Pkind == Mpole) && (Cell[k].Elem.M->n_design == Sext)) {
      Nsext++;

      if (Nsext > n_b3_max) {
        printf("get_SQ: max no of sextupoles exceeded %d (%d)\n",
               Nsext, n_b3_max);
        exit_(1);
      }

      Ss[Nsext-1] = Cell[k].S; S_locs[Nsext-1] = k;

      for (j = 0; j <= 1; j++) {
	sb[j][Nsext-1] = Cell[k].Beta[j];
	sNu[j][Nsext-1] = Cell[k].Nu[j] - nu_0[j];
      }

      if (trace) {
	printf("%-8s %7.3f %8.5f %8.5f %8.5f %8.5f\n",
	       Cell[k].Elem.PName, Ss[Nsext-1],
	       sb[X_][Nsext-1], sNu[X_][Nsext-1]-nu_0[X_],
	       sb[Y_][Nsext-1], sNu[Y_][Nsext-1]-nu_0[Y_]);
	fprintf(outf, "%-8s %7.3f %8.5f %8.5f %8.5f %8.5f\n",
		Cell[k].Elem.PName, Ss[Nsext-1],
		sb[X_][Nsext-1], sNu[X_][Nsext-1]-nu_0[X_],
		sb[Y_][Nsext-1], sNu[Y_][Nsext-1]-nu_0[Y_]);
      }
    }
  }

  if (trace) fclose(outf);

  // Number of sexts in the ring
  printf("No of sextupoles = %d\n", Nsext);

  if (trace) {
    // Get Twiss params in quads
    printf("\n");
    printf("Lattice functions at quadrupoles:\n");

    outf = file_write("latfunQ.out");

    fprintf(outf, "s name betax nux betay nuy\n");
  }

  for (k = 0; k < Nquad; k++) {
    Sq[k] = Cell[quad_prms[k]].S;
    for (j = 0; j <= 1; j++) {
      // does not work for machine file (get_twiss_3)...
//       if (Cell[quad_prms[k]].Elem.M->Pthick == thick) {
// 	get_twiss3(quad_prms[k], alpha3, beta3, nu3, eta3, etap3);
// 	qb[j][k] = beta3[Y_][j]; qNu[j][k] = nu3[Y_][j] - nu_0[j];
//       } else {
	qb[j][k] = Cell[quad_prms[k]].Beta[j];
	qNu[j][k] = Cell[quad_prms[k]].Nu[j] - nu_0[j];
//       }
    }

    if (trace) {
      printf("%-8s %7.3f %8.5f %8.5f %8.5f %8.5f\n",
	     Cell[quad_prms[k]].Elem.PName, Sq[k], qb[X_][k],
	     qNu[X_][k], qb[Y_][k], qNu[Y_][k]);

      fprintf(outf, "%-8s %7.3f %8.5f %8.5f %8.5f %8.5f\n",
	      Cell[quad_prms[k]].Elem.PName, Sq[k], qb[X_][k],
	      qNu[X_][k], qb[Y_][k], qNu[Y_][k]);
    }
  }

  if (trace) fclose(outf);

  // Number of quads in the ring
  printf("No of quads      = %d\n", Nquad);

  return true;
}


void id_corr::A_matrix(void)
{
  const string file_name = "AA.dat";

  int    k, j;
  double BtX, BtY, NuX, NuY;
  FILE   *outf;

  // Defining Twiss in undisturbed quads
  for (k = 0; k < Nquad; k++)
    for (j = 0; j <= 1; j++) {
      qb0[j][k] = qb[j][k]; qNu0[j][k] = qNu[j][k];
    }

  // Defining Twiss in undisturbed sexts
  for (k = 0; k < Nsext; k++)
    for (j = 0; j <= 1; j++)
      sNu0[j][k] = sNu[j][k];

  // Now creating matrix A in X=A*B2L
  for (k = 1; k <= Nsext; k++) {
    for (j = 1; j <= Nquad; j++) {
      BtX = Bet(qb0[X_][j-1], sNu0[X_][k-1], qNu0[X_][j-1], Nu_X0);
      NuX = -Nus(qb0[X_][j-1], sNu0[X_][k-1], qNu0[X_][j-1], Nu_X0);
      BtY = -Bet(qb0[Y_][j-1], sNu0[Y_][k-1], qNu0[Y_][j-1], Nu_Y0);
      NuY = Nus(qb0[Y_][j-1], sNu0[Y_][k-1], qNu0[Y_][j-1], Nu_Y0);
      A1[k][j] = scl_dbeta*BtX;
      A1[k+Nsext][j] = scl_dbeta*BtY;
      A1[k+2*Nsext][j] = scl_dnu*NuX;
      A1[k+3*Nsext][j] = scl_dnu*NuY;
    }
  }
  // Now adding 2 more constraints for global tunes
  for (j = 1; j <= Nquad; j++) {
    A1[4*Nsext+1][j] = -scl_nu*qb0[X_][j-1]/(4.0*M_PI);
    A1[4*Nsext+2][j] =  scl_nu*qb0[Y_][j-1]/(4.0*M_PI);
  }

  if (trace) {
    outf = file_write(file_name.c_str());
    fprintf(outf, "\n");
    fprintf(outf, "AA:\n");
    fprintf(outf, "\n");
    for (k = 1; k <= Nconstr; k++) {
      for (j = 1; j <= Nquad; j++)
	fprintf(outf, " %10.3e", A1[k][j]);
      fprintf(outf, "\n");
    }
    fclose(outf);
  }
}


void id_corr::X_vector(const bool first)
// Linear optics distortion vector [\delta\beta,\delta\mu,\delta\nu]
{
  int k;

  dnu0[X_] = globval.TotalTune[X_] - Nu_X0;
  dnu0[Y_] = globval.TotalTune[Y_] - Nu_Y0;

  if (first) {
    // Initial fill of X
    for (k = 1; k <= Nsext; k++) {
      Xsext0[k]         = sb[X_][k-1];  Xsext0[k+Nsext]   = sb[Y_][k-1];
      Xsext0[k+2*Nsext] = sNu[X_][k-1]; Xsext0[k+3*Nsext] = sNu[Y_][k-1];
    }
    Xsext0[4*Nsext+1] = 0.0; Xsext0[4*Nsext+2] = 0.0;
  } else {
    // Now substracting from X in X=A*B2L
    for (k = 1; k <= Nsext; k++) {
      Xsext[k]         = scl_dbeta*(Xsext0[k]-sb[X_][k-1])/sb[X_][k-1];
      Xsext[k+Nsext]   = scl_dbeta*(Xsext0[k+Nsext]-sb[Y_][k-1])/sb[Y_][k-1];
      Xsext[k+2*Nsext] = scl_dnu*(Xsext0[k+2*Nsext]-sNu[X_][k-1]+dnu0[X_]/2.0);
      Xsext[k+3*Nsext] = scl_dnu*(Xsext0[k+3*Nsext]-sNu[Y_][k-1]+dnu0[Y_]/2.0);
    }
    Xsext[4*Nsext+1] = scl_nu*(Nu_X0-globval.TotalTune[X_]);
    Xsext[4*Nsext+2] = scl_nu*(Nu_Y0-globval.TotalTune[Y_]);
  }

  if (trace) {
    printf("\n");
    printf("X:\n");
    printf("\n");
    if (first) {
      for (k = 1; k <= Nconstr; k++) {
	printf(" %10.3e", Xsext0[k]);
	if (k % 10 == 0)  printf("\n");
      }
      if (Nconstr % 10 != 0) printf("\n");
    } else {
      for (k = 1; k <= Nconstr; k++) {
	printf(" %10.3e", Xsext[k]);
	if (k % 10 == 0)  printf("\n");
      }
      if (Nconstr % 10 != 0) printf("\n");
    }
  }
}

// Initializing ID correction (NOT LOCO)
void id_corr::ini_ID_corr(const bool IDs, const int N_Fam, const int Q_Fam[])
{
  int k;

  if (IDs) {
    // store ID families
    get_IDs();

    // zero ID's
    set_IDs(0.0);
  }

  // Configuring quads (1 --> C means thin quads located in the middle of 1s)
  quad_config(N_Fam, Q_Fam);

  // Read Betas and Nus
  get_SQ();
  Nconstr = 4*Nsext + 2;

  // Note, allocated vectors and matrices are deallocated in the destructor.
  Xsext = dvector(1, Nconstr); Xsext0 = dvector(1, Nconstr);
  b2Ls_ = dvector(1, Nquad); A1 = dmatrix(1, Nconstr, 1, Nquad);
  U1 = dmatrix(1, Nconstr, 1, Nquad); w1 = dvector(1, Nquad);
  V1 = dmatrix(1, Nquad, 1, Nquad);

  for (k = 1; k <= Nquad; k++)
    b2Ls_[k] = 0.0;

  // shift zero point to center of ID
  //  nu_0[X_] = Cell[id_loc].Nu[X_]; nu_0[Y_] = Cell[id_loc].Nu[Y_];
  nu_0[X_] = 0.0;
  nu_0[Y_] = 0.0;

  // Defining undisturbed tunes
  Nu_X0 = globval.TotalTune[X_]; Nu_Y0 = globval.TotalTune[Y_];

  // Set-up matrix A in X=A*b2Ls_
  A_matrix();

  // Now fill the X in X=A*b2Ls_
  X_vector(true);
}


void id_corr::W_diag(void)
{
  double bxf, byf, nxf, nyf, b2Lsum;
  int    k;

  bxf = 0.0; byf = 0.0; nxf = 0.0; nyf = 0.0;
  for (k = 1; k <= Nsext; k++) {
    bxf += sqr(Xsext[k]);
    byf += sqr(Xsext[k+Nsext]);
    nxf += sqr(Xsext[k+2*Nsext]);
    nyf += sqr(Xsext[k+3*Nsext]);
  }

  dnu0[X_] = globval.TotalTune[X_] - Nu_X0;
  dnu0[Y_] = globval.TotalTune[Y_] - Nu_Y0;

  b2Lsum = 0.0;
  for (k = 1; k <= Nquad; k++)
    b2Lsum += sqr(b2Ls_[k]);

  printf("\n");
  printf("Residuals: beta [%%], dnu : \n");
  printf("dbeta_x: %6.2f dbeta_y: %6.2f nu_x: %12.6e nu_y: %12.6e\n",
	 sqrt(bxf/Nsext)*1e2, sqrt(byf/Nsext)*1e2,
	 sqrt(nxf/Nsext), sqrt(nyf/Nsext));
  printf("tune shift: dnu_x = %8.5f, dnu_y = %8.5f\n", dnu0[X_], dnu0[Y_]);
  printf("Sum b2Ls_: %12.6e\n", sqrt(b2Lsum/Nquad));
}


bool id_corr::ID_corr(const int N_calls, const int N_steps,
		      const bool IDs, const int cnt, const int N_Fam,
		      const int Q_Fam[], const double ID_s_cut)
{
  int    i, j, k, Fnum;
  double b2L, a2L, L;
  FILE   *outf;
  char fname[30];

  a2L=b2L=L=0.;

  printf("\n");
  printf("ID matching begins!\n");


  snprintf(fname, sizeof(fname), "ID_corr_%d.out",cnt);
  outf = file_write(fname);

  for (i = 1; i <= N_steps; i++) { //This brings ID strength in steps
    if (IDs) set_IDs((double)i/(double)N_steps);

    get_SQ();                               // Read Betas and Nus
    X_vector(false);                        // Fill in dX in dX=A*db2Ls_
    W_diag();                               // Get statistics
    for (j = 1; j <= N_calls; j++) {
      SVD(Nconstr, Nquad, A1, Xsext, b2Ls_, j == 1, ID_s_cut);

      if ((i == N_steps) && (j == N_calls)) fprintf(outf, "#b_2:\n");

      // add quad strengths (db2Ls_)
      for (k = 1; k <= Nquad; k++) {
	set_dbnL_design_elem(Cell[quad_prms[k-1]].Fnum,
			     Cell[quad_prms[k-1]].Knum, Quad,
			     -ID_step*b2Ls_[k], 0.0);

	if ((i == N_steps) && (j == N_calls)) {
	  Fnum = Cell[quad_prms[k-1]].Fnum; L = Cell[quad_prms[k-1]].Elem.PL;
	  get_bnL_design_elem(Fnum, Cell[quad_prms[k-1]].Knum, Quad, b2L, a2L);
	  // ElemFam not defined for flat file.
	  // fprintf(outf, "%10s %6.2f %3d %8.5f\n",
	  // 	  Cell[quad_prms[k-1]].Elem.PName, Cell[quad_prms[k-1]].S, k,
	  // 	  b2L-ElemFam[Fnum-1].ElemF.M->PBpar[HOMmax+Quad]*L);
	  fprintf(outf, "%10s %6.2f %3d %8.5f\n",
		  Cell[quad_prms[k-1]].Elem.PName, Cell[quad_prms[k-1]].S, k,
		  b2L);
	}
      }

      printf("\n");
      printf("Iteration: %2d\n", j);
      if (get_SQ()) {
	X_vector(false);                    // Fill in dX in dX=A*db2Ls_
	W_diag();                           // Get statistics

	printglob();
      } else {
	printf("ID_corr: correction failed\n");
	// restore lattice
	if (IDs) set_IDs(0.0);
	reset_quads(N_Fam, Q_Fam);
	return false;
      }
    }
  }
  fclose(outf);

  snprintf(fname, sizeof(fname), "ID_corr_res_%d.out",cnt);
  outf = file_write(fname);

  fprintf(outf, "# dbeta_x/beta_x  dbeta_y/beta_y  dnu_x  dnu_y\n");
  fprintf(outf, "#      [%%]             [%%]\n");
  fprintf(outf, "#\n");
  for (k = 1; k <= Nsext; k++)
    fprintf(outf, "%6.1f %6.2f %6.2f %10.3e %10.3e\n",
	    Ss[k-1], 1e2*Xsext[k], 1e2*Xsext[k+Nsext],
	    Xsext[k+2*Nsext], Xsext[k+3*Nsext]);
  fclose(outf);

  printf("\n");
  printf("ID matching ends!\n");

  return true;
}

}  // namespace corr
