/* Tracy-5

   J. Bengtsson, 2026.

*/


void Elem_Pass(const long i, tps &sigma)
{

  switch (Cell[i].Elem.Pkind) {
    case drift:
      Drift_Pass(Cell[i], sigma);
      break;
    case Mpole:
      Mpole_Pass(Cell[i], sigma);
      break;
    case Wigl:
      // Wiggler_Pass(Cell[i], sigma);
      break;
    case FieldMap:
      // FieldMap_Pass(Cell[i], sigma);
      break;
    case Insertion:
      // Insertion_Pass(Cell[i], sigma);
      break;
    case Cavity:
      Cav_Pass(Cell[i], sigma);
      break;
    case marker:
      Marker_Pass(Cell[i], sigma);
      break;
    case Spreader:
      break;
    case Recombiner:
      break;
    case Solenoid:
      // Solenoid_Pass(Cell[i], sigma);
      break;
    case Map:
      // Map_Pass(Cell[i], sigma);
      break;
    default:
      printf("Elem_Pass: *** undefined type - i = %ld type = %d\n",
	     i, Cell[i].Elem.Pkind);
      exit(1);
      break;
  }
}


void Cell_Pass(const long i0, const long i1, tps &sigma, long &lastpos)
{
  // Note: Sigma_k+1 = M_k Sigma_k M_k^T = (M_k (M_k Sigma_k)^T)^T
  const int n = 9;

  long int     jj[n][nv_tps];
  ss_vect<tps> Id, A;

  const double  deps = 1e-20;

#if 0

  Id.identity();

  map = Id + globval.CODvect;
  Cell_Pass(0, i0, map, lastpos);

  if (lastpos == i0) {
    map = Id + map.cst();
    Cell_Pass(i0, i1, map, lastpos);

    if (lastpos == i1) {
      // x_1 = zeta(x_0) => f_1(x) = f_0(zeta^-1(x))

      // Deterministic part.
      sigma = sigma*Inv(map-map.cst());
    }

  } else {
    printf("\nCell_Pass: particle lost at element %ld", lastpos);
    exit(1);
  }

#else

  for (auto i = i0; i <= i1; i++)
    Elem_Pass(i, sigma);
      
#endif

  if (globval.emittance) {
    // Stochastic part.

    for (auto i = 0; i < n; i++)
      for (auto j = 0; j < nv_tps; j++)
	jj[i][j] = 0;

    jj[0][x_]     = 2;
    jj[1][x_]     = 1;
    jj[1][px_]    = 1;
    jj[2][px_]    = 2;
    jj[3][y_]     = 2;
    jj[4][y_]     = 1;
    jj[4][py_]    = 1;
    jj[5][py_]    = 2;
    jj[6][ct_]    = 2;
    jj[7][ct_]    = 1;
    jj[7][delta_] = 1;
    jj[8][delta_] = 2;

    A = putlinmat(6, globval.Ascr);
    sigma = sigma*A;

    for (auto i = 0; i < 3; i++) {
      if (globval.eps[i] > deps) {
	sigma.pook(jj[3*i], sigma[jj[3*i]]-globval.D_rad[i]/2.0);
	sigma.pook(jj[3*i+2], sigma[jj[3*i+2]]-globval.D_rad[i]/2.0);
      }
    }

    sigma = sigma*Inv(A);
  }
}
