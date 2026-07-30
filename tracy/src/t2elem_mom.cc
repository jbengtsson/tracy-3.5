/* Tracy-5

   J. Bengtsson, 2026.

*/


void GtoL
(tps &ps, const Vector2 &S, const  Vector2 &R, const double c0, const double c1,
 const double s1)
{
  ss_vect<tps> M;

  M.identity();
  GtoL(M, S, R, c0, c1, s1);
  ps = ps*Inv(M);
}


void LtoG
(tps &ps, const Vector2 &S, const Vector2 &R, const double c0, const double c1,
 const double s1)
{
  ss_vect<tps> M;

  M.identity();
  LtoG(M, S, R, c0, c1, s1);
  ps = ps*Inv(M);
}


void Drift(const double L, tps &ps)
{
  ss_vect<tps> M;

  M.identity();
  Drift(L, M);
  ps = ps*Inv(M);
}

void Drift_Pass(const CellType &Cell, tps &ps)
{
  Drift(Cell.Elem.PL, ps);
}


void thin_kick
(CellType &Cell, const int Order, const double MB[], const double L,
 const double h_bend, const double h_ref, tps &ps)
{
  ss_vect<tps> M;

  M.identity();
  thin_kick(Cell, Order, MB, L, h_bend, h_ref, M);
  ps = ps*Inv(M);
}


void EdgeFocus(const double irho, const double phi, const double gap, tps &ps)
{
  ss_vect<tps> M;

  M.identity();
  EdgeFocus(irho, phi, gap, M);
  ps = ps*Inv(M);
}


void p_rot(const double phi, tps &ps)
{
  ss_vect<tps> M;

  M.identity();
  p_rot(phi, M);
  ps = ps*Inv(M);
}


void bend_fringe(const double hb, tps &ps)
{
  ss_vect<tps> M;

  M.identity();
  bend_fringe(hb, M);
  ps = ps*Inv(M);
}


void quad_fringe(const double b2, tps &ps)
{
  ss_vect<tps> M;

  M.identity();
  quad_fringe(b2, M);
  ps = ps*Inv(M);
}


void Mpole_Pass(CellType &Cell, tps &ps)
{
  const elemtype*  elemp = &Cell.Elem;
  const MpoleType* M     = elemp->M;

  int
    seg = 0, i;
  double
    dL = 0e0, dL1 = 0e0, dL2 = 0e0, dkL1 = 0e0, dkL2 = 0e0, h_ref = 0e0;

  GtoL(ps, Cell.dS, Cell.dT, M->Pc0, M->Pc1, M->Ps1);

  switch (M->Pmethod) {

  case Meth_Fourth:
    // Fringe fields.
    if (globval.quad_fringe && (M->PB[Quad+HOMmax] != 0e0))
      quad_fringe(M->PB[Quad+HOMmax], ps);
    if (!globval.Cart_Bend) {
      if (M->Pirho != 0e0) EdgeFocus(M->Pirho, M->PTx1, M->Pgap, ps);
    } else {
      p_rot(M->PTx1, ps);
      bend_fringe(M->Pirho, ps);
    }

    if (M->Pthick == thick) {
      if (!globval.Cart_Bend) {
	// Polar coordinates.
	h_ref = M->Pirho;
	dL = elemp->PL/M->PN;
      } else {
	// Cartesian coordinates.
	h_ref = 0e0;
	if (M->Pirho == 0e0)
	  dL = elemp->PL/M->PN;
	else
	  dL = 2e0/M->Pirho*sin(elemp->PL*M->Pirho/2e0)/M->PN;
      }

      dL1  = c_1*dL;
      dL2  = c_2*dL;
      dkL1 = d_1*dL;
      dkL2 = d_2*dL;

      for (seg = 1; seg <= M->PN; seg++) {
	Drift(dL1, ps);
	thin_kick(Cell, M->Porder, M->PB, dkL1, M->Pirho, h_ref, ps);
	Drift(dL2, ps);
	thin_kick(Cell, M->Porder, M->PB, dkL2, M->Pirho, h_ref, ps);

	Drift(dL2, ps);
	thin_kick(Cell, M->Porder, M->PB, dkL1, M->Pirho, h_ref, ps);
	Drift(dL1, ps);
      }
    } else
      thin_kick(Cell, M->Porder, M->PB, 1e0, 0e0, 0e0, ps);

    // Fringe fields.
    if (!globval.Cart_Bend) {
      if (M->Pirho != 0e0) EdgeFocus(M->Pirho, M->PTx2, M->Pgap, ps);
    } else {
      bend_fringe(-M->Pirho, ps);
      p_rot(M->PTx2, ps);
    }
    if (globval.quad_fringe && (M->PB[Quad+HOMmax] != 0e0))
      quad_fringe(-M->PB[Quad+HOMmax], ps);
    break;

  default:
    printf("Mpole_Pass: Method not supported %10s %d\n",
	   Cell.Elem.PName, M->Pmethod);
    exit_(0);
    break;
  }

  LtoG(ps, Cell.dS, Cell.dT, M->Pc0, M->Pc1, M->Ps1);
}


void Marker_Pass(CellType &Cell, tps &ps)
{
}


void Cav_Pass(const CellType &Cell, tps &ps)
{
  ss_vect<tps> M;

  M.identity();
  Cav_Pass(Cell, M);
  ps = ps*Inv(M);
}
