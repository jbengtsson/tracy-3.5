  procedure fcoeff(j, k, l, m, p : integer; fact : double; var cr, sr : double);

  { Compute coefficients in the Fourier expanded Hamiltonian }

  var	i, mx, my, nx, ny, n				: integer;
	b3, betaxi, betayi, nuxi, nuyi, nux, nuy	: double;
	C, si						: double;
	cell						: celltype;

  function nok(n, k : integer) : integer;

  var	i	: integer;
	x	: double;

  begin
    if k = 0 then
      nok := 1
    else
    begin
      x := n/k;
      for i:=1 to k-1 do
        x := x*(n-i)/i;
      nok := round(x);
    end;
  end;

  begin
    mx := j + k; nx := j - k; my := l + m; ny := l - m;
    n := j + k + l + m;
    nux := globval.totaltune[1]; nuy := globval.totaltune[2];
    getelem(globval.cell_nloc, cell); C := cell.s;
    cr := 0.0; sr := 0.0;
    for i:=1 to globval.cell_nloc do
    begin
      getelem(i, cell);
      if (cell.elem.Pkind = mpole) and (cell.elem.M.Porder >= sext) then
      begin
        b3 := cell.elem.M.Pb[sext];
        nuxi := cell.nu[1]; nuyi := cell.nu[2];
        betaxi := cell.beta[1]; betayi := cell.beta[2];
        si := cell.s;
	cr := cr + b3*pwr(betaxi, mx/2)*pwr(betayi, my/2)
	      *cos(2.0*pi*(nx*(nuxi-nux*si/C)+ny*(nuyi-nuy*si/C)+p*si/C));
	sr := sr + b3*pwr(betaxi, mx/2)*pwr(betayi, my/2)
	      *sin(2.0*pi*(nx*(nuxi-nux*si/C)+ny*(nuyi-nuy*si/C)+p*si/C));
      end;
    end;
    cr := nok(mx, j)*nok(my, l)/(2*pi*pwr(2, n/2))*cr*fact;
    sr := nok(mx, j)*nok(my, l)/(2*pi*pwr(2, n/2))*sr*fact;
  end;
