  const	simpvectmax = 100;

  type	simpvect	= array [1..simpvectmax] of double;

  procedure amoeba(var x, dx : simpvect; %ref ndim : integer; var fx : double;
		   %ref ftol : double; %ref niter : integer;
		   var ipos : integer); external;
