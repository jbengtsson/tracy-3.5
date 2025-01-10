	{ GRAPHX routines,  external binding, implies passing addresses }

	procedure graope; external;

	procedure graph0(var x, y : graphvect; %ref nel : integer); external;

	procedure graph1(var x, y : graphvect; %ref nel, linas, imar : integer;
			 %ref hmar : real); external;

	procedure gracls; external;

	procedure graclx(%ref ifdmp, ifptr, ifmet, ifkee : integer); external;

	procedure graheds(%descr title : strbuf; %ref len : integer;
			  %ref th : real); external;

	procedure gravwp(%ref gx0, gy0, gxlen, gylen : real); external;

	procedure grawnd(%ref gx0, gy0, gxlen, gylen : real); external;

	procedure gracoms(%ref h : real;
			  %descr xcom : strbuf; %ref xlen, nxpos : integer;
			  %descr ycom : strbuf; %ref ylen, nyelm : integer
			 ); external;

	procedure graxi1(%ref ifxax, nxdc, nxxp : integer; %ref pyaxp : real);
			external;

	procedure grayi1(%ref ifyax, nydc, nyxp : integer; %ref pxaxp : real);
			external;

	procedure graxi2(%ref pxdel : real; %ref mxdiv, nxcen : integer);
			 external;

	procedure grayi2(%ref pydel : real; %ref mydiv, nycen : integer);
			external;

	procedure grabeg; external;

	procedure graend; external;

	procedure gravcs(%ref gxmi, gymi, gxma, gyma : real); external;

	procedure gravcf(%ref ifklip : integer); external;

	procedure gravec(%ref gx1, gy1, gx2, gy2 : real); external;

	procedure grafpl(var x, y : graphvect; %ref nel : integer); external;

	procedure gracoi(%ref icol : integer); external;

	procedure gratyp(%ref ityp : integer); external;

	procedure grapat(%ref ipat : integer); external;

	procedure gratexs(%ref xst, yst, zh : real; %descr string : strbuf;
			  %ref strlen : integer; %ref win : real;
			  var xe, ye : real; %ref jfpl : integer); external;

	procedure gramar(%ref gx, gy : real; %ref imar : integer;
			 %ref ghmar : real); external;
