  { Interface }

  procedure ETY(n, low, igh : integer; var a : matrix; var ort : vector); external;

  procedure ETYT(n, low, igh : integer; var a : matrix; var ort : vector;
		 var z : matrix); external;

  procedure ety2(n, low, igh : integer; var h : matrix; var wr, wi : vector;
		 var z : matrix; var ierr : integer); external;
