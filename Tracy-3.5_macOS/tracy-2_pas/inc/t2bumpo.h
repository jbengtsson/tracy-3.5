  { Interface }

  const	bpmbumpmax=10; bumpmax=200; bpmmax=200;

  type	bpmarray	= array [1..bpmbumpmax] of double;

	bumprec	= record
		    corrind	: array [1..3] of integer;
		    corrcoeff	: array [1..3] of double;
		    nbpm	: integer;
		    bpmind	: array [1..bpmbumpmax] of integer;
		    bpmcoeff	: bpmarray;
		    bpmcoeffsum	: double;
		  end;

	bumparray	= array [1..bumpmax] of bumprec;

{  var   Nbumps  : [external] ivector2;
        bumps   : [external] array [horizontal..vertical] of bumparray;
        bumpf   : [external] text;
        cod     : [external] array [1..bpmmax, horizontal..vertical] of double;
}

  procedure getbumprec(ncorr : integer; var corr : fitvect;
				plane, corr1, corr2, corr3 : integer;
				var bump : bumprec); external;

  procedure SetUpBump(ncorr : integer; var corr : fitvect;
			       dnumin : double; plane : integer); external;

  procedure setlocbump(plane : integer; maxkick : double;
		       theta : double; var bump : bumprec; lastpos : integer);
		      external;

  procedure execbump(MaxKick : double; lastpos : integer); external;

  procedure InitBUMP(var ncorr : ivector2; var hcorr, vcorr : fitvect;
			      dnuhmin, dnuvmin : double); external;

  procedure getbpm(Fnum, Knum : integer; var x, y : double;
		   lastpos : integer); external;
