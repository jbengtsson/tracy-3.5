{ Interface }

const	NameLength = 15;   { maximum length of identifiers }
	blankname = '               '; maxincl = 5; maxfil = 10;
	bigvectmax = 4096;

type
        str255	= packed array [1..255] of char;
        str127	= packed array [1..127] of char;
        str80	= packed array [1..80] of char;
        str40	= packed array [1..40] of char;

        alfa	= packed array  [1..NameLength] of char;

	bigvect	= array [1..bigvectmax] of double;

VAR	finame  : [external] str80; { input  data file  }
	foname  : [external] str80; { output data file }
	fname   : [external] str80; { temp file name }

	fi	: [external] text; { lattice input  file  }
	fo	: [external] text; { lattice output file }
	psin	: [external] array [0..maxincl] of text; { program input file }
	psout	: [external] text; { program output file}
	prr	: [external] array [3..maxfil] of text; { prr[1] : input, prr[2] : output }

	ErrFlag	: [external] boolean;
