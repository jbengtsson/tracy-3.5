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

{ VAR	finame  : [external] str80;
	foname  : [external] str80;
	fname   : [external] str80;

	fi	: [external] text;
	fo	: [external] text;
	psin	: [external] array [0..maxincl] of text;
	psout	: [external] text;
	prr	: [external] array [3..maxfil] of text;

	ErrFlag	: [external] boolean;
}
