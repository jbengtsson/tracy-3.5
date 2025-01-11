module pascommon(input, output);

%include 'mathlib.def'
%include 'pascommono.def'

VAR	finame  : [global] str80; { input  data file  }
	foname  : [global] str80; { output data file }
	fname   : [global] str80; { temp file name }

	fi	: [global] text; { lattice input  file  }
	fo	: [global] text; { lattice output file }
	psin	: [global] array [0..maxincl] of text; { program input file }
	psout	: [global] text; { program output file}
	prr	: [global] array [3..maxfil] of text; { prr[1] : input, prr[2] : output }

	ErrFlag	: [global] boolean;

end.
