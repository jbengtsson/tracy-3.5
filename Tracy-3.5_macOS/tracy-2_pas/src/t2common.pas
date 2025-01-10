  module t2common(input, output);

  %include 'mathlib.def'
  %include 't2commono.def'

  var	globval			: [global] globvalrec;
  	status			: [global] statusrec;
	trace, cellconcat	: [global] boolean;

  end.
