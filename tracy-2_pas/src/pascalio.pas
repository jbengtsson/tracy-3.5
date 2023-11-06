  module pascalio(input, output);

  %include 'mathlib.def'
  %include 'dab.def'
  %include 'pascommon.def'
  %include 't2common.def'
  %include 'pascalio.def'

  { Interface }

  [global] procedure reset_(var fvar : text;
                   fname : packed array [low..high : integer] of char); forward;

  [global] procedure rewrite_(var fvar : text;
                     fname : packed array [low..high : integer] of char);
		     forward;


  [global] procedure t2init; forward;

  [global] procedure getglobv_(var globval1 : globvalrec); forward;

  [global] procedure putglobv_(var globval1 : globvalrec); forward;

  { Implementation }

  procedure reset_{var fvar : text;
                   fname : packed array [low..high : integer] of char};

  begin
    open(fvar, fname, history:=old); reset(fvar);
  end;


  procedure rewrite_{var fvar : text;
                     fname : packed array [low..high : integer] of char};

  begin
    open(fvar, fname, history:=new); rewrite(fvar);
  end;


  procedure getglobv_{var globval1 : globvalrec};

  begin
    globval1 := globval;
  end;


  procedure putglobv_{var globval1 : globvalrec};

  begin
    globval := globval1;
  end;


  procedure t2init;

  begin
    iniranf(0);
    pi := 4d0*arctan(1d0); { VAX }
{    pi := 4d0*arctan(1e0); SUN }
    cellconcat := false;
    DAini(no, nv, 0); Lieini(no, nv, nd2);
  end;

  end.
