  { Routine for Pascal-S version of input file }

  procedure init;

  begin
    pi := 4*arctan(1.0);
  end;

  function sngl(dble : double) : double;

  { Dummy routine for conversion from double to single precision }

  begin
    sngl := dble;
  end;
