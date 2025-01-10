  module b2perp(input, output);

  %include 'mathlib.def'
  %include 't2common.def'
  %include 'B2perp.def'

  [global] function B2perp(var irho : double; var B : vector3;
			   var x, xp, yp : double) : double;

  {             -   - 2        -
    Calculates |B x e| , where e is a unit vector in the direction of
    propagation } 
    
  var	xn	: double;
	e	: vector3;
      
  begin
    xn := 1d0/sqrt(sqr(1d0+x*irho)+sqr(xp)+sqr(yp));
    e[1] := xp*xn; e[2] := yp*xn; e[3] := (1d0+x*irho)*xn;
    B2perp := sqr(B[2]*e[3]-B[3]*e[2]) + sqr(B[1]*e[2]-B[2]*e[1])
	    + sqr(B[3]*e[1]-B[1]*e[3]);
  end;

  end.
