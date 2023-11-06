module svbksbmod(input, output);

{ Interface }

%include 'svd.def'
%include 'svbksb.def'

[global] PROCEDURE svbksb(u: svdmat; w: svdarray; v: svdmat;
			  m, n: integer; b: svdarray; VAR x: svdarray);
VAR
   jj,j,i: integer;
   s: double;
   tmp: svdarray;
BEGIN
   FOR j := 1 to n DO BEGIN
      s := 0d0;
      IF (w[j] <> 0d0) THEN BEGIN
         FOR i := 1 to m DO BEGIN
            s := s+u[i,j]*b[i]
         END;
         s := s/w[j]
      END;
      tmp[j] := s
   END;
   FOR j := 1 to n DO BEGIN
      s := 0d0;
      FOR jj := 1 to n DO BEGIN
         s := s+v[j,jj]*tmp[jj];
      END;
      x[j] := s
   END
END;

end.
