
%include 'svd.def'

[global] PROCEDURE svdcmp(VAR a: svdmat; m, n: integer;
			  VAR w: svdarray; VAR v: svdmat);

LABEL 1,2,3;

VAR
   nm,l,k,j,its,i: integer;
   z,y,x,scale,s,h,g,f,c,anorm: double;
   rv1: ARRAY [1..mnp] OF double;


FUNCTION sign(a,b: double): double;

   BEGIN
      IF (b >= 0d0) THEN sign := abs(a) ELSE sign := -abs(a)
   END;

FUNCTION max(a,b: double): double;

   BEGIN
      IF (a > b) THEN max := a ELSE max := b
   END;


BEGIN
   g := 0d0;
   scale := 0d0;
   anorm := 0d0;
   FOR i := 1 to n DO BEGIN
      l := i+1;
      rv1[i] := scale*g;
      g := 0d0;
      s := 0d0;
      scale := 0d0;
      IF (i <= m) THEN BEGIN
         FOR k := i to m DO BEGIN
            scale := scale+abs(a[k,i])
         END;
         IF (scale <> 0d0) THEN BEGIN
            FOR k := i to m DO BEGIN
               a[k,i] := a[k,i]/scale;
               s := s+a[k,i]*a[k,i]
            END;
            f := a[i,i];
            g := -sign(sqrt(s),f);
            h := f*g-s;
            a[i,i] := f-g;
            IF (i <> n) THEN BEGIN
               FOR j := l to n DO BEGIN
                  s := 0d0;
                  FOR k := i to m DO BEGIN
                     s := s+a[k,i]*a[k,j]
                  END;
                  f := s/h;
                  FOR k := i to m DO BEGIN
                     a[k,j] := a[k,j]+
                        f*a[k,i]
                  END
               END
            END;
            FOR k := i to m DO BEGIN
               a[k,i] := scale*a[k,i]
            END
         END
      END;
      w[i] := scale*g;
      g := 0d0;
      s := 0d0;
      scale := 0d0;
      IF ((i <= m) AND (i <> n)) THEN BEGIN
         FOR k := l to n DO BEGIN
            scale := scale+abs(a[i,k])
         END;
         IF (scale <> 0d0) THEN BEGIN
            FOR k := l to n DO BEGIN
               a[i,k] := a[i,k]/scale;
               s := s+a[i,k]*a[i,k]
            END;
            f := a[i,l];
            g := -sign(sqrt(s),f);
            h := f*g-s;
            a[i,l] := f-g;
            FOR k := l to n DO BEGIN
               rv1[k] := a[i,k]/h
            END;
            IF (i <> m) THEN BEGIN
               FOR j := l to m DO BEGIN
                  s := 0d0;
                  FOR k := l to n DO BEGIN
                     s := s+a[j,k]*a[i,k]
                  END;
                  FOR k := l to n DO BEGIN
                     a[j,k] := a[j,k]
                        +s*rv1[k]
                  END
               END
            END;
            FOR k := l to n DO BEGIN
               a[i,k] := scale*a[i,k]
            END
         END
      END;
      anorm := max(anorm,(abs(w[i])+abs(rv1[i])))
   END;
   FOR i := n DOWNTO 1 DO BEGIN
      IF (i < n) THEN BEGIN
         IF (g <> 0d0) THEN BEGIN
            FOR j := l to n DO BEGIN
               v[j,i] := (a[i,j]/a[i,l])/g
            END;
            FOR j := l to n DO BEGIN
               s := 0d0;
               FOR k := l to n DO BEGIN
                  s := s+a[i,k]*v[k,j]
               END;
               FOR k := l to n DO BEGIN
                  v[k,j] := v[k,j]+s*v[k,i]
               END
            END
         END;
         FOR j := l to n DO BEGIN
            v[i,j] := 0d0;
            v[j,i] := 0d0
         END
      END;
      v[i,i] := 1d0;
      g := rv1[i];
      l := i
   END;
   FOR i := n DOWNTO 1 DO BEGIN
      l := i+1;
      g := w[i];
      IF (i < n) THEN BEGIN
         FOR j := l to n DO BEGIN
            a[i,j] := 0d0
         END
      END;
      IF (g <> 0d0) THEN BEGIN
         g := 1d0/g;
         IF (i <> n) THEN BEGIN
            FOR j := l to n DO BEGIN
               s := 0d0;
               FOR k := l to m DO BEGIN
                  s := s+a[k,i]*a[k,j]
               END;
               f := (s/a[i,i])*g;
               FOR k := i to m DO BEGIN
                  a[k,j] := a[k,j]+f*a[k,i]
               END
            END
         END;
         FOR j := i to m DO BEGIN
            a[j,i] := a[j,i]*g
         END
      END ELSE BEGIN
         FOR j := i to m DO BEGIN
            a[j,i] := 0d0
         END
      END;
      a[i,i] := a[i,i]+1d0
   END;
   FOR k := n DOWNTO 1 DO BEGIN
      FOR its := 1 to 30 DO BEGIN
         FOR l := k DOWNTO 1 DO BEGIN
            nm := l-1;
            IF ((abs(rv1[l])+anorm) = anorm) THEN GOTO 2;
            IF ((abs(w[nm])+anorm) = anorm) THEN GOTO 1
         END;
1:       c := 0d0;
         s := 1d0;
         FOR i := l to k DO BEGIN
            f := s*rv1[i];
            IF ((abs(f)+anorm) <> anorm) THEN BEGIN
               g := w[i];
               h := sqrt(f*f+g*g);
               w[i] := h;
               h := 1d0/h;
               c := (g*h);
               s := -(f*h);
               FOR j := 1 to m DO BEGIN
                  y := a[j,nm];
                  z := a[j,i];
                  a[j,nm] := (y*c)+(z*s);
                  a[j,i] := -(y*s)+(z*c)
               END
            END
         END;
2:       z := w[k];
         IF (l = k) THEN BEGIN
            IF (z < 0d0) THEN BEGIN
               w[k] := -z;
               FOR j := 1 to n DO BEGIN
                  v[j,k] := -v[j,k]
               END
            END;
            GOTO 3
         END;
         IF (its = 30) THEN BEGIN
            writeln ('no convergence in 30 SVDCMP iterations'); readln
         END;
         x := w[l];
         nm := k-1;
         y := w[nm];
         g := rv1[nm];
         h := rv1[k];
         f := ((y-z)*(y+z)+(g-h)*(g+h))/(2d0*h*y);
         g := sqrt(f*f+1d0);
         f := ((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x;
         c := 1d0;
         s := 1d0;
         FOR j := l to nm DO BEGIN
            i := j+1;
            g := rv1[i];
            y := w[i];
            h := s*g;
            g := c*g;
            z := sqrt(f*f+h*h);
            rv1[j] := z;
            c := f/z;
            s := h/z;
            f := (x*c)+(g*s);
            g := -(x*s)+(g*c);
            h := y*s;
            y := y*c;
            FOR nm := 1 to n DO BEGIN
               x := v[nm,j];
               z := v[nm,i];
               v[nm,j] := (x*c)+(z*s);
               v[nm,i] := -(x*s)+(z*c)
            END;
            z := sqrt(f*f+h*h);
            w[j] := z;
            IF (z <> 0d0) THEN BEGIN
               z := 1d0/z;
               c := f*z;
               s := h*z
            END;
            f := (c*g)+(s*y);
            x := -(s*g)+(c*y);
            FOR nm := 1 to m DO BEGIN
               y := a[nm,j];
               z := a[nm,i];
               a[nm,j] := (y*c)+(z*s);
               a[nm,i] := -(y*s)+(z*c)
            END
         END;
         rv1[l] := 0d0;
         rv1[k] := f;
         w[k] := x
      END;
3: END
END;
