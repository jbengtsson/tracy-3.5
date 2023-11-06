{ dynacc.inc }

function ellipse(xp, yp: graphvect; dpp: double; m: integer): double;
{ get beam ellipse parameters from tracking}                                   
CONST	np=5;
TYPE	glnpbynp =	ARRAY [1..np,1..np] OF double;
	glnarray =	ARRAY [1..np] OF double;
	glindx = 	ARRAY [1..np] OF integer;
VAR	alpha, beta, eta, etap, emitta:		double;
	b:					glnarray;
	a:					glnpbynp;
	index:					glindx;
	i, j, n:				integer;
	d, x, y, x2, y2, Gamma, RedEmi, CurlyH:	double;

PROCEDURE ludcmp(VAR a: glnpbynp; n: integer; VAR indx: glindx; VAR d: double);
CONST                                                                          
   tiny=1.0e-20;                                                               
VAR                                                                            
   k,j,imax,i: integer;                                                        
   sum,dum,big: double;                                                        
   vv: glnarray;                                                               
BEGIN                                                                          
   d := 1.0;                                                                   
   FOR i := 1 to n DO BEGIN                                                    
      big := 0.0;                                                              
      FOR j := 1 to n DO IF (abs(a[i,j]) > big) THEN big := abs(a[i,j]);       
      IF (big = 0.0) THEN BEGIN                                                
         writeln('pause in LUDCMP - singular matrix'); readln                  
      END;                                                                     
      vv[i] := 1.0/big                                                         
   END;                                                                        
   FOR j := 1 to n DO BEGIN                                                    
      IF (j > 1) THEN BEGIN                                                    
         FOR i := 1 to j-1 DO BEGIN                                            
            sum := a[i,j];                                                     
            IF (i > 1) THEN BEGIN                                              
               FOR k := 1 to i-1 DO BEGIN                                      
                  sum := sum-a[i,k]*a[k,j]                                     
               END;                                                            
               a[i,j] := sum                                                   
            END                                                                
         END                                                                   
      END;                                                                     
      big := 0.0;                                                              
      FOR i := j to n DO BEGIN                                                 
         sum := a[i,j];                                                        
         IF (j > 1) THEN BEGIN                                                 
            FOR k := 1 to j-1 DO BEGIN                                         
               sum := sum-a[i,k]*a[k,j]                                        
            END;                                                               
            a[i,j] := sum                                                      
         END;                                                                  
         dum := vv[i]*abs(sum);                                                
         IF (dum > big) THEN BEGIN                                             
            big := dum;                                                        
            imax := i                                                          
         END                                                                   
      END;                                                                     
      IF (j <> imax) THEN BEGIN                                                
         FOR k := 1 to n DO BEGIN                                              
            dum := a[imax,k];                                                  
            a[imax,k] := a[j,k];                                               
            a[j,k] := dum                                                      
         END;                                                                  
         d := -d;                                                              
         vv[imax] := vv[j]                                                     
      END;                                                                     
      indx[j] := imax;                                                         
      IF (j <> n) THEN BEGIN                                                   
         IF (a[j,j] = 0.0) THEN a[j,j] := tiny;                                
         dum := 1.0/a[j,j];                                                    
         FOR i := j+1 to n DO BEGIN                                            
            a[i,j] := a[i,j]*dum                                               
         END                                                                   
      END                                                                      
   END;                                                                        
   IF (a[n,n] = 0.0) THEN a[n,n] := tiny                                       
END;                                                                           
                                                                               
PROCEDURE lubksb(a: glnpbynp; n: integer; indx: glindx; VAR b: glnarray);      
VAR                                                                            
   j,ip,ii,i: integer;                                                         
   sum: double;                                                                
BEGIN                                                                          
   ii := 0;                                                                    
   FOR i := 1 to n DO BEGIN                                                    
      ip := indx[i];                                                           
      sum := b[ip];                                                            
      b[ip] := b[i];                                                           
      IF  (ii <> 0) THEN BEGIN                                                 
         FOR j := ii to i-1 DO BEGIN                                           
            sum := sum-a[i,j]*b[j]                                             
         END                                                                   
      END ELSE IF (sum <> 0.0) THEN BEGIN                                      
         ii := i                                                               
      END;                                                                     
      b[i] := sum                                                              
   END;                                                                        
   FOR i := n DOWNTO 1 DO BEGIN                                                
      sum := b[i];                                                             
      IF (i < n) THEN BEGIN                                                    
         FOR j := i+1 to n DO BEGIN                                            
            sum := sum-a[i,j]*b[j]                                             
         END                                                                   
      END;                                                                     
      b[i] := sum/a[i,i]                                                       
   END                                                                         
END;                                                                           
                                                                               
BEGIN                                                                          
  n:=np;
  for i:=1 to n do
  begin
    b[i]:=0;
    for j:=1 to n do a[i,j]:=0;
  end;
  for i:=1 to m do
  begin
    x:=xp[i]; x2:=Sqr(x); y:=yp[i]; y2:=Sqr(y);
    a[1,1]:=a[1,1]+Sqr(x2); a[1,2]:=a[1,2]+x*x2;   a[1,3]:=a[1,3]+x2*y2;
    a[1,4]:=a[1,4]+x2*y;    a[1,5]:=a[1,5]+x*x2*y; b[1]:=b[1]+x2;
    a[2,3]:=a[2,3]+x*y2;    a[2,4]:=a[2,4]+x*y;    b[2]:=b[2]+x;
    a[3,3]:=a[3,3]+Sqr(y2); a[3,4]:=a[3,4]+y*y2;   a[3,5]:=a[3,5]+x*y*y2;
    b[3]:=b[3]+y2;          b[4]:=b[4]+y;
  end;
  a[2,2]:=b[1]; a[2,5]:=a[1,4]; a[4,4]:=b[3]; a[4,5]:=a[2,3]; a[5,5]:=a[1,3];
  b[5]  :=a[2,4];
  for i:=1 to n do
    for j:=i+1 to n do
      a[j,i]:=a[i,j];
  LUDCMP(a,n,index,d); LUBKSB(a,n,index,b);
  RedEmi:=sign(b[1])/Sqrt(abs(b[1]*b[3]-0.25*Sqr(b[5])));
  Gamma :=b[1]*RedEmi; Beta :=b[3]*RedEmi; Alpha :=0.5*b[5]*RedEmi;
  eta :=0.5*RedEmi*(Alpha*b[4]-Beta *b[2]);
  etap :=0.5*RedEmi*(Alpha*b[2]-Gamma*b[4]);
  CurlyH:=Gamma*Sqr(eta)+2*Alpha*eta*etap+Beta*Sqr(etap);
  Emitta:=RedEmi+CurlyH;
{  writeln('  Parameters for ellipse-fit of boundary (rel. to closed orbit) :');
  writeln('  RedEmi = ',RedEmi:10,'  pi rad m');
  writeln('  Emitta = ',Emitta:10,'  pi rad m');
  writeln('  Beta   = ',Beta:10  ,'  m');
  writeln('  Alpha  = ',Alpha:10 );
  writeln('  Gamma  = ',Gamma:10 ,'  1/m');
}
  if dpp<>0 then
  begin
    eta:=eta/dpp;
    etap:=etap/dpp;
    CurlyH:=CurlyH/Sqr(dpp);
{    writeln('  eta    = ', eta:10  ,  ' m        ( = eta  + NLO*dp  )');
    writeln('  etap   = ', etap:10 ,  ' rad      ( = etap + NLO''*dp )');
    writeln('  CurlyH = ', CurlyH:10, ' pi rad m ( = H + NLO-Emi ) ');
}
  end
  else
  begin
{    writeln('  COD     = ', eta:10  ,  ' m       ( = NLO (nonlinear offset) )');
    writeln('  COD''   = ', etap:10 ,  ' rad');
    writeln('  COD-Emi = ', CurlyH:10, ' pi rad m');
}
  end;
  ellipse := emitta;
end;                                                                         

  procedure dynacc(var r : double; delta, toler : double; nturn : integer;
		   var Ax, Ay : double);
  var	lastn, lastpos, n, i:	integer;
	phi:			double;
	x, px, y, py, xa, ya:	graphvect;
 begin
  getcod(delta, lastpos); 
  if lastpos=globval.cell_nloc then begin
    getglobv_(globval);
    phi := pi-1d-6; 
    getdynap(r, phi, delta, toler, nturn, false);
    xa[1] := r*cos(phi); ya[1] := r*sin(phi);
    writeln;
    writeln('  Dynamical Aperture [mm]:         ',
	    1d3*xa[1]:10, ' ', 1d3*ya[1]:10);
    phi := pi/2.0-1d-6; getdynap(r, phi, delta, toler, nturn, false);
    xa[2] := r*cos(phi); ya[2] := r*sin(phi);
    writeln('  Dynamical Aperture [mm]:         ',
	    1d3*xa[2]:10, ' ', 1d3*ya[2]:10);
    if xa[1]*ya[2]<>0.0 then begin
      writeln('tracking with delta = ',delta:10:4);
      track(xa[1], 0.0, ya[1], 0.0, delta, nturn, lastn, lastpos, false);
      gettrack(n, x, px, y, py); 
      Ax := ellipse(x, px, delta, n);
      track(xa[2], 0.0, ya[2], 0.0, delta, nturn, lastn, lastpos, false);
      gettrack(n, x, px, y, py); 
      Ay := ellipse(y, py, delta, n);
    end
    else begin
      writeln('no tracking');
      ax:=0.0; ay:=0.0;
    end; 
  end
  else begin
    ax:=0.0; ay:=0.0;
  end;
  writeln('  Dynamical Acceptance [mmxmrad]:  ', 1d6*Ax:10, ' ', 1d6*Ay:10);
 end;
