  procedure tdiag(order, psuper : integer; Jx, Jy : double;
                  deltanu : double; grout : boolean);

{       Program to plot Tune-diagram and resonance widths.}

{ PROGRAM QRESMIN.FOR (QDIAG.FOR AT LBL, JOHAN BENGTSSON)
  PSI, S.ADAM, W.JOHO   14.8.1992
  PLOTS RESONANCE LINES  NX*QX+NY*QY=P
  SPECIAL VERSION WHICH SUPPRESSES RESONANCES FROM SKEW TERMS }

{ Transferred from Fortran to Pascal-s to become TRACY2 input file
  17.3.1995, A.Streun, PSI}

{ Assumes that "gettwiss" has been called. }

  %include 'fcoeff.inc'

type    double2 = array[1..2] of double;

const	{Labeling of resonance lines}		lblng= true;
        xorig=2.0; yorig=2.0; width=20.0; height=20.0; 
        nvrmax = 100;

var     NX,NY,nvr				: integer;
        P,KAPA,DE,Radius, cr, sr, Ex, Ey	: DOUBLE;
        nvrx, nvry, nvrp			: array[1..nvrmax] of integer;
        nux0, nuy0, fracnux0, fracnuy0, xe, ye	: double;
        sham, cham, degrad			: double;
        intnux0, intnuy0, i, ii			: integer;
 	cell					: celltype;
	str1, str2, str3, str			: tstring;
	ox, oy					: array[1..4] of double;
	xmin, xmax, ymin, ymax			: double;

procedure inc(var i:integer);
begin
  i:=succ(i);
end;


function Atan2(x,y: double): double;
var  a: double;
const eps=1E-10;
begin
  if Abs(x)<Eps then a:=Pi/2 
  else begin
    a:=ArcTan(Abs(y/x));
    if x<0 then a:=Pi-a;
  end;
  if y<0 then a:=-a;
  atan2:=a;
  end;


procedure RESWID (KAPA: double; NX,NY: integer; var DE: double);

var      NX1,NY1,N1              : integer;

begin
  N1:=ABS(NX)+ABS(NY)-2;
  NX1:=ABS(NX)-2;
  NY1:=ABS(NY)-2;

  IF (NY >=0) THEN
{   Sum resonance}
    DE:=2*KAPA*pwr((Radius/(2*PI)),(N1/2))*pwr(EX,(NX1/2))
         *pwr(EY,(NY1/2))*(sqr(NX)*EY+sqr(NY)*EX)
  ELSE
{   Difference resonance }
    DE:=2*KAPA*pwr((Radius/(2*PI)),(N1/2))*pwr(EX,(NX1/2))
         *pwr(EY,(NY1/2))*(abs(NX)*EY+abs(NY)*EX);
end;


procedure WRTINT (var STR: tstring; VAL: integer);

var       I          : integer;
          STR5       : tstring;

begin
  getint(str5,val,5);
  for I:=1 to 5 do
  begin
    IF (STR5.str[I] <> ' ') THEN
    begin
      inc(str.len);
      str.str[str.len]:=str5.str[i];
    end;
  end;
end;


procedure GETSTRing (NX,NY,P: integer; var STR: tstring);

begin
  str.len:=0;
  IF (NX <> 0) THEN
  begin
    IF (NX = -1) THEN
    begin;
      inc(str.len);
      str.str[str.len]:='-';
    end
    ELSE if not (nx=1) then
      WRTINT (STR,NX);
    concat (STR,'Qh');
  end;
  IF (NY <> 0) THEN
  begin
    IF ((NY = 1) AND (NX <> 0)) THEN
    begin
      inc(str.len);
      str.str[str.len]:='+';
    end
    ELSE begin
      IF (NY = -1) THEN
      begin
        inc(str.len);
        str.str[str.len]:='-';
      end
      ELSE if ((ny>1) and (nx<>0)) then
      begin
        inc(str.len);
        str.str[str.len]:='+';
      end;  
      WRTINT(STR,NY);
    end;
    Concat(str,'Qv');
  end;
  inc(str.len);
  str.str[str.len]:='=';
  WRTINT (STR,P);
end;

procedure WRTRES(var STR: tstring; NX,NY,P: integer);

var     ANG, X1, Y1, DX, DY, XE, YE, XH, PLF : double;
        xbuf, ybuf: array[1..2] of double;

begin
  XBUF[1]:=1.0*P*NX/((NX*NX + NY*NY)) + 200.0 * NY;
  XBUF[2]:=1.0*P*NX/((NX*NX + NY*NY)) - 200.0 * NY;
  YBUF[1]:=1.0*P*NY/((NX*NX + NY*NY)) - 200.0 * NX;
  YBUF[2]:=1.0*P*NY/((NX*NX + NY*NY)) + 200.0 * NX;
  IF ((NX <> 0) OR (NY <> 0)) THEN
  begin
    PLF := width/(XMAX-XMIN);
    IF ( height >  width ) then PLF := height/(YMAX-YMIN);
    DX:=XBUF[2]-XBUF[1];
    DY:=YBUF[2]-YBUF[1];
    ANG:=ATAN2(dx,dy);
    IF( (ANG  >  PI/2 + 0.0001)) THEN
    begin
      XH:=XBUF[1];
      XBUF[1]:=XBUF[2];
      XBUF[2]:=XH;
      XH:=YBUF[1];
      YBUF[1]:=YBUF[2];
      YBUF[2]:=XH;
      DX := -DX;
      DY := -DY;
      ANG := ANG - PI;
    end;
    IF(ABS(ANG)  <=  0.26*PI) THEN
    begin
      IF ((NX+NY)  <= 1) THEN
        X1:=2.3 + (abs(NX)+abs(NY)-1)*2.0
      ELSE
        X1:=22.0 - (abs(NX)+abs(NY)-1)*2.0;
      Y1 := 2.0 + PLF *(YBUF[1]+(DY*(XMIN-XBUF[1])/ DX) - YMIN) +
                (X1-2.0)*TAN_(ANG) + 0.15/COS(ANG);
      IF (Y1  <  2.05) THEN
        begin
        IF (ANG  >  0.0) THEN
          X1:=22.0 - 2.0
        ELSE
          X1:=2.3 + 2.0;
      Y1 := 2.0 + PLF *(YBUF[1]+(DY*(XMIN-XBUF[1])/ DX)-YMIN)+
                (X1-2.0)*TAN_(ANG) + 0.15/COS(ANG);
      IF(Y1  <  2.05) OR ( Y1  >  20.8) THEN
      begin
        Y1:=-20; X1:=-20;
      end;
    end;
    IF (Y1  >  20.8) THEN
    begin
      IF (ANG  >  0.0) THEN
         X1:=2.3 + 2.0
       ELSE
         X1:=22.0 - 2.0;
     end;
     X1:=2.3 + 2.0;
     Y1 := 2.0 + PLF *(YBUF[1]+(DY*(XMIN-XBUF[1])/ DX)-YMIN)+
              (X1-2.0)*TAN_(ANG) + 0.15/COS(ANG);
     IF (Y1  <  2.05) or ( Y1  >  20.8) THEN
     begin
       Y1:=-20; X1:=-20;
     end;
   end
   ELSE begin
{C  HIGH ANGLE CASES}
     IF ((NX+NY)  <= 1) THEN
      Y1:=2.3 + (abs(NX)+abs(NY)-1)*2.0 + 0.2*SIN(ANG)
    ELSE
      Y1:=22.0 - (abs(NX)+abs(NY)-1)*2.0+ 0.2*SIN(ANG);
      X1 := 2.0 + PLF *(XBUF[1]+(DX*(YMIN-YBUF[1])/ DY) - XMIN) +
              (Y1 - 2.0)*DX/DY - 0.15/SIN(ANG);
      IF (X1  <  2.34) THEN
      begin
        IF (ANG  >  0.0) THEN
          Y1:=22.0 - 2.0+ 0.2*SIN(ANG)
        ELSE
          Y1:=2.3 + 2.0 + 0.2*SIN(ANG);
        X1 := 2.0 + PLF *(XBUF[1]+(DX*(YMIN-YBUF[1])/ DY)
                  - XMIN) + (Y1 - 2.0)*DX/DY - 0.15/SIN(ANG);
        IF(X1  <  2.34 ) OR ( X1  >  20.8) THEN
        begin
          Y1:=-20; X1:=-20;
        END;
      END;
    IF (X1  >  20.8) THEN
    begin
      IF (ANG  >  0.0) THEN
        Y1:=2.3 + 2.0 + 0.2*SIN(ANG)
      ELSE
        Y1:=22.0 - 2.0+ 0.2*SIN(ANG);
      X1 := 2.0 + PLF *(XBUF[1]+(DX*(YMIN-YBUF[1])/ DY)
               - XMIN) + (Y1 - 2.0)*DX/DY - 0.15/SIN(ANG);
      IF(X1  <  2.34 ) OR ( X1  >  20.8) THEN
      begin
        Y1:=-20; X1:=-20;
      end;
    end;
  end;
  IF ((DX <> 0.0) or (DY <> 0.0)) THEN
    GRATEXs (X1,Y1, 0.25, STR.str, str.len, ATAN2(dx,dy)/DEGRAD,XE,YE,1);
  end;
end;



procedure caption;
var     STR, STR1, STR2, STR3           : tstring;
	XE, YE				: double;
	I				: integer;
begin
  gracoi(1);
  str.len:=0; str1.len:=0; str2.len:=0; str3.len:=0;
  getstr(STR,'Q-DIAGRAM FOR RESONANCES UP TO ');
  getstr(STR2,'SUPERPERIOD =');
  getint(STR1,order,3);
  str3.len:=0;
  getint(STR3,PSUPER,2);
  for i:=1 to str3.len do
  begin
    inc(str2.len);
    STR2.str[str2.len]:=STR3.str[i];
  end;
  for I:=1 to 3 do
  begin
    IF (STR1.str[I] <> ' ') THEN
    begin
      inc(str.len);
      STR.str[str.len]:=STR1.str[i];
    end;
  end;
  IF (ORDER = 2) THEN
    concat(STR,'nd')
  ELSE begin
    IF (ORDER = 3) THEN
      concat(STR,'rd')
    ELSE
      concat(STR,'th');
  end;
  concat(str,' ORDER');
  concat(STR,' (no skew term)');

  str1.len:=0; str3.len:=0;
  for i:=15 to 31 do begin
    inc(str1.len);
    str1.str[str1.len]:=str.str[i];
    end;
  for i:=32 to str.len do begin
    inc(str3.len);
    str3.str[str3.len]:=str.str[i];
    end;
  str.len:=14;
  GRATEXs (23.5, 21.0, 0.5, STR.str, str.len, 0.0, XE,YE,1);
  GRATEXs (23.5, 19.5, 0.5, STR1.str,str1.len, 0.0, XE,YE,1);
  GRATEXs (23.5, 18.0, 0.5, STR3.str,str3.len, 0.0, XE,YE,1);
  GRATEXs (23.5, 16.5, 0.5, STR2.str, str2.len, 0.0, XE,YE,1);
  getstr(STR1,'Resonance widths');
  getstr(str2,'calculated for');
  GRATEXs (23.5, 9.0, 0.5, STR1.str, str1.len, 0.0, XE,YE,1);
  GRATEXs (23.5, 7.8, 0.5, STR2.str, str2.len, 0.0, XE,YE,1);
  str2.len:=0;
  getreale(str2,ex,11,3);
  str.len:=0;
  getstr(str,'2J#il#x#cl# = ');
  for i:=1 to str2.len do
  begin
    inc(str.len);
    str.str[str.len]:=str2.str[i];
  end;
  str3.len:=0;
  getreale(str3,ey,11,3);
  str1.len:=0;
  getstr(str1,'2J#il#y#cl# = ');
  for i:=1 to str3.len do
  begin
    inc(str1.len);
    str1.str[str1.len]:=str3.str[i];    
  end;
  GRATEXs(24.0, 6.6, 0.5, STR.str, str.len , 0.0, XE,YE,1);
  GRATEXs(24.0, 5.4, 0.5, STR1.str,str1.len, 0.0, XE,YE,1);
end;

procedure findres(ORDER, PSUPER: integer);

var     NX, NY, I, START, P, P1, P2	: integer;
        PS				: double;

begin

  nvr:=0;
  for NX:=0 to ORDER do
  begin
    IF (NX = 0) THEN
      START:=1
    ELSE
      START:=-ORDER;
    for NY:=START to ORDER do
    begin
      I:=abs(NX)+abs(NY);
{             ELIMINATION OF SKEW TERMS}
      if not ( (ABS(NY) = 3) or ( (ABS(NY) = 1) and (ORDER >=2) ) ) then
      begin
        IF ((I <= ORDER) AND (I > 0)) THEN
        begin
          IF (NY >=0) THEN
          begin
{                limiting p's for sum resonance} {and one more...}
            P1:=trunc(NX*XMIN+NY*YMIN)  -1;
            P2:=trunc(NX*XMAX+NY*YMAX)+1+1;
          end
          ELSE begin
{               limiting p's for difference resonance}
            P1:=trunc(NX*XMIN+NY*YMAX)  -1;
            P2:=trunc(NX*XMAX+NY*YMIN)+1+1;
          end;
          for P:=P1 to P2 do
          begin
            PS:=P/PSUPER;
            IF( (PS = trunc(PS)) or ((I <= 2) and((NX = 0) OR (NY = 0))) ) THEN
            begin
	      inc(nvr);
	      nvrx[nvr]:=nx; nvry[nvr]:=ny; nvrp[nvr]:=p;
            end;
          end; {for P}
        END;
      end; {if not}
    end; {for NY}
  end; {for NX}
END;


procedure baktra (var x, y :  double);
{physical to geometrical coordinates}
begin
  x:=xorig+width *(x-xmin)/(xmax-xmin);
  y:=yorig+height*(y-ymin)/(ymax-ymin);
end;

procedure MakeLine (a, b: integer; pd: double; var x1, y1, x2, y2: double);
var  xu, yu, xo, yo: double;
begin
  xu:=xmin; xo:=xmax; yu:=ymin; yo:=ymax;
  if (a=0) or (b=0) then 
  begin
    if a=0 then 
    begin
      x1 := xu;   x2 := xo;
      y1 := pd/b; y2:=y1;
    end
    else 
    begin {b=0}
      y1 := yu;   y2 := yo;
      x1 := pd/a; x2:=x1;
    end;
  end
  else  {a,b <> 0}
  begin
    y1 := (pd-a*xu)/b;
    y2 := (pd-a*xo)/b;
    x1 := xu; x2 := xo;
  end;
end;

procedure PloBand (a, b, p: integer;  d: double);
var i: integer;
    x1, y1, x2, y2: double;
begin
  if grout then
  begin
    gracoi(13); 
    grapat(2);
    makeline(a,b,p+d/2,ox[1],oy[1],ox[2],oy[2]);
    makeline(a,b,p-d/2,ox[4],oy[4],ox[3],oy[3]);
    for i:=1 to 4 do baktra(ox[i],oy[i]); 
    grafpl(ox,oy,4);
  end
  else 
    writeln('nx, ny, p, de : ',a:3,b:4,p:4,'  ',d:10);
end;


BEGIN
  degrad := pi/180.0;

  Ex := 2.0*Jx; Ey := 2.0*Jy;

  getelem(globval.cell_nloc,cell);
  nux0:=cell.nu[1]; intnux0:=Trunc(nux0); fracnux0:=nux0-intnux0;
  nuy0:=cell.nu[2]; intnuy0:=Trunc(nuy0); fracnuy0:=nuy0-intnuy0;
  Radius:=cell.s/2/pi;
  xmin:=Trunc(10*(nux0-deltanu/2)  )*0.1;
  xmax:=Trunc(10*(nux0+deltanu/2)+1)*0.1;
  ymin:=Trunc(10*(nuy0-deltanu/2)  )*0.1;
  ymax:=Trunc(10*(nuy0+deltanu/2)+1)*0.1;

  if grout then 
  begin
    GRAOPE;
    GRAVWP(xorig, yorig, width, height);
    GRAWND(XMIN,XMAX,YMIN,Ymax);
    GRAXI1(-4,999999,999999,1E36);
    GRAYI1(-4,999999,999999,1E36);
  end;

  findres (ORDER,PSUPER);
  if grout then 
  begin
    gravcf(2);
    gravcs(xorig,yorig,xorig+width,yorig+height);
  end
  else  
    gracls;
  for i:=1 to nvr do 
  begin
    if (nvrx[i]=1) and (nvry[i]= 0) and not grout then
    begin
      fcoeff(2, 1, 0, 0, nvrp[i], 1/3, cr, sr);
      kapa:=sqrt(sqr(cr)+sqr(sr));
        writeln('fcoeff: c, s, kapa: ',cr:12,' ',sr:12,' ',kapa:12);
      fcoeff(1, 0, 1, 1, nvrp[i], -1, cr, sr) ;
      kapa:=sqrt(sqr(cr)+sqr(sr));
        writeln('fcoeff: c, s, kapa: ',cr:12,' ',sr:12,' ',kapa:12);
    end                                           
    else if (nvrx[i]=3) and (nvry[i]= 0) then
    begin
      fcoeff(3, 0, 0, 0, nvrp[i], 1/3, cr, sr);
      kapa:=sqrt(sqr(cr)+sqr(sr));
      if not grout then 
        writeln('fcoeff: c, s, kapa: ',cr:12,' ',sr:12,' ',kapa:12);
      RESWID (KAPA, NvrX[i],NvrY[i],DE);
      PloBand(NvrX[i],Nvry[i],NVRP[i],DE);
    end                                           
    else if (nvrx[i]=1) and (nvry[i]= 2) then
    begin
      fcoeff(1, 0, 0, 2, nvrp[i], -1, cr, sr);
      kapa:=sqrt(sqr(cr)+sqr(sr));
      if not grout then 
        writeln('fcoeff: c, s, kapa: ',cr:12,' ',sr:12,' ',kapa:12);
      RESWID (KAPA, NvrX[i],NvrY[i],DE);
      PloBand(NvrX[i],Nvry[i],NVRP[i],DE);
    end                                           
    else if (nvrx[i]=1) and (nvry[i]=-2) then
    begin
      fcoeff(1, 0, 2, 0, nvrp[i], -1, cr, sr);
      kapa:=sqrt(sqr(cr)+sqr(sr));
      if not grout then 
        writeln('fcoeff: c, s, kapa: ',cr:12,' ',sr:12,' ',kapa:12);
      RESWID (KAPA, NvrX[i],NvrY[i],DE);
      PloBand(NvrX[i],Nvry[i],NVRP[i],DE);
    end ;
  end;
  if grout then 
  begin
    for i:=1 to nvr do 
    begin
      gracoi(1);
      makeline(nvrx[i], nvry[i], nvrp[i],ox[1],oy[1],ox[2],oy[2]);
      graph1(ox,oy,2,1,0,0.0);
      GETSTRing (NvrX[i],Nvry[i],nvrP[i],STR);
      IF LBLNG then WRTRES(STR,NvrX[i],NvrY[i],nvrP[i]);
    end;
  end;                                      
  if grout then begin
    caption;
    ox[1]:=nux0; oy[1]:=nuy0; ox[2]:=ox[1]; oy[2]:=oy[1];
    gracoi(1);
    graph1(ox,oy,2,0,2,0.3);
    GRACLS;
  end;
END;
