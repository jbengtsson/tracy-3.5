var	cinf, coutf	: text;

procedure cogir(seed : integer; gdxrms, gdzrms : double);
{  Closed Orbit with correlated magnet misalignments, i.e. assuming,
   that elements are mounted on girders. 
     
   Girders are assumed to be connected by joints, if they are close 
   to each other, or if the "fargird" option is activated. 
   
   Three kinds of displacement errors: 
   1. Girder joints (and loose ends resp. if no connection found)
   2. Tolerance of the joint itself (joint play)
   3. Relative error on elements on girders.

   BPMs may be included in girders or read out absolute. The name of
   the BPMs in the lattice file has to be "MON" !
     
   Girders end and start at a marker whose name begins with "G".
   Do not use with split elements, otherwise the COD is calculated wrong.      }

const  ngrmax=300; igrmax=600; 
       nbfammax=10;
       maxloomax=1000;
       dgeps=1E-10; {max. distance [i].s1-[i-1].s2 to connect girders}

type  girdertype = record 
        ngr : integer;
        f   : array[1..ngrmax] of integer;
        k   : array[1..ngrmax] of integer;
        s1, s2, x1, x2, z1, z2: double;
      end;    
        
var   cell	: celltype;
      k, k1, k2, maxloo, lastpos, icha, ip, iloo, ib, nbfam, ipos, 
        nbpm,j, igr, i: integer;
      gir  : array[1..igrmax] of girdertype;
      memxc, memyc: array[1..maxloomax] of double;
      gsetout, esetout, found, includeMON, fargird: boolean;
      a: char;
      mx, my, ampx, ampy, meanx, meany, xmax, ymax, xclob, yclob, 
        jdxrms,jdzrms,edxrms,edzrms :double;
      str: tstring;

procedure GirAssemble (includemon: boolean);
var itismon: boolean; i, numel: integer;
begin
  igr:=0;
  for i:=1 to globval.cell_nloc do begin
    getelem(i,cell);
    if (cell.elem.pkind=marker) and (cell.elem.pname[1]='g') then begin
      if igr=0 then igr:=igr+1
      else begin 
        gir[igr].s2:=cell.s; {end old gir}
        if gir[igr].ngr>0 then igr:=igr+1;
      end;
      numel:=0;
      gir[igr].s1:=cell.s;
    end
    else if (cell.elem.pkind=Mpole) then begin
      itismon := (cell.elem.pname[1]='m') and (cell.elem.pname[2]='o') and
              (cell.elem.pname[3]='n');
      if (itismon and includemon) or not itismon then begin
        numel:=numel+1;
        gir[igr].ngr:=numel;
        gir[igr].f[numel]:=cell.fnum;
        gir[igr].k[numel]:=cell.knum;
      end;
    end; 
  end; 
  if gir[igr].ngr=0 then igr:=igr-1;
end;

procedure GirShowComp;
var i, j, ieol: integer;
begin
  writeln('Number of girders = ',igr:3);
  for i:=1 to igr do begin
    writeln;
    writeln('Girder ',i:3,' from ',gir[i].s1:8:3, ' to ', gir[i].s2:8:3,
                ' contains ',gir[i].ngr:3,' elements:');
    ieol:=0;
    for j:=1 to gir[i].ngr do begin
      ipos:=elem_getpos(gir[i].f[j],gir[i].k[j]);
      getelem(ipos,cell);
      write(cell.elem.pname[1],cell.elem.pname[2],cell.elem.pname[3],
          cell.elem.pname[4],cell.knum:3);
      ieol:=ieol+1; 
      if (ieol=8)or(j=gir[i].ngr) then begin
        ieol:=0;
        writeln;
      end
      else write(' | '); 
    end;
  end;
end;

procedure GirSet (gsetout, esetout, fargird: boolean);
var sj: double;
begin
  if gsetout then writeln('Girder displacement error set up:');
  gir[1].x1:=gdxrms*normranf;  
  gir[1].x2:=gdxrms*normranf;
  gir[1].z1:=gdzrms*normranf;  
  gir[1].z2:=gdzrms*normranf;
  for i:=2 to igr do begin
    {assume a joint if girders are close, otherwise loose end:}
    if (abs(gir[i].s1-gir[i-1].s2)<dgeps) or fargird then begin 
      gir[i].x1:=gir[i-1].x2+jdxrms*normranf;
      gir[i].z1:=gir[i-1].z2+jdzrms*normranf;
    end
    else begin {start new}
      gir[i].x1:=gdxrms*normranf;
      gir[i].z1:=gdzrms*normranf;
    end;
    gir[i].x2:=gdxrms*normranf;
    gir[i].z2:=gdzrms*normranf;
  end;
  {close the circle:}
  if (abs(gir[1].s1-gir[igr].s2)<dgeps) or fargird then begin
    gir[1].x1:=gir[igr].x2+jdxrms*normranf;
    gir[1].z1:=gir[igr].z2+jdzrms*normranf;
  end;
 
  for i:=1 to igr do begin
    if gsetout then begin
      writeln('Girder ',i:3,': s1=',gir[i].s1:8:3,' m. ', 
              'dx1=',gir[i].x1*1E6:6:1,',   dz1=',gir[i].z1*1E6:6:1,' um');
      writeln(    '          : s2=',gir[i].s2:8:3,' m. ',
              'dx2=',gir[i].x2*1E6:6:1,',   dz2=',gir[i].z2*1E6:6:1,' um');
    end;

    {Elements on girders:}
    for j:=1 to gir[i].ngr do begin
      ipos:=elem_getpos(gir[i].f[j],gir[i].k[j]);
      getelem(ipos,cell);
      sj:=cell.s-cell.elem.pl/2; {centre of element}
      cell.dS[1]:=(gir[i].x1*(gir[i].s2-sj)+gir[i].x2*(sj-gir[i].s1))/
                  (gir[i].s2-gir[i].s1) + edxrms*normranf;
      cell.dS[2]:=(gir[i].z1*(gir[i].s2-sj)+gir[i].z2*(sj-gir[i].s1))/
                  (gir[i].s2-gir[i].s1) + edzrms*normranf;
      if esetout then
        writeln('   ',cell.elem.pname[1],cell.elem.pname[2],cell.elem.pname[3],
          cell.elem.pname[4],'/',gir[i].k[j]:2,': sj=',sj:8:3,' m. ',
          ' dx=',cell.dS[1]*1D6:6:1,' um, dz=',cell.dS[2]*1D6:6:1,' um');
      putelem(ipos,cell);
    end;
  end;
end;

procedure ascodstat (var xmean, ymean, xsigma, ysigma, xmax, ymax: double);
{supplements the no longer compatible procedure codstat}
var x, y, sx, sy, sx2, sy2: double;
    k: integer;
begin
  getcod(0.0,lastpos);
  getglobv_(globval);
  cell_pass(0, globval.cell_nloc, globval.codvect, lastpos);
  nbpm := getnkid(globval.bpm);
  sx:=0; sy:=0; sx2:=0; sy2:=0; xmax:=0; ymax:=0;
  for k:=1 to nbpm do begin
    getElem(elem_getpos(globval.bpm, k), cell);
    x:=cell.beampos[1]; y:=cell.beampos[3];
    sx :=sx +    x ;    sy :=sy +    y;
    sx2:=sx2+sqr(x);    sy2:=sy2+sqr(y);
    if abs(x)>xmax then xmax:=abs(x); if abs(y)>ymax then ymax:=abs(y);
  end;
  xmean :=sx/nbpm; ymean:=sy/nbpm;
  xsigma:=sqrt(sx2/nbpm-sqr(xmean)); ysigma:=sqrt(sy2/nbpm-sqr(ymean));
end;

{.........................................................................}

begin
  getstr(str,'cogir.dat'); reset_(cinf,str.str);

  writeln;
  writeln('globval.cell_nloc    = ',globval.cell_nloc); 
  writeln('globval.elem_nfam = ',globval.elem_nfam);

  writeln('Apply errors also to BPMs (name=MON)? (j/n) ');
  readln(cinf, a);
  includeMON:= not ((a='n')or(a='N'));
  if includeMON then writeln('BPMs with errors')
  else writeln('BPMs without errors');

  GirAssemble (includemon);
 
  writeln('Show girder composition (Y/n) ?');
  readln(cinf, a);
  if not ( (a='n')or(a='N') ) then GirShowComp;
  writeln;
  writeln('Also connect girders OVER empty girders (y/N) ?');
  readln(cinf, a);
  fargird:= ( (a='y')or(a='Y') );
  writeln;
  writeln('Input of error amplitudes for gaussian errors, cut 2 sigma');
  writeln('----------------------------------------------------------');
  writeln('Give rms errors for displacements in micron,',
               ' horizontal and vertical:');
{  writeln('__ Absolute displacement of girder joints and ends :');
  readln(cinf, gdxrms,gdzrms);}
  writeln('__ Relative displacement WITHIN girder joints (joint play):');
  readln(cinf, jdxrms,jdzrms);
  writeln('__ Relative displacement of elements ON a girder:');
  readln(cinf, edxrms,edzrms);
  writeln;
  writeln('Number of seeds for calculation of mean rms COD :');
  writeln;
  readln(cinf, maxloo);
  if maxloo>maxloomax then maxloo:=maxloomax;

  writeln('rms Girder error  : dx=',1e6*gdxrms:5:0,' um, dz=',
	  1e6*gdzrms:5:0,' um');
  writeln('rms Joint  error  : dx=',1e6*jdxrms:5:0,' um, dz=',
	  1e6*jdzrms:5:0,' um');
  writeln('rms Element error : dx=',1e6*edxrms:5:0,' um, dz=',
	  1e6*edzrms:5:0,' um');

{  gdxrms:=gdxrms*1E-6;  gdzrms:=gdzrms*1E-6; 
  jdxrms:=jdxrms*1E-6;  jdzrms:=jdzrms*1E-6; 
  edxrms:=edxrms*1E-6;  edzrms:=edzrms*1E-6; }

  writeln('Extended output on girder/element error setting [0 0]...[1 1] ?');
  readln(cinf, k1,k2);
  gsetout:= k1=1;  esetout:= k2=1;

{  writeln('Creating plotfile cogir.plt');
  getstr(str,'cogir.plt'); rewrite_(coutf,str.str);
  writeln(coutf,gdxrms:20,' ',gdzrms:20,' 0.0');
  writeln(coutf,jdxrms:20,' ',jdzrms:20,' 0.0');
  writeln(coutf,edxrms:20,' ',edzrms:20,' 0.0');
  writeln(coutf,maxloo:20);}
  
  iniranf(seed); setrancut(2.0);
  
  mx:=0; my:=0;
  for iloo:=1 to maxloo do begin
    GirSet (gsetout, esetout, fargird);
    ascodstat(meanx,meany,xclob,yclob,xmax,ymax);
    writeln(iloo:5,' CODmean  : x=',meanx*1000:8:3,' mm; y=',meany*1000:8:3,' mm');
    writeln('      CODsigma : x=',xclob*1000:8:3,' mm; y=',yclob*1000:8:3,' mm');
    writeln('      CODmax   : x=',xmax*1000:8:3,' mm; y=',ymax*1000:8:3,' mm');
    writeln;
    mx:=mx + xclob;
    my:=my + yclob;
    memxc[iloo]:=xclob*1D3;
    memyc[iloo]:=yclob*1D3;
{    writeln(coutf,memxc[iloo]:20:8,memyc[iloo]:20:8);}
  end;

  mx:=mx/maxloo; my:=my/maxloo;
  writeln('X mean sigma = ',mx*1000:10:4,' mm');
  writeln('Z mean sigma = ',my*1000:10:4,' mm');
  writeln;
  close(cinf);{ close(coutf);}
end;
