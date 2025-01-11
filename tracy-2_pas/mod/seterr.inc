{   as-080198

follows bcosys: now systematic multipole errors to be set for 
different magnet types independently

   Closed Orbit with correlated magnet misalignments, i.e. assuming,
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
   Do not use with split elements, otherwise the COD is calculated wrong.

   Girder marker "K2" defines a dependent girder, sitting ON the two
   adjacent girder and treated like an element.

}




const
   ngrmax=300; 
   igrmax=1000; 
   maxloomax=1000;

type  
   girdertype = record 
     ngr, bset, kind : integer;
     f   : array[1..ngrmax] of integer;
     k   : array[1..ngrmax] of integer;
     s1, s2, x1, x2, z1, z2: double;
   end;    

var 
   cinf, outf, gout: text;
   
   gir  : array[1..igrmax] of girdertype;
   igr: integer;
   gsetout, esetout, found, includeMON, fargird: boolean;
   cell	: celltype;
   maxloo: integer;
   dbbrms, gdxrms,gdzrms,jdxrms,jdzrms,edxrms,edzrms,
   mp11, mp12, mp13, mp24, mp26, mp36, mp39, circumference: double;
   str: tstring;

{.....................................................................}

procedure GirSet (gfrac, jfrac, efrac: double);
                 
const
    dgeps=1.0E-10; {max. distance [i].s1-[i-1].s2 to connect girders}
var 
    dbb, sj: double;
    i, j, io, iv, ipos: integer;
    k2flag: boolean;
begin
  if gsetout then writeln('Girder displacement error set up:');
  gir[1].x1:=gdxrms*normranf*gfrac;  
  gir[1].x2:=gdxrms*normranf*gfrac;
  gir[1].z1:=gdzrms*normranf*gfrac;  
  gir[1].z2:=gdzrms*normranf*gfrac;
  iv:=1;
  k2flag:=false;
  for i:=2 to igr do begin
    if gir[i].kind<2 then begin
      {assume a joint if girders are close, otherwise loose end:}
      if (abs(gir[i].s1-gir[iv].s2)<dgeps) or fargird then begin 
        gir[i].x1:=gir[iv].x2+jdxrms*normranf*jfrac;
        gir[i].z1:=gir[iv].z2+jdzrms*normranf*jfrac;
      end
      else begin {start new}
        gir[i].x1:=gdxrms*normranf*gfrac;
        gir[i].z1:=gdzrms*normranf*gfrac;
      end;
      gir[i].x2:=gdxrms*normranf*gfrac;
      gir[i].z2:=gdzrms*normranf*gfrac;
      iv:=i;
    end
    else k2flag:=true;
  end;
  {close the circle:}
  if (abs(gir[1].s1+circumference-gir[igr].s2)<dgeps) or fargird then begin
    gir[1].x1:=gir[igr].x2+jdxrms*normranf*jfrac;
    gir[1].z1:=gir[igr].z2+jdzrms*normranf*jfrac;
  end;

{2nd loop: set errors from adjacent girders}
  for i:=2 to igr-1 do begin
    if gir[i].kind=2 then begin
      sj:=gir[i].s1;
      gir[i].x1:= (gir[i-1].x1*(gir[i-1].s2-sj)+gir[i-1].x2*(sj-gir[i-1].s1))/
                  (gir[i-1].s2-gir[i-1].s1) + edxrms*normranf*efrac;
      gir[i].z1:= (gir[i-1].z1*(gir[i-1].s2-sj)+gir[i-1].z2*(sj-gir[i-1].s1))/
                  (gir[i-1].s2-gir[i-1].s1) + edzrms*normranf*efrac;
      sj:=gir[i].s2;
      gir[i].x2:= (gir[i+1].x1*(gir[i+1].s2-sj)+gir[i+1].x2*(sj-gir[i+1].s1))/
                  (gir[i+1].s2-gir[i+1].s1) + edxrms*normranf*efrac;
      gir[i].z2:= (gir[i+1].z1*(gir[i+1].s2-sj)+gir[i+1].z2*(sj-gir[i+1].s1))/
                  (gir[i+1].s2-gir[i+1].s1) + edzrms*normranf*efrac;
    end;
  end;

{
  getstr(str,'girdset.dat'); rewrite_(gout,str.str);
  for i:=1 to igr do begin
    writeln(gout,gir[i].s1:15:3,gir[i].x1*1e6:15:1,gir[i].z1*1e6:15:1);
    writeln(gout,gir[i].s2:15:3,gir[i].x2*1e6:15:1,gir[i].z2*1e6:15:1);
  close(gout);
}
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
                  (gir[i].s2-gir[i].s1) + edxrms*normranf*efrac;
      cell.dS[2]:=(gir[i].z1*(gir[i].s2-sj)+gir[i].z2*(sj-gir[i].s1))/
                  (gir[i].s2-gir[i].s1) + edzrms*normranf*efrac;
      if (cell.elem.pkind = mpole) then begin
        cell.elem.m.porder:=3*cell.elem.m.porder;
        if cell.elem.m.pthick = thick then begin{thick elements}
          cell.elem.m.pbsys[1]:= mp11*cell.elem.m.pirho; {dipole error}
          cell.elem.m.pbsys[2]:= mp12*cell.elem.m.pirho;
          cell.elem.m.pbsys[3]:= mp13*cell.elem.m.pirho;
        end
        else begin
          cell.elem.m.pbsys[dip] := 0.0;
        end;
        cell.elem.m.pbsys[4]:= mp24*cell.elem.m.pbpar[quad];
        cell.elem.m.pbsys[6]:= mp26*cell.elem.m.pbpar[quad]+
                               mp36*cell.elem.m.pbpar[sext];
        cell.elem.m.pbsys[9]:= mp36*cell.elem.m.pbpar[sext];
      end;
 {
      if esetout then begin
        writeln('   ',cell.elem.pname[1],cell.elem.pname[2],cell.elem.pname[3],
          cell.elem.pname[4],'/',gir[i].k[j]:2,': sj=',sj:8:3,' m. ',
          ' dx=',cell.dS[1]*1D6:6:1,' um, dz=',cell.dS[2]*1D6:6:1,' um');
        for io:= 1 to cell.elem.m.porder do begin
          writeln(io:1,'-pole, error : ',cell.elem.m.pbsys[io]:15:8);
        end;
      end;
}
      putelem(ipos,cell);
      write('   ',cell.elem.pname[1],cell.elem.pname[2],cell.elem.pname[3],
        cell.elem.pname[4],'/',gir[i].k[j]:2,': ');
      for io:=dip to cell.elem.m.porder do begin
	mpole_setpb(cell.fnum, cell.knum, io);
        write(cell.elem.m.pbsys[io]:10);
      end;
      writeln;
    end;
  end;
end;

{.....................................................................}

procedure GirAssemble (includemon: boolean);
var itismon: boolean; 
    i, j, numel: integer;
begin
  getelem(globval.cell_nloc, cell);
  circumference:=cell.s;
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
      gir[igr].bset:=0;
      gir[igr].kind:=0;
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
    end
    else if (cell.elem.pkind=marker) THEN BEGIN
      IF (cell.elem.pname[1]='x') then begin
        if cell.elem.pname[2]='1' then gir[igr].bset:=1;
        if cell.elem.pname[2]='2' then gir[igr].bset:=2;
      END;
      IF (cell.elem.pname[1]='k') then begin
        if cell.elem.pname[2]='2' then gir[igr].kind:=2; 
      END;
    END;
  end; 
  if gir[igr].ngr=0 then igr:=igr-1;
end;

procedure GirShowComp;
var 
   i, j, ipos, ieol: integer;
begin
  writeln('Number of girders = ',igr:3); 
  for i:=1 to igr do begin
    writeln;
    writeln('Girder ',i:3,' (bset/kind : ',gir[i].bset:1,'/',gir[i].kind:1
                      ,') from ',gir[i].s1:8:3, ' to ', gir[i].s2:8:3,
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


procedure ascodstat (var xmean, ymean, xsigma, ysigma, xmax, ymax: double);
{supplements the no longer compatible procedure codstat}
var 
   x, y, sx, sy, sx2, sy2: double; 
   nbpm, k, lastpos: integer;
begin
  getcod(0.0,lastpos);
  getglobv_(globval);
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

procedure seterr_in;
var 
  k1, k2, i: integer;
  a: char;

begin
  init;
  getstr(str,'seterr.dat'); reset_(cinf,str.str);
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
  writeln('__ Absolute displacement of girder joints and ends :');
  readln(cinf, gdxrms,gdzrms);
  writeln(gdxrms,gdzrms);
  writeln('__ Relative displacement WITHIN girder joints (joint play):');
  readln(cinf, jdxrms,jdzrms);
  writeln(jdxrms,jdzrms);
  writeln('__ Relative displacement of elements ON a girder:');
  readln(cinf, edxrms,edzrms);
  writeln(edxrms,edzrms);
  writeln('-- Dipole-dipole error: ');
  readln(cinf,mp11);
  writeln(mp11);
  writeln('-- Dipole-quadrupole error: ');
  readln(cinf,mp12);
  writeln(mp12);
  writeln('-- Dipole-sextupole error: ');
  readln(cinf,mp13);
  writeln(mp13);
  writeln('-- Quadrupole-octupole error:');
  readln(cinf,mp24);
  writeln(mp24);
  writeln('-- Quadrupole-dodekapole error: ');
  readln(cinf,mp26);
  writeln(mp26);
  writeln('-- Sextupole-dodekapole error: ');
  readln(cinf,mp36);
  writeln(mp36);
  writeln('-- Sextupole-18pol error: ');
  readln(cinf,mp39);
  writeln(mp39);
  writeln;
  writeln('Number of seeds for calculation of mean rms COD :');
  writeln;
  readln(cinf,maxloo);
  writeln(maxloo);

  if maxloo>maxloomax then maxloo:=maxloomax;

  writeln('rms Girder error  : dx=',gdxrms:5:0,' um, dz=',gdzrms:5:0,' um');
  writeln('rms Joint  error  : dx=',jdxrms:5:0,' um, dz=',jdzrms:5:0,' um');
  writeln('rms Element error : dx=',edxrms:5:0,' um, dz=',edzrms:5:0,' um');
  gdxrms:=gdxrms*1E-6;  gdzrms:=gdzrms*1E-6; 
  jdxrms:=jdxrms*1E-6;  jdzrms:=jdzrms*1E-6; 
  edxrms:=edxrms*1E-6;  edzrms:=edzrms*1E-6; 

  writeln('Extended output on girder/element error setting [0 0]...[1 1] ?');
  readln(cinf,k1,k2);
  gsetout:= k1=1;  esetout:= k2=1;
  close(cinf);
end;

procedure BGirSCod(iloo: integer; gf, jf, ef: double; var xclob, yclob: double);
var
  meanx, meany, xmax, ymax: double;
begin
  GirSet (gf, jf, ef);
  ascodstat(meanx,meany,xclob,yclob,xmax,ymax);
  writeln(iloo:5,' CODmean  : x=',meanx*1000:8:3,' mm; y=',meany*1000:8:3,' mm');
  writeln('      CODsigma : x=',xclob*1000:8:3,' mm; y=',yclob*1000:8:3,' mm');
  writeln('      CODmax   : x=',xmax*1000:8:3,' mm; y=',ymax*1000:8:3,' mm');
  writeln;
end;

procedure seterr;

var
  lastpos, iloo: integer;
  mx, my, xclob, yclob: double;
  memxc, memyc: array[1..maxloomax] of double;

begin
  Ring_gettwiss(true, 0.0);
  getglobv_(globval); printglob;
  writeln;
  writeln('globval.cell_nloc    = ',globval.cell_nloc); 
  writeln('globval.elem_nfam = ',globval.elem_nfam);

  seterr_In;

  iniranf(1); setrancut(2);

  if maxloo>1 then begin
    writeln('Creating plotfile seterr.plt');
    getstr(str,'seterr.plt');
    rewrite_(outf,str.str);
    writeln(outf,gdxrms:20,' ',gdzrms:20,' 0.0');
    writeln(outf,jdxrms:20,' ',jdzrms:20,' 0.0');
    writeln(outf,edxrms:20,' ',edzrms:20,' 0.0');
    writeln(outf,maxloo:20);
  
    mx:=0; my:=0;
    for iloo:=1 to maxloo do begin
      BGirSCod(iloo, 1.0, 1.0, 1.0, xclob, yclob);
      mx:=mx + xclob;
      my:=my + yclob;
      memxc[iloo]:=xclob*1D3;
      memyc[iloo]:=yclob*1D3;
      if maxloo>1 then writeln(outf,memxc[iloo]:20:8,memyc[iloo]:20:8);
    end;
   
    mx:=mx/maxloo; my:=my/maxloo;
    writeln('X mean sigma = ',mx*1000:10:4,' mm');
    writeln('Z mean sigma = ',my*1000:10:4,' mm');
    writeln;
    if maxloo>1 then   close(outf);
  end {maxloo>1}


  else begin {maxloo=0}
  end;


end;
