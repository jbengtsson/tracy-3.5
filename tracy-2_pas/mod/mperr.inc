type
  mparrtype= array[1..3,1..10] of double;

var
  mparr0: mparrtype;
{  cell: celltype;
  cinf: text;
}
{.....................................................................}

  procedure MPerr_Set (mp: mparrtype);
  var 
    io, ipos: integer;
  begin
    for ipos:=1 to globval.cell_nloc do begin
      getelem(ipos,cell);
      if (cell.elem.pkind = mpole) then begin
        for io:=1 to 10 do begin
          if cell.elem.m.pthick = thick then 
               cell.elem.m.pbsys[io]:= mp[1,io]*cell.elem.m.pirho
          else cell.elem.m.pbsys[io]:= 0.0;
          cell.elem.m.pbsys[io]:= cell.elem.m.pbsys[io]+
                                 mp[2,io]*cell.elem.m.pbpar[quad]+
                                 mp[3,io]*cell.elem.m.pbpar[sext];
          if cell.elem.m.pbsys[io]<>0.0 then cell.elem.m.porder:=io;
        end; 
        putelem(ipos,cell);
        for io:=dip to cell.elem.m.porder do 
          mpole_setpb(cell.fnum, cell.knum, io);
{
        writeln('   ',cell.elem.pname[1],cell.elem.pname[2],cell.elem.pname[3],
           cell.elem.pname[4]);
        for io:=dip to cell.elem.m.porder do 
          writeln('io=',io:2,' ',cell.elem.m.pbsys[io]:10,
                           ' ',cell.elem.m.pbpar[io]:10);
        writeln;
}
      end;
    end;
  end;
 

{.........................................................................}

procedure mperr_in;
var
  i,ipos: integer;
  str: tstring;
begin
  init;
  getstr(str,'mperr.dat'); reset_(cinf,str.str);
  for i:=1 to 10 do  begin
    readln(cinf,mparr0[1,i]);
    writeln(i:2,'-pole/dipole ratio   : ',mparr0[1,i]);
  end;
  for i:=1 to 10 do  begin
    readln(cinf,mparr0[2,i]);
    writeln(i:2,'-pole/quad.pole ratio: ',mparr0[2,i]);
  end;
  for i:=1 to 10 do  begin
    readln(cinf,mparr0[3,i]);
    writeln(i:2,'-pole/quad.pole ratio: ',mparr0[3,i]);
  end;
  writeln;
  close(cinf);
end;
