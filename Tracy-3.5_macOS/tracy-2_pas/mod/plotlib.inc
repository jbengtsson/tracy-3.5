  { General routines for plotting }

  const	conid=1; crtdev=70;

  procedure inigks(wkid : integer; metafile : boolean);

  { Initialize GKS }

  var	i	: integer;
	lasf	: lasftyp;
	str	: tstring;

  begin
    { Open GKS }
    gopks(6, 0);

    { Set connection identifier }
    getstr(str, 'metafile.dat');
    if metafile then guesc050s(wkid, str.len, str.str);

    { Open workstation }
    gopwk(wkid, conid, crtdev);
    { Set deferral state }
    gsds(wkid, gbnil, gsuppd);
    { Activate workstation }
    gacwk(wkid);

    { Set text font (Roman) and precision }
    gstxfp(4, gstrkp);
    { Set aspect source flags to individual }
    for i:=1 to lasfnb do
      lasf[i] := gindiv;
{    gsasf(lasf);}
  end;

  procedure closegks(wkid : integer);

  { Close GKS }

  begin
    { Deactivate workstation }
    gdawk(wkid);
    { Close workstation }
    gclwk(wkid);

    { Close GKS }
    gclks;
  end;


  procedure getminmax(n : integer; var x : graphvect; var xmin, xmax : double);

  { Procedure to find min and max values in an array }

  var	i	: integer;

  begin
    xmin := 1d30; xmax := -1d30;
    for i:=1 to n do
    begin
      xmin := min_(x[i], xmin); xmax := max_(x[i], xmax);
    end;
  end;


  procedure getscale(n : integer; var x, y : graphvect; wn : integer);

  { Get coordinate system with appropiate scale }

  var	xmin, xmax, ymin, ymax		: double;
	rxmin, rxmax, rymin, rymax	: real;

  begin
    getminmax(n, x, xmin, xmax); getminmax(n, y, ymin, ymax);
    rxmin := sngl(xmin); rxmax := sngl(xmax);
    rymin := sngl(ymin); rymax := sngl(ymax);
    cframe(wn, rxmin, rxmax, rymin, rymax);
  end;


  procedure gethead(var str : tstring);

  { Get heading for coordinate system }

  begin
    { problems with VAX GKS }
    gschh(0.03); gstxal(gacent, gabase);
    gtxs(0.5, 1.0-0.06, str.len, str.str);
    gstxal(galeft, gabase);
  end;


  procedure getlabel(var xstr, ystr : tstring);

  { Get heading for coordinate system }

  begin
    gschh(0.02); gstxal(gacent, gabase);
    gtxs(0.5, 0.035, xstr.len, xstr.str);
    gschup( -1.0, 0.0 );
    gtxs(0.05, 0.5, ystr.len, ystr.str);
    gschup( 0.0, 1.0 );
    gstxal(galeft, gabase);
  end;


  procedure plotbars(n : integer; var px, py : graphvect);

  { Plot data as bars }

  var	i		: integer;
{	px1, py1	: array [1..2] of double; for compability}
	px1, py1	: graphvect;

  begin
    for i:=1 to n do
    begin
      px1[1] := px[i]; py1[1] := 0;
      px1[2] := px[i]; py1[2] := py[i];
      gpl(2, px1, py1);
    end;
  end;
