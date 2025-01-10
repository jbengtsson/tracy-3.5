%{ Module to generate a lattice flat file.

   Type codes:

     marker     -1
     drift	 0
     multipole   1
     cavity      2
     thin kick   3
     wiggler     4
     kick_map    6
     map         7

   Integration methods:

     fourth order symplectic integrator   4

   Format:

     name, family no, kid no, element no
     type code, integration method, no of integration steps, reverse
     apertures: xmin, xmax, ymin, ymax

   The following lines follows depending on element type.

     type

     drift:	 L

     multipole:  hor., ver. displacement, roll angle (design),
                                          roll angle (error)
                 L, 1/rho, entrance angle, exit angle
                 apertures[4]
		 no of nonzero multipole coeff., n design
		 n, b , a
		     n   n
		    ...

     wiggler:    L [m], Lambda [m]
                 no of harmonics
                 harm no, kxV [1/m], BoBrhoV [1/m], kxH, BoBrhoH, phi
                    ...

     cavity:	 L [m], cavity voltage/beam energy [eV], omega/c,
                 beam energy [eV], phi

     thin kick:	 hor., ver. displacement, roll angle (total)
		 no of nonzero multipole coeff.
		 n, b , a
		     n   n
		    ...

     kick_map:   order <file name>

     map:
%}


% Numerical type codes
marker_   = -1;
drift_    =  0;
mpole_    =  1;
cavity_   =  2;
thinkick_ =  3;
wiggler_  =  4;
kick_map_ =  6;
map_      =  7;

function prtName(fp, i, type, method, N, reverse)
    fprintf(fp, '%-15s %4d %4d %4d\n', Cell(i).Elem.PName, Cell(i).Fnum,
	    Cell(i).Knum, i);
    fprintf(fp, ' %3d %3d %3d %4d\n', type, method, N, reverse);
    fprintf(fp, ' %23.16e %23.16e %23.16e %23.16e\n', Cell(i).maxampl(1,1),
	    Cell(i).maxampl(1,2), Cell(i).maxampl(2,1), Cell(i).maxampl(2,2));
end

function prtHOM(fp, n_design, PB, Order)
    nmpole = 0;
    for i = 1:Order
        if (PB(Order-i+1) ~= 0.0) || (PB(Order+i+1) ~= 0.0)
            nmpole = nmpole + 1;
        end
    end
    fprintf(fp, '  %2d %2d\n', nmpole, n_design);
    for i = 1:Order
        if (PB(Order-i+1) ~= 0.0) || (PB(Order+i+1) ~= 0.0)
            fprintf(fp, '%3d %23.16e %23.16e\n', i, PB(Order+i+1),
		    PB(Order-i+1));
        end
    end
end

function prtmfile(mfile_dat)
    n_ps = 6;
    mfile = file_write(mfile_dat);
    for i = 1:globval.Cell_nLoc
        switch Cell(i).Elem.Pkind
            case 'drift'
                prtName(mfile, i, 'drift_', 0, 0, 0);
                fprintf(mfile, ' %23.16e\n', Cell(i).Elem.PL);
            case 'Mpole'
                if Cell(i).Elem.PL ~= 0.0
                    prtName(mfile, i, 'mpole_', Cell(i).Elem.M.Pmethod,
			    Cell(i).Elem.M.PN, Cell(i).Elem.Reverse);
                    fprintf(mfile, ' %23.16e %23.16e %23.16e %23.16e\n',
			    Cell(i).dS(1), Cell(i).dS(2), Cell(i).Elem.M.PdTpar,
			    Cell(i).Elem.M.PdTsys
			    + Cell(i).Elem.M.PdTrms*Cell(i).Elem.M.PdTrnd);
                    fprintf(mfile, ' %23.16e %23.16e %23.16e %23.16e %23.16e\n',
			    Cell(i).Elem.PL, Cell(i).Elem.M.Pirho,
			    Cell(i).Elem.M.PTx1, Cell(i).Elem.M.PTx2,
			    Cell(i).Elem.M.Pgap);
                    prtHOM(mfile, Cell(i).Elem.M.n_design, Cell(i).Elem.M.PB,
			   Cell(i).Elem.M.Porder);
                else
                    prtName(mfile, i, 'thinkick_', Cell(i).Elem.M.Pmethod,
			    Cell(i).Elem.M.PN, Cell(i).Elem.Reverse);
                    fprintf(mfile, ' %23.16e %23.16e %23.16e\n', Cell(i).dS(1),
			    Cell(i).dS(2),
			    Cell(i).Elem.M.PdTsys
			    + Cell(i).Elem.M.PdTrms*Cell(i).Elem.M.PdTrnd);
                    prtHOM(mfile, Cell(i).Elem.M.n_design, Cell(i).Elem.M.PB,
			   Cell(i).Elem.M.Porder);
                end
            case 'Wigl'
                prtName(mfile, i, 'wiggler_', Cell(i).Elem.W.Pmethod,
			Cell(i).Elem.W.PN, Cell(i).Elem.Reverse);
                fprintf(mfile, ' %23.16e %23.16e\n', Cell(i).Elem.PL,
			Cell(i).Elem.W.Lambda);
                fprintf(mfile, '%2d\n', Cell(i).Elem.W.n_harm);
                for j = 1:Cell(i).Elem.W.n_harm
                    fprintf(mfile,
			    '%2d %23.16e %23.16e %23.16e %23.16e %23.16e\n',
			    Cell(i).Elem.W.harm(j), Cell(i).Elem.W.kxV(j),
			    Cell(i).Elem.W.BoBrhoV(j), Cell(i).Elem.W.kxH(j),
			    Cell(i).Elem.W.BoBrhoH(j), Cell(i).Elem.W.phi(j));
                end
            case 'Cavity'
                prtName(mfile, i, 'cavity_', 0, 0, 0);
                fprintf(mfile, ' %23.16e %23.16e %23.16e %d %23.16e %23.16e\n',
			Cell(i).Elem.PL, Cell(i).Elem.C.V_RF
			/(1e9*globval.Energy),
			2.0*pi*Cell(i).Elem.C.f_RF/c0, Cell(i).Elem.C.harm_num,
			1e9*globval.Energy, Cell(i).Elem.C.phi_RF);
            case 'marker'
                prtName(mfile, i, 'marker_', 0, 0, 0);
            case 'Insertion'
                prtName(mfile, i, 'kick_map_', Cell(i).Elem.ID.Pmethod,
			Cell(i).Elem.ID.PN, Cell(i).Elem.Reverse);
                if Cell(i).Elem.ID.firstorder
                    fprintf(mfile, ' %3.1f %1d %s\n', Cell(i).Elem.ID.scaling,
			    1, Cell(i).Elem.ID.fname1);
                else
                    fprintf(mfile, ' %3.1f %1d %s\n', Cell(i).Elem.ID.scaling,
			    2, Cell(i).Elem.ID.fname2);
                end
            case 'Map'
                prtName(mfile, i, 'map_', 0, 0, 0);
                for j = 1:n_ps
                    for k = 1:n_ps
                        fprintf(mfile, ' %23.16e', Cell(i).Elem.Map.M(j,k));
                    end
                    fprintf(mfile, '\n');
                end
            otherwise
                fprintf('prtmfile: unknown type %d\n', Cell(i).Elem.Pkind);
                return;
        end
    end
    fclose(mfile);
end
