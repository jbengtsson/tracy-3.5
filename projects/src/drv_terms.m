classdef drv_terms_type
    properties (Access = private)
        h_label
        h_scl
        h_c
        h_s
        h
    end
    
    methods
        function obj = drv_terms_type()
            obj.h_label = {};
            obj.h_scl = [];
            obj.h_c = [];
            obj.h_s = [];
            obj.h = [];
        end
        
        function get_h_scl(obj, scl_ksi_1, scl_h_3, scl_h_3_delta, scl_h_4, scl_ksi_2, scl_chi_2, scl_chi_delta_2, scl_ksi_3)
            % Implementation of get_h_scl
        end
        
        function get_h(obj, delta_eps, twoJ, delta, twoJ_delta)
            % Implementation of get_h
        end
        
        function prt_h(obj, outf, str, i0, i1)
            % Implementation of prt_h
        end
        
        function prt_h_abs(obj, outf, str, i0, i1)
            % Implementation of prt_h_abs
        end
        
        function print(obj)
            % Implementation of print
        end
    end
end

function A = propagate_drift(L, A)
    % Drift propagator.
    
    Id = eye(size(A)); % Identity matrix
    M = eye(size(A)); % Identity matrix
    
    M(1, 2) = L * Id(1, 2); % Assuming px_ corresponds to column 2
    M(3, 4) = L * Id(3, 4); % Assuming py_ corresponds to column 4
    A = M * A;
end

function A = propagate_thin_kick(L, rho_inv, b2, A)
    % Thin-kick propagator.
    
    Id = eye(size(A)); % Identity matrix
    M = eye(size(A)); % Identity matrix
    
    M(2, 2) = M(2, 2) - (b2 + (rho_inv^2)) * L * Id(1, 1); % Assuming px_ corresponds to column 2
    M(4, 4) = M(4, 4) + b2 * L * Id(3, 3); % Assuming py_ corresponds to column 4
    M(2, 5) = M(2, 5) + rho_inv * L * Id(6, 6); % Assuming delta_ corresponds to column 5
    
    A = M * A;
end


function propagate_fringe_field(L, rho_inv, phi, A)
    % Dipole harde-edge fringe field propagator.

    Id = eye(4);

    if phi ~= 0
        for k = 1:4
            A(px_) = A(px_) + rho_inv * tan(phi * pi / 180) * A(x_, k) * Id(k);
            A(py_) = A(py_) - rho_inv * tan(phi * pi / 180) * A(y_, k) * Id(k);
        end
    end
end

function [h_c, h_s] = get_quad(L, rho_inv, phi1, phi2, b2, alpha0, beta0, nu0, eta0, etap0, m_x, m_y, n_x, n_y, m)
    % Quad propagator: second order symplectic integrator.
    % Use Simpson's rule to integrate for thick magnets.

    h_c = 0;
    h_s = 0;
    n_step = 3;

    alpha = alpha0;
    beta = beta0;
    nu = nu0;
    eta = eta0;
    etap = etap0;

    A1 = get_A(alpha, beta, eta, etap);
    A1 = get_A_CS(2, A1, dnu);

    propagate_fringe_field(L, rho_inv, phi1, A1);

    h = L / n_step;

    for j = 1:n_step
        get_ab(A1, alpha, beta, dnu, eta, etap);
        nu = nu + dnu;
        A1 = get_A_CS(2, A1, dnu);

        phi = 2 * pi * (n_x * nu(X_) + n_y * nu(Y_));
        r = beta(X_)^(m_x / 2) * beta(Y_)^(m_y / 2) * eta(X_)^(m - 1);
        h_c = h_c + r * cos(phi);
        h_s = h_s - r * sin(phi);

        propagate_drift(h / 2, A1);

        if (rho_inv ~= 0) || (b2 ~= 0)
            propagate_thin_kick(h, rho_inv, b2, A1);
        end

        get_ab(A1, alpha, beta, dnu, eta, etap);
        nu = nu + dnu;
        A1 = get_A_CS(2, A1, dnu);

        phi = 2 * pi * (n_x * nu(X_) + n_y * nu(Y_));
        r = beta(X_)^(m_x / 2) * beta(Y_)^(m_y / 2) * eta(X_)^(m - 1);
        h_c = h_c + 4 * r * cos(phi);
        h_s = h_s - 4 * r * sin(phi);

        propagate_drift(h / 2, A1);

        get_ab(A1, alpha, beta, dnu, eta, etap);
        nu = nu + dnu;
        A1 = get_A_CS(2, A1, dnu);
    end

    h_c = h_c * b2 * L / (6 * n_step);
    h_s = h_s * b2 * L / (6 * n_step);

    propagate_fringe_field(L, rho_inv, phi2, A1);

    get_ab(A1, alpha, beta, dnu, eta, etap);
    nu = nu + dnu;
    A1 = get_A_CS(2, A1, dnu);
end


function get_h_ijklm(i, j, k, l, m, A, nu, c, s)
    k1 = 0;
    alpha = zeros(1, 2);
    beta = zeros(1, 2);
    dnu = zeros(1, 2);
    eta = zeros(1, 2);
    etap = zeros(1, 2);
    phi = 0;
    ampl = 0;

    prt = false;

    [alpha, beta, dnu, eta, etap] = get_ab(A);
    A = get_A_CS(2, A, dnu);
    nu = nu + dnu;

    if prt
        fprintf('  %7.3f %6.3f %5.3f %6.3f %6.3f %7.3f %6.3f %5.3f\n', ...
            alpha(X_), beta(X_), nu(X_), eta(X_), etap(X_), ...
            alpha(Y_), beta(Y_), nu(Y_));
    end

    phi = 2 * pi * ((i - j) * nu(X_) + (k - l) * nu(Y_));
    ampl = beta(X_)^((i + j) / 2) * beta(Y_)^((k + l) / 2);
    
    if m >= 1
        ampl = ampl * -eta(X_)^m;
        if m == 1
            ampl = ampl * 2;
        end
    end

    c = c + ampl * cos(phi);
    s = s - ampl * sin(phi);
end

function n = ind(k)
    n = 0;

    if (0 <= k) && (k <= globval.Cell_nLoc)
        n = k;
    elseif k == -1
        n = globval.Cell_nLoc;
    else
        fprintf('ind: incorrect index %ld\n', k);
    end
end

function get_twiss(n, alpha, beta, eta, etap, nu)
    for k = 1:2
        alpha(k) = Cell(n).Alpha(k);
        beta(k) = Cell(n).Beta(k);
        nu(k) = Cell(n).Nu(k);
        eta(k) = Cell(n).Eta(k);
        etap(k) = Cell(n).Etap(k);
    end
end

function get_lin_map(n_step, delta, Elem, M1, M2, M)
    if Elem.PL ~= 0
        M1 = get_edge_lin_map(Elem.M.Pirho, Elem.M.PTx1, Elem.M.Pgap, delta);
        M2 = get_edge_lin_map(Elem.M.Pirho, Elem.M.PTx2, Elem.M.Pgap, delta);
        M = get_sbend_lin_map(Elem.PL / n_step, Elem.M.Pirho, Elem.M.PB(Quad + HOMmax), delta);
    else
        M1 = eye(size(M1)); % Identity matrix
        M2 = eye(size(M2)); % Identity matrix
        M = get_thin_kick_lin_map(Elem.PL * Elem.M.PB(Quad + HOMmax), delta);
    end
end

function get_mult(loc, n_step, delta, i, j, k, l, m, h_c, h_s)
    k1 = 0;
    nu = zeros(1, 2);
    M1 = zeros(2); % Initialize M1
    M2 = zeros(2); % Initialize M2
    M = zeros(2);  % Initialize M
    A = zeros(2);  % Initialize A

    prt = false;

    A = get_A(Cell(ind(loc - 1)).Alpha, Cell(ind(loc - 1)).Beta, ...
              Cell(ind(loc - 1)).Eta, Cell(ind(loc - 1)).Etap);
    nu = Cell(ind(loc - 1)).Nu;

    [M1, M2, M] = get_lin_map(n_step - 1, delta, Cell(loc).Elem);

    if prt
        fprintf('\n');
    end
    h_c = 0;
    h_s = 0;
    get_h_ijklm(i, j, k, l, m, A, nu, h_c, h_s);
    A = M1 * A;
    
    for k1 = 1:n_step - 2
        A = M * A;
        get_h_ijklm(i, j, k, l, m, A, nu, h_c, h_s);
    end
    
    A = M2 * M * A;
    get_h_ijklm(i, j, k, l, m, A, nu, h_c, h_s);
    h_c = h_c / n_step;
    h_s = h_s / n_step;
end


function sxt_1(scl, twoJ, delta, i, j, k, l, m, incl_b3)
    % First order generators.

    n = 0;
    h_c = 0;
    h_s = 0;
    m_x = i + j;
    m_y = k + l;
    n_x = i - j;
    n_y = k - l;
    scl1 = scl * (twoJ(1)^(m_x/2)) * (twoJ(2)^(m_y/2)) * (delta^m);
    
    for n = 0:globval.Cell_nLoc
        if Cell(n + 1).Elem.Pkind == Mpole
            if (Cell(n + 1).Elem.M.Pirho ~= 0 || Cell(n + 1).Elem.M.Porder == Quad) && (m >= 1)
                [b2, a2] = get_bn_design_elem(Cell(n + 1).Fnum, 1, Quad);
                [alpha0, beta0, eta0, etap0, nu0] = get_twiss(ind(n), alpha0, beta0, eta0, etap0, nu0);
                [c, s] = get_quad(Cell(n + 1).Elem.PL, Cell(n + 1).Elem.M.Pirho, ...
                                  Cell(n + 1).Elem.M.PTx1, Cell(n + 1).Elem.M.PTx2, ...
                                  b2, alpha0, beta0, nu0, eta0, etap0, ...
                                  m_x, m_y, n_x, n_y, m);
                h_c = h_c + scl1 * c; 
                h_s = h_s + scl1 * s;
            elseif incl_b3 && (Cell(n + 1).Elem.M.Porder >= Sext)
                [b3L, a3L] = get_bnL_design_elem(Cell(n + 1).Fnum, 1, Sext);
                [c, s] = get_mult(n, 3, 0, i, j, k, l, m);
                h_c = h_c + scl1 * b3L * c; 
                h_s = h_s + scl1 * b3L * s;
            end
        end
    end
end

function sxt_2(scl, i1, j1, k1, l1, i2, j2, k2, l2)
    % Second order generators.

    h_c = 0;
    h_s = 0;
    m_x = [i1 + j1, i2 + j2];
    m_y = [k1 + l1, k2 + l2];
    n_x = [i1 - j1, i2 - j2];
    n_y = [k1 - l1, k2 - l2];

    for n1 = 0:globval.Cell_nLoc
        if (Cell(n1 + 1).Elem.Pkind == Mpole) && (Cell(n1 + 1).Elem.M.Porder >= Sext)
            [b3L1, a3L1] = get_bnL_design_elem(Cell(n1 + 1).Fnum, 1, Sext);
            cp = Cell(ind(n1));
            dnu = n_x(1) * cp.Nu(1) + n_y(1) * cp.Nu(2);
            A1 = scl * b3L1 * (cp.Beta(1)^(m_x(1)/2)) * (cp.Beta(2)^(m_y(1)/2));
            c = 0; 
            s = 0;
            for n2 = 0:globval.Cell_nLoc
                if (Cell(n2 + 1).Elem.Pkind == Mpole) && (Cell(n2 + 1).Elem.M.Porder >= Sext)
                    [b3L2, a3L2] = get_bnL_design_elem(Cell(n2 + 1).Fnum, 1, Sext);
                    cp = Cell(ind(n2));
                    phi = 2 * pi * (dnu + n_x(2) * cp.Nu(1) + n_y(2) * cp.Nu(2));
                    A2 = b3L2 * (cp.Beta(1)^(m_x(2)/2)) * (cp.Beta(2)^(m_y(2)/2));
                    if n2 < n1
                        s = s + A2 * cos(phi); 
                        c = c + A2 * sin(phi);
                    elseif n2 > n1
                        s = s - A2 * cos(phi); 
                        c = c - A2 * sin(phi);
                    end
                end
            end
            h_s = h_s + A1 * s;
            h_c = h_c + A1 * c;
        end
    end
end


