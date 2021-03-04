function [x_guess] = inguess_ecos(N, nrev, inguess, auxdata)
%
% inguess_ecos
% 
% FUNCTION DESCRIPTION: this function computes the initial trajectory
% guess.
% 
% INPUTS:
%
% N              Number of nodes
%
% nrev           Number of revolutions
% 
% inguess        Type of initial guess: 1 for Cubic SB, 2 for FFS
%
% OUTPUTS:                          
% 
% x_guess:       Initial guess
% 
% AUTHOR: Andrea Carlo Morelli
% 
% DATE: 25/02/2021
%

% Initial and Final BCs
x0 = auxdata.x0;
xf = auxdata.xf;

% Cubic - based initial guess
[x_guess_SB, x_guess_cyl_SB] = SB_guess(N, x0, xf, nrev, auxdata);

% FFS - based initial guess
if inguess == 2
    ncoeff = 10;
    [x_guess_FFS] = FFS_guess(N, ncoeff, x_guess_cyl_SB, auxdata);
    x_guess_FFS(:,7) = x_guess_SB(:,7);
end

% Choice of the initial guess
if inguess == 1
    x_guess = x_guess_SB;
else
    x_guess = x_guess_FFS;
end

function [x_guess, x_guess_cyl] = SB_guess(m_guess, x0, xf, nrev, auxdata)
%
% SB_guess
% 
% FUNCTION DESCRIPTION: this function computes the initial guess based on
% the shape-based method.
% 
% INPUTS:
% 
% m_guess:      Number of time instants                         [1x1]
%
% x0:           Vector of the initial BCs                       [1x7]
%
% xf:           Vector of the final BCs                         [1x7]
%
% nrev:         Number of revolutions of the initial 
%               guess trajectory                                [1x1]
%
% OUTPUTS:
%
% x_guess:      Matrix of initial guess in cartesian 
%               coordinates computed with the shape method        [m_guess x 11]
% x_guess:      Matrix of initial guess in cylindrical
%               coordinates computed with the shape method        [m_guess x 11]
%
% AUTHOR: Andrea Carlo Morelli
% 
% DATE: 19/09/2020
%

c = auxdata.engine.c;
tf = auxdata.tf;
ve = auxdata.engine.ve;
V0 = auxdata.units.V0;

% Normalized time vector & time step
t_ad = linspace(0,1, m_guess);
t_aux = linspace(0,tf, m_guess);
h = t_aux(2) - t_aux(1);

% Initial and final states
x_0 = x0(1); y_0 = x0(2); vx_0 = x0(4); vy_0 = x0(5);
vw_0 = x0(6); z_0 = x0(7);

x_f = xf(1); y_f = xf(2); vx_f = xf(4); vy_f = xf(5);
vw_f = xf(6);

% Initial and final positions in cylindrical coordinates
[th_0,r_0,w_0] = cart2pol(x0(1),x0(2),x0(3));
[th_f,r_f,w_f] = cart2pol(xf(1),xf(2),xf(3));
th_f = th_f + nrev*2*pi;

% Initial and final velocities in cylindrical coordinates
rdot_0 = tf*(x_0*vx_0 + y_0*vy_0)/sqrt(x_0^2 + y_0^2);
rdot_f = tf*(x_f*vx_f + y_f*vy_f)/sqrt(x_f^2 + y_f^2);

thdot_0 = tf*(x_0*vy_0 - y_0*vx_0)/(x_0^2 + y_0^2);
thdot_f = tf*(x_f*vy_f - y_f*vx_f)/(x_f^2 + y_f^2);

% Coefficients matrix
M = [0 0 0 1;...
    1 1 1 1;...
    0 0 1 0;...
    3 2 1 0];

% Constant terms vectors
br = [r_0 r_f rdot_0 rdot_f]';
bth = [th_0 th_f thdot_0 thdot_f]';
bw = [w_0 w_f vw_0*tf vw_f*tf]';

% Coefficient of the polynomials
coeff_r = M\br;
coeff_th = M\bth;
coeff_w = M\bw;

% Expressions of the spatial coordinates
r = coeff_r(1)*t_ad.^3 + coeff_r(2)*t_ad.^2 + ...
     + coeff_r(3)*t_ad.^2 + coeff_r(4);
th = coeff_th(1)*t_ad.^3 + coeff_th(2)*t_ad.^2 + ...
    coeff_th(3)*t_ad + coeff_th(4);
w = coeff_w(1)*t_ad.^3 + coeff_w(2)*t_ad.^2 + ...
    coeff_w(3)*t_ad + coeff_w(4);

% Expressions of the velocities
rdot = (3*coeff_r(1)*t_ad.^2 + 2*coeff_r(2)*t_ad + ...
    coeff_r(3))/tf;
thdot = (3*coeff_th(1)*t_ad.^2 + 2*coeff_th(2)*t_ad + ...
    coeff_th(3))/tf;
vw = (3*coeff_w(1)*t_ad.^2 + 2*coeff_w(2)*t_ad + ...
    coeff_w(3))/tf;

% Expressions of the accelerations
r2dot = (6*coeff_r(1)*t_ad + 2*coeff_r(2))/tf^2;
th2dot = (6*coeff_th(1)*t_ad + 2*coeff_th(2))/tf^2;
aw = (6*coeff_w(1)*t_ad + 2*coeff_w(2))/tf^2;

% Position in Cartesian coordinates
[xvar,y,w] = pol2cart(th,r,w);

% Velocity in Cartesian coordinates
vx = (rdot.*cos(th) - r.*thdot.*sin(th));
vy = (rdot.*sin(th) + r.*thdot.*cos(th));

% Acceleration in Cartesian coordinates
ax = (r2dot.*cos(th) - rdot.*thdot.*...
    sin(th) - rdot.*thdot.*sin(th) - ...
    r.*th2dot.*sin(th) - r.*thdot.^2.*...
    cos(th));

ay = (r2dot.*sin(th) + rdot.*thdot.*...
    cos(th) + rdot.*thdot.*cos(th) + ...
    r.*th2dot.*cos(th) - r.*thdot.^2.*...
    sin(th));

% Control variables first guess
r_tot = sqrt(xvar.^2 + y.^2 + w.^2);
taux = (ax + xvar./r_tot.^3);
tauy = (ay + y./r_tot.^3);
tauw = (aw + w./r_tot.^3);

tau = sqrt(taux.^2 + tauy.^2 + tauw.^2);

% Mass history first guess
z = zeros(1,m_guess);
z(1) = z_0;
for i = 1 : m_guess-1
    z(i+1) = z(i) + h*(-c/(ve/V0))*tau(i);
end

% Cartesian initial guess
x_guess = [xvar' y' w' vx' vy' vw' z' taux' tauy' tauw' tau'];

% Cylindrical initial guess
x_guess_cyl = [r' th' w' rdot' thdot' vw'];

end
function [x_guess] = FFS_guess(m_guess, ncoeff, x_guess_cyl, auxdata)
%
% FFS_guess
% 
% FUNCTION DESCRIPTION: this function computes the initial guess based on
% the Finite Fourier Series method.
% 
% INPUTS:
% 
% m_guess:      Number of time instants                         [1x1]
%
% ncoeff:       Number of harmonical terms                      [1x1]
%
% x_guess_cyl:  Matrix of initial guess in cylindrical 
%               coordinates computed with a shape-based method  [m_guess x 11]
%
% OUTPUTS:
%
% x_guess:      Matrix of initial guess in cartesian 
%               coordinates computed with the FFS method        [m_guess x 11]
%
% AUTHOR: Andrea Carlo Morelli
% 
% DATE: 19/09/2020
%

tf = auxdata.tf;
ve = auxdata.engine.ve;
V0 = auxdata.units.V0;
Ta_max = auxdata.Ta_max;

% Auxiliary functions definition
function A = A_FFS(n, m, tau)
%
% A_FFS
% 
% FUNCTION DESCRIPTION: Matrix of coefficients for the Finite Fourier
% Series method for finding the initial trajectory guess.
% 
% INPUTS:
% 
% n:    Number of armonic functions to be considered    [1x1]
%
% m:    Number of time instants                         [1x1]
%
% tau:  Vector of nondimensional, scaled time           [1 x m]
%
% OUTPUTS:
%
% A:    Matrix of coefficients
%
% AUTHOR: Andrea Carlo Morelli
% 
% DATE: 19/09/2020
%

A = zeros(m,2*n - 3);

Ca0 = 0.5*(1 - cos(2*pi*tau));
A(:,1) = Ca0;

Can = zeros(m, n - 2);
Cbn = zeros(m, n - 2);
for j = 3 : n
    if mod(j,2) == 0
        Can(:,j-2) = cos(j*pi*tau) - cos(2*pi*tau);
        Cbn(:,j-2) = sin(j*pi*tau) - 0.5*j*sin(2*pi*tau);
    else
        Can(:,j-2) = cos(j*pi*tau) - cos(pi*tau);
        Cbn(:,j-2) = sin(j*pi*tau) - j*sin(pi*tau);
    end
end

for j = 2 : size(A,2)
    
    if mod(j,2) == 0
        A(:,j) = Can(:, 0.5*j);
    else
        A(:,j) = Cbn(:, 0.5*(j - 1));
    end
end

end
function Ap = Ap_FFS(n, m, tau)
%
% A_gauss
% 
% FUNCTION DESCRIPTION: First derivative of the matrix of coefficients 
% for the Finite Fourier Series method for finding the initial trajectory guess.
% 
% INPUTS:
% 
% n:    Number of armonic functions to be considered    [1x1]
%
% m:    Number of time instants                         [1x1]
%
% tau:  Vector of nondimensional, scaled time           [1 x m]
%
% OUTPUTS:
%
% Ap:    Matrix of coefficients
%
% AUTHOR: Andrea Carlo Morelli
% 
% DATE: 19/09/2020
%

Ap = zeros(m,2*n - 3);

Ca0 = pi*sin(2*pi*tau);
Ap(:,1) = Ca0;

Can = zeros(m, n - 2);
Cbn = zeros(m, n - 2);
for j = 3 : n
    if mod(j,2) == 0
        Can(:,j-2) = -j*pi*sin(j*pi*tau) + 2*pi*sin(2*pi*tau);
        Cbn(:,j-2) = j*pi*cos(j*pi*tau) - j*pi*cos(2*pi*tau);
    else
        Can(:,j-2) = -j*pi*sin(j*pi*tau) + pi*sin(pi*tau);
        Cbn(:,j-2) = j*pi*cos(j*pi*tau) - j*pi*cos(pi*tau);
    end
end

for j = 2 : size(Ap,2)
    
    if mod(j,2) == 0
        Ap(:,j) = Can(:, 0.5*j);
    else
        Ap(:,j) = Cbn(:, 0.5*(j - 1));
    end
end

end
function Ap2 = Ap2_FFS(n, m, tau)
%
% A_gauss
% 
% FUNCTION DESCRIPTION: Second derivative of the matrix of coefficients 
% for the Finite Fourier Series method for finding the initial trajectory guess.
% 
% INPUTS:
% 
% n:    Number of armonic functions to be considered    [1x1]
%
% m:    Number of time instants                         [1x1]
%
% tau:  Vector of nondimensional, scaled time           [1 x m]
%
% OUTPUTS:
%
% Ap2:    Matrix of coefficients
%
% AUTHOR: Andrea Carlo Morelli
% 
% DATE: 19/09/2020
%

Ap2 = zeros(m,2*n - 3);

Ca0 = 2*pi^2*cos(2*pi*tau);
Ap2(:,1) = Ca0;

Can = zeros(m, n - 2);
Cbn = zeros(m, n - 2);
for j = 3 : n
    if mod(j,2) == 0
        Can(:,j-2) = -(j*pi)^2*cos(j*pi*tau) + 4*pi^2*cos(2*pi*tau);
        Cbn(:,j-2) = -(j*pi)^2*sin(j*pi*tau) + 2*j*pi^2*sin(2*pi*tau);
    else
        Can(:,j-2) = -(j*pi)^2*cos(j*pi*tau) + pi^2*cos(pi*tau);
        Cbn(:,j-2) = -(j*pi)^2*sin(j*pi*tau) + j*pi^2*sin(pi*tau);
    end
end

for j = 2 : size(Ap2,2)
    
    if mod(j,2) == 0
        Ap2(:,j) = Can(:, 0.5*j);
    else
        Ap2(:,j) = Cbn(:, 0.5*(j - 1));
    end
end

end
function Ta = obj_FFS(x, Ar, Aw, Ar_p, Ath_p, Ar_p2, ...
    Ath_p2, Aw_p2, Fr, Fth_p, Fw, Fr_p, Fr_p2, Fth_p2, Fw_p2)
%
% obj_FFS
%
% FUNCTION DESCRIPTION: this function computes the objective to be
% minimized to find the coefficient for the FFS initial guess.
% 
% INPUTS:
% 
% x      Vector of states   
%
% Ar     Matrix of coefficients realted to the radius
%
% Aw     Matrix of coefficients realted to w
%
% Ar_p   First derivative matrix of coefficients realted to the radius
% 
% Ath_p  First derivative matrix of coefficients realted to theta
%
% Ar_p2  Second derivative matrix of coefficients realted to the radius
%
% Ath_p2 Second derivative matrix of coefficients realted to theta
% 
% Aw_p2  First derivative matrix of coefficients realted to w
%
% Fr     Vector of coefficients realted to the radius
%
% Fth_p  First derivative vector of coefficients realted to theta
%
% Fw     Vector of coefficients realted to w
%
% Fr_p   First derivative vector of coefficients realted to the radius
%
% Fr_p2  Second derivative vector of coefficients realted to the radius
%
% Fth_p2 Second derivative vector of coefficients realted to theta
%
% Fw_p2  Second derivative vector of coefficients realted to w
%
% OUTPUTS:
%
% Ta:    Objective function to be minimized
%
% AUTHOR: Andrea Carlo Morelli
% 
% DATE: 19/09/2020
%

xr = x(:,1);
xth = x(:,2);
xw = x(:,3);
s_FFS = sqrt((Ar*xr + Fr).^2 + (Aw*xw + Fw).^2);

fr = (Ar_p2*xr + Fr_p2)/tf^2 - (Ar*xr + Fr).*((Ath_p*xth + Fth_p)/tf).^2 ...
    + (Ar*xr + Fr)./s_FFS.^3;

ftheta = 2*(Ar_p*xr + Fr_p)/tf + (Ar*xr + Fr).*(Ath_p2*xth + ...
    Fth_p2)/(tf^2);

fw = (Aw_p2*xw + Fw_p2)/tf^2 + (Aw*xw + Fw)./s_FFS.^3;

Ta = sum(sqrt(fr.^2 + ftheta.^2 + fw.^2));

end
function [con ceq] = Ta_FFS(x, Ar, Aw, Ar_p, Ath_p, Ar_p2, ...
    Ath_p2, Aw_p2, Fr, Fth_p, Fw, Fr_p, Fr_p2, Fth_p2, Fw_p2)
%
% obj_FFS
%
% FUNCTION DESCRIPTION: this function defines the inequality constraints to
% be respected while solving the NLP realted to the FFS algorithm.

% INPUTS:
% 
% x      Vector of states   
%
% Ar     Matrix of coefficients realted to the radius
%
% Aw     Matrix of coefficients realted to w
%
% Ar_p   First derivative matrix of coefficients realted to the radius
% 
% Ath_p  First derivative matrix of coefficients realted to theta
%
% Ar_p2  Second derivative matrix of coefficients realted to the radius
%
% Ath_p2 Second derivative matrix of coefficients realted to theta
% 
% Aw_p2  First derivative matrix of coefficients realted to w
%
% Fr     Vector of coefficients realted to the radius
%
% Fth_p  First derivative vector of coefficients realted to theta
%
% Fw     Vector of coefficients realted to w
%
% Fr_p   First derivative vector of coefficients realted to the radius
%
% Fr_p2  Second derivative vector of coefficients realted to the radius
%
% Fth_p2 Second derivative vector of coefficients realted to theta
%
% Fw_p2  Second derivative vector of coefficients realted to w
%
% OUTPUTS:
%
% con:    Inequality constraints vector
% 
% ceq:    Equality constraints vector
%
% AUTHOR: Andrea Carlo Morelli
% 
% DATE: 19/09/2020
%

xr = x(:,1);
xth = x(:,2);
xw = x(:,3);
s = sqrt((Ar*xr + Fr).^2 + (Aw*xw + Fw).^2);

fr = (Ar_p2*xr + Fr_p2)/tf^2 - (Ar*xr + Fr).*((Ath_p*xth + Fth_p)/tf).^2 ...
    + (Ar*xr + Fr)./s.^3;

ftheta = 2*(Ar_p*xr + Fr_p)/tf + (Ar*xr + Fr).*(Ath_p2*xth + ...
    Fth_p2)/(tf^2);

fw = (Aw_p2*xw + Fw_p2)/tf^2 + (Aw*xw + Fw)./s.^3;

con = [];
ceq = [];

n = length(fr);
con = zeros(1,n);

for j = 1 : n
    
    con(j) = sqrt(fr(j).^2 + ftheta(j).^2 + fw(j).^2) - Ta_max;
    
end


end

% Normalized time vector & time step
t_ad = linspace(0,1,m_guess);
t_aux = linspace(0,tf,m_guess);
h = t_aux(2) - t_aux(1);

% BCs settings
r_0 = x_guess_cyl(1,1);
th_0 = x_guess_cyl(1,2);
w_0 = x_guess_cyl(1,3);
rdot_0 = x_guess_cyl(1,4)*tf;
thdot_0 = x_guess_cyl(1,5)*tf;
vw_0 = x_guess_cyl(1,6)*tf;

r_f = x_guess_cyl(end,1);
th_f = x_guess_cyl(end,2);
w_f = x_guess_cyl(end,3);
rdot_f = x_guess_cyl(end,4)*tf;
thdot_f = x_guess_cyl(end,5)*tf;
vw_f = x_guess_cyl(end,6)*tf;

% State first guess from SB approach
r_guess = x_guess_cyl(:,1);
th_guess = x_guess_cyl(:,2);
w_guess = x_guess_cyl(:,3);

% Fr, Ftheta and  Fw

% Fr
Fr = (0.5*(r_0 - r_f)*cos(pi*t_ad) + 1/(2*pi)*(rdot_0 - rdot_f)* ...
    sin(pi*t_ad) + 0.5*(r_0 + r_f)*cos(2*pi*t_ad) + ...
    1/(4*pi)*(rdot_0 + rdot_f)*sin(2*pi*t_ad))';

Fr_p = (-0.5*(r_0 - r_f)*pi*sin(pi*t_ad) + 0.5*(rdot_0 - rdot_f)* ...
    cos(pi*t_ad) - (r_0 + r_f)*pi*sin(2*pi*t_ad) + ...
    0.5*(rdot_0 + rdot_f)*cos(2*pi*t_ad))';

Fr_p2 = (-0.5*pi^2*(r_0 - r_f)*cos(pi*t_ad) - 0.5*pi*(rdot_0 - rdot_f)* ...
    sin(pi*t_ad) - 2*pi^2*(r_0 + r_f)*cos(2*pi*t_ad) - ...
    pi*(rdot_0 + rdot_f)*sin(2*pi*t_ad))';

% Ftheta
Fth = (0.5*(th_0 - th_f)*cos(pi*t_ad) + 1/(2*pi)*(thdot_0 - thdot_f)* ...
    sin(pi*t_ad) + 0.5*(th_0 + th_f)*cos(2*pi*t_ad) + ...
    1/(4*pi)*(thdot_0 + thdot_f)*sin(2*pi*t_ad))';

Fth_p = (-0.5*(th_0 - th_f)*pi*sin(pi*t_ad) + 0.5*(thdot_0 - thdot_f)* ...
    cos(pi*t_ad) - (th_0 + th_f)*pi*sin(2*pi*t_ad) + ...
    0.5*(thdot_0 + thdot_f)*cos(2*pi*t_ad))';

Fth_p2 = (-0.5*pi^2*(th_0 - th_f)*cos(pi*t_ad) - 0.5*pi*(thdot_0 - thdot_f)* ...
    sin(pi*t_ad) - 2*pi^2*(th_0 + th_f)*cos(2*pi*t_ad) - ...
    pi*(thdot_0 + thdot_f)*sin(2*pi*t_ad))';

% Fw
Fw = (0.5*(w_0 - w_f)*cos(pi*t_ad) + 1/(2*pi)*(vw_0 - vw_f)* ...
    sin(pi*t_ad) + 0.5*(w_0 + w_f)*cos(2*pi*t_ad) + ...
    1/(4*pi)*(vw_0 + vw_f)*sin(2*pi*t_ad))';

Fw_p = (-0.5*(w_0 - w_f)*pi*sin(pi*t_ad) + 0.5*(vw_0 - vw_f)* ...
    cos(pi*t_ad) - (w_0 + w_f)*pi*sin(2*pi*t_ad) + ...
    0.5*(vw_0 + vw_f)*cos(2*pi*t_ad))';

Fw_p2 = (-0.5*pi^2*(w_0 - w_f)*cos(pi*t_ad) - 0.5*pi*(vw_0 - vw_f)* ...
    sin(pi*t_ad) - 2*pi^2*(w_0 + w_f)*cos(2*pi*t_ad) - ...
    pi*(vw_0 + vw_f)*sin(2*pi*t_ad))';

% A matrices
Ar = A_FFS(ncoeff, m_guess,t_ad);
Ar_p = Ap_FFS(ncoeff, m_guess,t_ad);
Ar_p2 = Ap2_FFS(ncoeff, m_guess,t_ad);

Ath = A_FFS(ncoeff, m_guess,t_ad);
Ath_p = Ap_FFS(ncoeff, m_guess,t_ad);
Ath_p2 = Ap2_FFS(ncoeff, m_guess,t_ad);

Aw = A_FFS(ncoeff, m_guess,t_ad);
Aw_p = Ap_FFS(ncoeff, m_guess,t_ad);
Aw_p2 = Ap2_FFS(ncoeff, m_guess,t_ad);

% Initial guess for the coefficients
Xr_guess = Ar\(r_guess - Fr);
Xth_guess = Ath\(th_guess - Fth);
Xw_guess = Aw\(w_guess - Fw);

x0_FFS = [Xr_guess Xth_guess Xw_guess];

% FMINCON 
options = optimset('Algorithm','interior-point' , 'LargeScale', 'on' ,...
'Display' , 'iter' ,...
'TolX' , 1e-5 , 'TolFun' , 1e-5 , 'TolCon' , 1e-5 ,...
'MaxIter', 5e3 , 'MaxFunEvals' , 5e6);

optimal = fmincon(@(x)obj_FFS(x, Ar, Aw, Ar_p, Ath_p, Ar_p2, ...
    Ath_p2, Aw_p2, Fr, Fth_p, Fw, Fr_p, Fr_p2, Fth_p2, Fw_p2) , x0_FFS,...
[],[],[],[],[],[],@(x)Ta_FFS(x,  Ar, Aw, Ar_p, Ath_p, Ar_p2, ...
    Ath_p2, Aw_p2, Fr, Fth_p, Fw, Fr_p, Fr_p2, Fth_p2, Fw_p2),options);

% Initial guess reconstruction

% Positions
r = Ar*optimal(:,1) + Fr;
th = Ath*optimal(:,2) + Fth;
w = Aw*optimal(:,3) + Fw;

% Velocities
rdot = (Ar_p*optimal(:,1) + Fr_p)/tf;
thdot = (Ath_p*optimal(:,2) + Fth_p)/tf;
wdot = (Aw_p*optimal(:,3) + Fw_p)/tf;

% Accelerations
r2dot = (Ar_p2*optimal(:,1) + Fr_p2)/tf^2;
th2dot = (Ath_p2*optimal(:,2) + Fth_p2)/tf^2;
w2dot = (Aw_p2*optimal(:,3) + Fw_p2)/tf^2;

% Position in Cartesian coordinates
[xvar,y,w] = pol2cart(th,r,w);

% Velocity in Cartesian coordinates
vx = (rdot.*cos(th) - ...
    r.*thdot.*sin(th));
vy = (rdot.*sin(th) + ...
    r.*thdot.*cos(th));
vw = wdot;

% Acceleration in Cartesian coordinates
ax = (r2dot.*cos(th) - rdot.*thdot.*...
    sin(th) - rdot.*thdot.*sin(th) - ...
    r.*th2dot.*sin(th) - r.*thdot.^2.*...
    cos(th));
ay = (r2dot.*sin(th) + rdot.*thdot.*...
    cos(th) + rdot.*thdot.*cos(th) + ...
    r.*th2dot.*cos(th) - r.*thdot.^2.*...
    sin(th));
aw = w2dot;

% % Control variables first guess
 s_tau = sqrt(r.^2 + w.^2);

Tr = (r2dot - r.*(thdot).^2 + r./s_tau.^3);
Tth = (2*rdot + r.*th2dot);
Tw = (w2dot + w./s_tau.^3);

Tx = Tr.*cos(th) - Tth.*sin(th);
Ty = Tr.*sin(th) + Tth.*cos(th);

taux = Tx;
tauy = Ty;
tauw = Tw;

tau = sqrt(taux.^2 + tauy.^2 + tauw.^2);

z = zeros(1,m_guess);
z(1) = 0;
for i = 1 : m_guess - 1
    z(i+1) = z(i) + h*(-1/(ve/V0))*tau(i);
end

x_guess = [xvar y w vx vy vw z' taux tauy tauw tau];

end

% First guess plot
figure
hold all
plot3(x_guess(:,1),x_guess(:,2),x_guess(:,3),'b','LineWidth',1)
xlabel('[AU]')
ylabel('[AU]')
zlabel('[AU]')
grid on


end

