% ======================================================================= %
% DRIVER

function Xdot = TwoBodyOptTrapDriver(t, x, a, Force, Stresslet)

% initialize LHS
Xdot = zeros(6,1);

% extract position vector components
X1 = x(1); Y1 = x(2); Z1 = x(3);
X2 = x(4); Y2 = x(5); Z2 = x(6);

% compute separation vector
rhox = X1 - X2; rhoy = Y1 - Y2; rhoz = Z1 - Z2;

% compute magnitude of separation vector
rho = sqrt(rhox.*rhox + rhoy.*rhoy + rhoz.*rhoz);

% compute direction cosines
nrhox = rhox./rho; nrhoy = rhoy./rho; nrhoz = rhoz./rho;

% % compute Rn \cdot \hat{\rho} ofr n = 1,2
% R1dotrho = X1*nrhox + Y1*nrhoy + Z1*nrhoz;
% R2dotrho = X2*nrhox + Y2*nrhoy + Z2*nrhoz;

% compute Faxen factors
switch HI
    case 'on'
        preFax = (2*a)/(3*rho);
        Fax1 = preFax * (1 + (2*a^2)/(3*rho^2));
        Fax2 = preFax * (1 - (2*a^2)/(rho^2));
    case 'off'
        Fax1 = 0; Fax2 = 0;
end

% ODE system
Xdot(1) = - trapFactor*X1 + Fax1*Xdot(4) ...
    + Fax2*(Xdot(4)*nrhox + Xdot(5)*nrhoy + Xdot(6)*nrhoz)*nrhox;

Xdot(2) = - trapFactor*Y1 + Fax1*Xdot(5) ...
    + Fax2*(Xdot(4)*nrhox + Xdot(5)*nrhoy + Xdot(6)*nrhoz)*nrhoy;

Xdot(3) = - trapFactor*Z1 + Fax1*Xdot(6) ...
    + Fax2*(Xdot(4)*nrhox + Xdot(5)*nrhoy + Xdot(6)*nrhoz)*nrhoz;

Xdot(4) = Fax1*Xdot(1) ...
    + Fax2*(Xdot(1)*nrhox + Xdot(2)*nrhoy + Xdot(3)*nrhoz)*nrhox;

Xdot(5) = Fax1*Xdot(2) ...
    + Fax2*(Xdot(1)*nrhox + Xdot(2)*nrhoy + Xdot(3)*nrhoz)*nrhoy;

Xdot(6) = Fax1*Xdot(3) ...
    + Fax2*(Xdot(1)*nrhox + Xdot(2)*nrhoy + Xdot(3)*nrhoz)*nrhoz;

end