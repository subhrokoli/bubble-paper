% ======================================================================= %
% DRIVER

function Rdot = TwoBodyOptTrapDriver(t, R, a, Force, Stresslet)

% initialize LHS
Rdot = zeros(6,1);

% COMPUTE MOBILITY MATRIX ELEMENTS

% % extract position vector components
% X1 = R(1); Y1 = R(2); Z1 = R(3);
% X2 = R(4); Y2 = R(5); Z2 = R(6);

% compute separation vector
% Rho(1) = X1 - X2; Rho(2) = Y1 - Y2; Rho(3) = Z1 - Z2;
Rho = [R(1) - R(4); R(2) - R(5); R(3) - R(6)];

% compute magnitude of separation vector
rho = norm(Rho);

% compute direction cosines
hatRho = Rho / rho;

% compute Faxen factors for J^VF (includes the prefactor 3a/4\rho)
FaxVF1 = ((3*a)/(4*rho)) * (1 + (2*a^2)/(3*rho^2));
FaxVF2 = ((3*a)/(4*rho)) * (1 - (2*a^2)/(rho^2));

% compute Faxen factors for J^VS (includes the prefactor 3a/4\rho^2)
FaxVS1 = ((3*a)/(4*rho^2)) * ((16*a^2)/(5*rho^2));
FaxVS2 = ((3*a)/(4*rho^2)) * (3 - (8*a^2)/(rho^2));

% assemble two-body J^VF (symmetry allows for 5 independent elements)
J_VF(1, 1) = FaxVF2 * hatRho(1) * hatRho(1) + FaxVF1; 
J_VF(1, 2) = FaxVF2 * hatRho(1) * hatRho(2);
J_VF(1, 3) = FaxVF2 * hatRho(1) * hatRho(3); 
J_VF(2, 2) = FaxVF2 * hatRho(2) * hatRho(2) + FaxVF1; 
J_VF(2, 3) = FaxVF2 * hatRho(2) * hatRho(3); 

% assemble two-body J^SF (symmetry allows for 10 independent elements)
J_VS(1, 1, 1) = - FaxVS2 * hatRho(1) * hatRho(1) * hatRho(1) - FaxVS1 * hatRho(1);
J_VS(1, 1, 2) = - FaxVS2 * hatRho(1) * hatRho(1) * hatRho(2) - FaxVS1 * hatRho(2);
J_VS(1, 1, 3) = - FaxVS2 * hatRho(1) * hatRho(1) * hatRho(3) - FaxVS1 * hatRho(3);
J_VS(1, 2, 2) = - FaxVS2 * hatRho(1) * hatRho(2) * hatRho(2);
J_VS(1, 1, 3) = - FaxVS2 * hatRho(1) * hatRho(1) * hatRho(3);
J_VS(1, 3, 3) = - FaxVS2 * hatRho(1) * hatRho(3) * hatRho(3);
J_VS(2, 2, 2) = - FaxVS2 * hatRho(2) * hatRho(2) * hatRho(2) - FaxVS1 * hatRho(2);
J_VS(2, 2, 3) = - FaxVS2 * hatRho(2) * hatRho(2) * hatRho(3) - FaxVS1 * hatRho(3);
J_VS(2, 3, 3) = - FaxVS2 * hatRho(2) * hatRho(3) * hatRho(3);
J_VS(3, 3, 3) = - FaxVS2 * hatRho(3) * hatRho(3) * hatRho(3) - FaxVS1 * hatRho(3);





 

% ODE system
Rdot(1) = - trapFactor*X1 + Fax1*Rdot(4) ...
    + Fax2*(Rdot(4)*nrhox + Rdot(5)*nrhoy + Rdot(6)*nrhoz)*nrhox;

Rdot(2) = - trapFactor*Y1 + Fax1*Rdot(5) ...
    + Fax2*(Rdot(4)*nrhox + Rdot(5)*nrhoy + Rdot(6)*nrhoz)*nrhoy;

Rdot(3) = - trapFactor*Z1 + Fax1*Rdot(6) ...
    + Fax2*(Rdot(4)*nrhox + Rdot(5)*nrhoy + Rdot(6)*nrhoz)*nrhoz;

Rdot(4) = Fax1*Rdot(1) ...
    + Fax2*(Rdot(1)*nrhox + Rdot(2)*nrhoy + Rdot(3)*nrhoz)*nrhox;

Rdot(5) = Fax1*Rdot(2) ...
    + Fax2*(Rdot(1)*nrhox + Rdot(2)*nrhoy + Rdot(3)*nrhoz)*nrhoy;

Rdot(6) = Fax1*Rdot(3) ...
    + Fax2*(Rdot(1)*nrhox + Rdot(2)*nrhoy + Rdot(3)*nrhoz)*nrhoz;

end