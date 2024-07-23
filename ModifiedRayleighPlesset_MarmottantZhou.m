% -----------------------------------------------------------------------------
% Copyright (C) 2024 Marco Cattaneo
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the Free
% Software Foundation, either version 3 of the License, or (at your option)
% any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of  MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
% more details.
%
% You should have received a copy of the GNU General Public License along with
% this program.  If not, see <http://www.gnu.org/licenses/>.
% -----------------------------------------------------------------------------

function dxdt = ModifiedRayleighPlesset_MarmottantZhou(t_x,x,Lambda,y1,y2,alpha1,alpha2,alpha3,beta1,beta2,betas,M0,R_g,Ts,p_t,pa,p0,chi,kappas,sigma_R0,R0,gamma,c,mu,rho,Sw)

R(1) = x(1);
R(2) = x(2);
p = x(3);
M1 = x(4);
M3 = x(5);

M2 = M0 - M1 - M3;
rho1 = M1/(alpha1*R(1)^3);
rho2 = M2/(alpha2*R(1)^3);
rho3 = M3/(alpha3*R(1)^3);
T1 = p/(rho1*R_g);
T2 = p/(rho2*R_g);
T3 = p/(rho3*R_g);
G1 = (T2 - T1)/(beta1*R(1));
G2 = (T3 - T2)/(beta2*R(1));
Gs = (Ts - T3)/(betas*R(1));
u1rel = Lambda/p*(G1 - y1*Gs);
u2rel = Lambda/p*(G2 - y2*Gs);

if u1rel >= 0
    rho1uw = rho1;
else
    rho1uw = rho2;
end

if u2rel >= 0
    rho2uw = rho2;
else
    rho2uw = rho3;
end



Rbuck = R0/sqrt((sigma_R0/chi)+1); % buckling radius
Rrupt = Rbuck*(1+(Sw)/chi)^(1/2);  % rupture radius

if R(1) < Rbuck % Buckled state
   sigma_R = 0;
elseif R(1) >= Rbuck && R(1) <= Rrupt % Elastic regime
   sigma_R = chi*((R(1)^2/Rbuck^2)-1);
elseif R(1) > Rrupt % Ruptured state
   sigma_R = Sw;
else
   sigma_R = 0;
end



dxdt  = [R(2);
         (((  p + R(1)/c*(3*gamma/R(1)*(Lambda*Gs - p*R(2)))   -(2*sigma_R/R(1))-(4*mu*(R(2)/R(1)))-(4*kappas*(R(2)/(R(1)^2)))-(p0+interp1(p_t,pa,t_x)))/rho)-(3/2)*R(2)^2)/R(1);
         3*gamma/R(1)*(Lambda*Gs - p*R(2));
         -y1^2*rho1uw*u1rel*R(1)^2;
         y2^2*rho2uw*u2rel*R(1)^2;
        ];



