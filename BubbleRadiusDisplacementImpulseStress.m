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

clear all; close all; clc

%% Initial parameters
% ------- Medium Parameters ------------------------ %
p0 = 102200      ; % [Pa]        ambient pressure 
c = 1481         ; % [m/s]       speed of sound 
mu = 0.000954    ; % [Pas]       medium dynamic viscosity
rho = 997.8      ; % [kg/m^3]    water density
Sw = 0.0728      ; % [N/m]       water-air surface tension
rhoGC = 1.225    ; % [kg/m^3]    gas core density

% ------- Bubble parameters ------------------------ %
R0 = 2.95e-6     ; % [m]         
Gamma = 1.4      ; % [-]         polytropic exponent 

% ------- Marmottant parameters for the shell ------ %
chi = 0.15        ; % [N/m]       shell elasticity
kappas = 5e-9    ; % [kg/s]      shell viscosity
sigma_R0 = 0.0   ; % [N/m]       initial surface tension of shell

% ------- Ultrasound transmission parameters ------- %
f = 1000e3       ; % [Hz]        driving frequency
NCP = 40         ; %             number of cycles
pa = 160e3       ; % [Pa]        pressure amplitude

%% Initialization

% Prepare temporal discretisation
NCS = 10; % Number cycles before pulse
NCE = 10; % Number cycles after pulse

NIxC = 400; % Number of intervals per cycle
n = (NCS + NCP + NCE)*NIxC; % Total number of intervals

T_starting_buffer = NCS/f; % time before burst
T_burst = NCP/f; % burst time
T_ending_buffer = NCE/f; % time after burst

dt = (1/f)/NIxC; % time step
t = linspace(-T_starting_buffer,T_burst + T_ending_buffer,n+1); % time values

% Prepare driving pulse
p_NIxC = 400; % Number of intervals per cycle for the pressure pulse
p_n = (NCS + NCP + NCE)*p_NIxC; % number of intervals
p_dt = (1/f)/p_NIxC; % time step
p_t = linspace(-T_starting_buffer,T_burst + T_ending_buffer,p_n+1); % pressure time values

p = zeros(1,p_n+1); % initialize pressure values
i_burst_first = find(p_t >= 0, 1, 'first'); % index for start of burst
i_burst_last = find(p_t <= T_burst, 1, 'last'); % index for end of burst
p(i_burst_first:i_burst_last) = pa.*(sin(2*pi*f*p_t(i_burst_first:i_burst_last))); % pressure values
   
% Multiply the driving pulse by the experiemntally recorded pressure envelope
[t_env, up, lo] = ExperimentalPressureEnvelope;
p(p>0) = p(p>0).*interp1(t_env,up,p_t(p>0),'linear', 0);
p(p<0) = p(p<0).*interp1(t_env,lo,p_t(p<0),'linear', 0);


%% Solve radial dynamics

% Computation time
T_start = -T_starting_buffer; 
T_end = T_burst + T_ending_buffer;

% Initialisation parameters for Zhou model
lambda = 0.0259;
R_g = 287;
Ts = 295;
Lambda = (Gamma-1)*lambda/Gamma;

eta = 7;
y1 = eta^2/(eta^2+eta+1);
y2 = (eta^2 + eta)/(eta^2+eta+1);

alpha1 = 1/3*y1^3;
alpha2 = 1/3*(y2^3-y1^3);
alpha3 = 1/3*(1-y2^3);

c1 = y1^4/(4*alpha1);
c2 = (y2^4-y1^4)/(4*alpha2);
c3 = (1 - y2^4)/(4*alpha3);

beta1 = c2 - c1;
beta2 = c3 - c2;
betas = 1 - c3;

M0 = (p0+2*sigma_R0/R0)/(R_g*Ts)*R0^3/3;
M1 = 3*alpha1*M0;
M3 = 3*alpha3*M0;

% Solving the dynamics
opts = odeset('MaxStep',1e-9);
[t_x,x] = ode15s(@(t_x,x) ModifiedRayleighPlesset_MarmottantZhou(t_x,x,Lambda,y1,y2,alpha1,alpha2,alpha3,beta1,beta2,betas,M0,R_g,Ts,p_t,p,p0,chi,kappas,sigma_R0,R0,Gamma,c,mu,rho,Sw), [T_start T_end], [R0 0 p0+2*sigma_R0/R0 M1 M3],opts);
R = real(interp1(t_x, x(:,1), t));
R_dot = real(interp1(t_x, x(:,2), t));
p_g = real(interp1(t_x, x(:,3), t));
M3 = real(interp1(t_x, x(:,5), t));

% Retrieve R_dotdot
for i = 1:length(R)

   rho3 = M3(i)/(alpha3*R(i)^3);
   T3 = p_g(i)/(rho3*R_g);
   Gs = (Ts - T3)/(betas*R(i));
   
   Rbuck = R0/sqrt((sigma_R0/chi)+1); % buckling radius
   Rrupt = Rbuck*(1+(Sw)/chi)^(1/2);  % rupture radius
   
   if R(i) < Rbuck % Buckled state
      sigma_R = 0;
   elseif R(i) >= Rbuck && R(i) <= Rrupt % Elastic regime
      sigma_R = chi*((R(i)^2/Rbuck^2)-1);
   elseif R(i) > Rrupt % Ruptured state
      sigma_R = Sw; 
   else
      sigma_R = 0;
   end
   
   R_dotdot(i) = (((  p_g(i) + R(i)/c*(3*Gamma/R(i)*(Lambda*Gs - p_g(i)*R_dot(i)))   -(2*sigma_R/R(i))-(4*mu*(R_dot(i)/R(i)))-(4*kappas*(R_dot(i)/(R(i)^2)))-(p0+interp1(p_t,p,t(i))))/rho)-(3/2)*R_dot(i)^2)/R(i);

end

%% Solving the translational dynamics

% Initialisation
x = zeros(n+1,2);
dx = zeros(n+1,2);
Re = zeros(n+1,1);
bf = zeros(n+1,1);
amf = zeros(n+1,1);
inf = zeros(n+1,1);
rf = zeros(n+1,1);
df = zeros(n+1,1);
sbf = zeros(n+1,1);
dp1 = zeros(n+1,1);
dp2 = zeros(n+1,1);
cf = zeros(n+1,1);

x(1,1) = 0; % Initial bubble position
x(1,2) = 0; % Initial bubble velocity

t_s = T_start;

L = 12e-6; % Thickness cell layer

beta = 0.8; % Fractional derivative cell
c_beta = 1; % Prefactor fractional unit cell

constant_gradp = 1;

% Solver
for i = 1:n
       Re(i) = 2*R(i)*rho*abs(x(i,2))/mu;
       CD = 24/Re(i) + 6/(1+sqrt(Re(i))) + 0.4;
       if (Re(i) == 0) 
         CD = 0;
       end 
       [~,k] = mink(abs(p_t-(t_s-x(i,1)/c)),2);
       k = sort(k);
       pdot = (p(k(2))-p(k(1)))/p_dt;

       d_beta_x = 0;
       for j = 2:i
           b = (i-j+1)^(1-beta) - (i-j)^(1-beta);
           d_beta_x = d_beta_x + b*(x(j,1) - x(j-1,1));
       end
       d_beta_x = dt^(-beta)/gamma(2-beta)* d_beta_x;

       dx(i,1) = x(i,2);
       dx(i,2) = 1/(0.5*rho*4/3*pi*R(i)^3 + rhoGC*4/3*pi*R0^3)*(4/3*pi*R(i)^3/c*pdot - 0.25*pi*CD*Re(i)*R(i)*mu*x(i,2) - 0.5*rho*x(i,2)*4*pi*R(i)^2*R_dot(i) - 4/3*pi*R(i)^3*rho/(2*(L-x(i,1)))^2*(2*R(i)*R_dot(i)^2 + R(i)^2*R_dotdot(i)) -pi*R(i)^2*c_beta*d_beta_x/L);  

       % Compute displacement
       x(i+1,1) = (x(i,1) + dt * dx(i,1)); 
       x(i+1,2) = (x(i,2) + dt * dx(i,2)); 

       % Compute forces 
       bf(i) = 4/3*pi*R(i)^3/c*pdot;
       amf(i) = -0.5*rho*4/3*pi*R(i)^3*dx(i,2) - 0.5*rho*x(i,2)*4*pi*R(i)^2*R_dot(i);
       inf(i) = -rhoGC*4/3*pi*R0^3*dx(i,2);
       df(i) = -0.25*pi*CD*Re(i)*R(i)*mu*x(i,2);
       sbf(i) = -4/3*pi*R(i)^3*rho/(2*(L-x(i,1)))^2*(2*R(i)*R_dot(i)^2 + R(i)^2*R_dotdot(i));
       rf(i) = -pi*R(i)^2*c_beta*d_beta_x/L;

       % Compute pressure gradients 
       dp1(i) = 1/c*pdot;
       dp2(i) = -rho/(2*(L-x(i,1)))^2*(2*R(i)*R_dot(i)^2 + R(i)^2*R_dotdot(i));
       
       % Compute virtual force for virtual impulse
       cf(i) = -4/3*pi*R(i)^3*constant_gradp;
       
   t_s = t_s + dt;
end

i_burst_first = find(t >= 0, 1, 'first'); % index for start of burst
i_burst_last = find(t <= 10.4e-6, 1, 'last'); % index for end of burst
Ret = mean(2*R(i_burst_first:i_burst_last).*abs(x(i_burst_first:i_burst_last,2))*rho/mu);
Rev = mean(2*R(i_burst_first:i_burst_last).*abs(R_dot(i_burst_first:i_burst_last))*rho/mu);


%% Dimensionless impulses calculation

nv = (R_dot<0);
start1n = strfind([0,nv==1],[0 1]);
end1n = strfind([nv==1,0],[1 0]);
start1n = start1n(3:end);
end1n = end1n(3:end);

pv = (R_dot>0);
start1p = strfind([0,pv==1],[0 1]);
end1p = strfind([pv==1,0],[1 0]);
start1p = start1p(2:end);
end1p = end1p(2:end);

I_bf = zeros(length(start1n),1);
I_sbf = zeros(length(start1n),1);

% Dimensional impulses calculation
for i = 1:length(start1n)
    I_bf(i) = dt*trapz(bf(start1p(i):end1n(i)));
    I_sbf(i) = dt*trapz(sbf(start1p(i):end1n(i)));
end

% Virtual impulse calculation
I_cf = zeros(length(start1n),1);

for i = 1:length(start1n)
    I_cf(i) = dt*trapz(cf(start1p(i):end1n(i)));
end

% Equivalent pressure driver calculation
deltap = (4.789.*R(start1n)'.^4.*sqrt(rho).*constant_gradp./I_cf).^2;

% Dimensionless impulses calculation
XI_bf = I_bf./(4.789.*R(start1n)'.^3.*sqrt(deltap*rho));
XI_sbf = I_sbf./(4.789.*R(start1n)'.^3.*sqrt(deltap*rho));



%% calculate stress from different damage mechanisms

% Jet impact
u_jet = 60;
P_jet = 0.5*rho*c*u_jet; % Water hammer on substrate with characteristics of water

% Oscillation impact
P_pp = 0.5*rho.*sign(R_dot).*R_dot.^2;

% Shear stress
h = sqrt(2*mu/(2*pi*f*rho));
tau_AC = mu*(R_dot)/h;

% Bjerknes pressure
P_bjerknes = (bf + sbf)' ./ (pi.*R.^2);



%% Plots

t_delay = 1.57e-6; % Time between start of the video recording and start of the ultrasound pulse

% Plot time-radius-time evolution
f = figure(1);
f.Position = [0 0 1120 1120/4]; 
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex');
grid on
box on
hold on
plot(t+t_delay,R)
xlim([0,120e-7])
ax=gca; ax.XAxis.Exponent = -6;
xlabel('$t$ [s]','interpreter','latex')
ylabel('$R$ [m]','interpreter','latex')

% Plot pressure driving pulse
f = figure(2);
f.Position = [0 0 1120 1120/4]; 
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex');
grid on
box on
hold on
plot(t+t_delay,p/1000)
xlim([0,120e-7])
ax=gca; ax.XAxis.Exponent = -6;
xlabel('$t$ [s]','interpreter','latex')
ylabel('$p$ [kPa]','interpreter','latex')

% Plot time-displacement evolution
f = figure(3);
f.Position = [0 0 1120 1120/4]; 
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex');
grid on
box on
hold on
plot(t+t_delay,x(:,1))
xlim([0,120e-7])
ylim([0,120e-7])
ax=gca; ax.XAxis.Exponent = -6;
xlabel('$t$ [s]','interpreter','latex')
ylabel('$y$ [m]','interpreter','latex')

% Plot time-pressure gradients evolution 
f = figure(4);
f.Position = [0 0 1120 1120/4]; 
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex');
grid on
box on
hold on
plot(t+t_delay,dp1)
plot(t+t_delay,dp2)
xlim([0,120e-7])
ax=gca; ax.XAxis.Exponent = -6;
ylim([-8e9,8e9])
xlabel('$t$ [s]','interpreter','latex')
ylabel('$\nabla p$ [Pa/m]','interpreter','latex')
legend('From ultrasound', 'From substrate','interpreter','latex','Location','northwest')

% Plot time-dimensionless impulses evolution 
f = figure(5);
f.Position = [0 0 1120 1120/4]; 
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex');
grid on
box on
hold on
tI = 4e-6:1e-6:3e-6+length(XI_bf)*1e-6;
plot(tI,XI_bf,'o')
plot(tI,XI_sbf,'o')
plot(tI,XI_bf+XI_sbf,'o')
xlim([0,120e-7])
ax=gca; ax.XAxis.Exponent = -6;
xlabel('$t$ [s]','interpreter','latex')
ylabel('$\zeta$ [-]','interpreter','latex')
legend('From ultrasound', 'From substrate','Total','interpreter','latex','Location','northwest')


% Plot time-stress evolution
f = figure(6);
f.Position = [0 0 1120 1120/4]; 
set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',20)
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
semilogy(t+t_delay,abs(tau_AC))
hold on
grid on
yline(P_jet)
semilogy(t+t_delay,abs(P_pp))
semilogy(t+t_delay,abs(P_bjerknes))
ylim([1e0, 1e9])
xlim([0,120e-7])
ax=gca; ax.XAxis.Exponent = -6;
xlabel('Time [s]','interpreter','latex')
ylabel('Stress [Pa]','interpreter','latex')
legend('$\tau_{\rm oscillation}$','$p_{\rm jet}$','$p_{\rm oscillation}$','$p_{\rm radiation}$','interpreter','latex','Location','northeast','orientation','horizontal');