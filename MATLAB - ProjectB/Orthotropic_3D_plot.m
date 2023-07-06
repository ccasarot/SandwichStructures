%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%%%%%%%%%%%  41517 Stiffened Plates and Sandwich Constructions  %%%%%%%%%%%
%%%%%%%%%%%%%%%   DTU - TECHNICAL UNIVERSITY OF DENMARK    %%%%%%%%%%%%%%%%
%%%%%%%%%%%%     Project B: Design of a helicopter floor      %%%%%%%%%%%%%
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%   Copenhagen, Spring semester 2023   %%%%%%%%%%%
%                                                                         %
%                               Christian Casarotto - s223302             %
%                                    Irene Berganzo - s223230             %
%                                                                         %
%%%%%%%    This file is a copy of the orthotropic file, modified    %%%%%%%
%%%%%%%    with the only purpose to create 3D plots of stress on    %%%%%%%
%%%%%%%    the face sheets and tau on the core.                     %%%%%%%
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

% Collection of dimensions of plate panels from the hellicopter floor:
% a_dimension = [2060 2060 1630 670 980 1630 1030 1030 2320];
% b_dimension = [650  650  825  650 670 825  870  870  780];

%% Input data

% % % % % % % % % % % % PLATE PARAMETERS % % % % % % % % % % % % % % % % % 

a = 2060;         % [mm] Length for edge a
b = 650;          % [mm] Length for edge b
t_1 = 4;          % [mm] Thickness for face sheet 1
t_2 = t_1;        % [mm] Thickness for face sheet 2
t_c = 22;         % [mm] Thickness of the core
d   = t_c+t_1;    % [mm] distance between axis of face sheets

% % % % % % % % % % % % % % % % LOAD % % % % % % % % % % % % % % % % % % % 

q_dist = 35.3/1000;         % [MPa] Uniformely distributed load
P_local = 100;              % [N] Local load
A_local = 100;              % [mm^2] Area where the local load acts
q_conc = P_local/A_local;   % [MPa] concentrated load

% % % % % % % % % % % % MATERIAL PARAMETERS % % % % % % % % % % % % % % % % 

% Upper sheet (1)           % Lower sheet (2)
E_x1  = 4.0593e+04;         E_x2  = 4.0593e+04;   % [MPa] Young modulus along x
E_y1  = E_x1;               E_y2  = E_x2;         % [MPa] Young modulus along y
G_xy1 = 1.5185e+04;         G_xy2 = 1.5185e+04;   % [MPa] Shear modulus 
% (metallic materials E_x1/(2*(1+0.3)) -> you might need to change it

v = 0.3366; % If metals (if non metallic face sheets, approx. are later given)   
% it is the same in the other direction  

% Core
G_cx = 320.5;      % [MPa] Shear modulus core in x
G_cy = 155;        % [MPa] Shear modulus core in y
E_cx = 0.305;      % [MPa] mostly used to check for weak core conditon
E_cy = 0.7395;     % [MPa] mostly used to check for weak core conditon
E_c  = 552;        % [MPa] (z!) mostly used to check for weak core conditon

% Compute the integration extremes - q_mn conc
ext_1_min = a/2 - sqrt(A_local)/2; % The are where the load is applied is 100x100 mm
ext_1_max = a/2 + sqrt(A_local)/2; % and centered in the plate
ext_2_min = b/2 - sqrt(A_local)/2;
ext_2_max = b/2 + sqrt(A_local)/2;

%% TYPE OF MATERIAL - CONTROL PANEL

% FACE SHEETS
Metallic_Face_Sheets = 'n';

%% Calculating rigidities

% Stiffnesses
D_x = (E_x1*t_1*E_x2*t_2*d^2) / (E_x1*t_1 + E_x2*t_2);
D_y = (E_y1*t_1*E_y2*t_2*d^2) / (E_y1*t_1 + E_y2*t_2);

% Orthotropic
D_xy = (2*G_xy1*t_1*G_xy2*t_2*d^2) / (G_xy1*t_1 + G_xy2*t_2);

S_x = (G_cx*d^2) / t_c;   % Shear stiffness along x of the core
S_y = (G_cy*d^2) / t_c;   % Shear stiffness along y of the core

if Metallic_Face_Sheets == 'y'
    v_xy=v;
    v_yx=v;
else
    % Poisson ratios of the whole plate, approximations!
    v_xy = 0.25*sqrt(D_x/D_y);
    v_yx = v_xy*(D_x/D_y);
end

%% plot control panel

% define the number of iterations and number of points per edge
limit = 20;
np = 25;

% Creating grid of x and y coordinates
Cx1 = (linspace(0,a,np));
Cx2 = (linspace(0,b,np))';



% PROGRAM STARTS



%% Check for thin faces and weak core condtion

% Thin faces condition
if 3*(d/t_1)^2 < 100 || d/t_1 < 5.77 || 3*(d/t_2)^2 < 100 || d/t_2 < 5.77
    fprintf('You are NOT in thin faces condition!')
    return
end
% Weak core condition
if (6*E_x1*t_1*d^2)/(E_c*t_c^3) < 100 || (6*E_x2*t_2*d^2)/(E_c*t_c^3) < 100
    fprintf('You are NOT in weak core condition!')
    return
end

%% Calculating parameters: X_mn, Y_mn, Z_mn and W_mn

syms m n
syms x y

X_mn = (1/S_y) * ((1/2)*((m*pi/a)^5)*(D_x*D_xy)/(1 - v_xy*v_yx) + ...
    ((m*pi/a)^3)*((n*pi/b)^2)*((D_x*D_y)/(1 - v_xy*v_yx) - ...
    (D_xy*(v_xy*D_x + v_yx*D_y))/(2*(1 - v_xy*v_yx))) + ...
    (1/2)*((m*pi/a)*(n*pi/b)^4)*(D_y*D_xy)/(1 - v_xy*v_yx) ...
    + S_y*(m*pi/a)*(((m*pi/a)^2)*(D_x/(1 - v_xy*v_yx)) + ...
    ((n*pi/b)^2)*(D_xy + (v_yx*D_x)/(1 - v_xy*v_yx))));

Y_mn = (1/S_x) * (-(1/2)*((n*pi/b)^5)*(D_y*D_xy)/(1 - v_xy*v_yx) - ...
    ((m*pi/a)^2)*((n*pi/b)^3)*((D_x*D_y)/(1 - v_xy*v_yx) - ...
    (D_xy*(v_xy*D_x + v_yx*D_y))/(2*(1 - v_xy*v_yx))) - ...
    (1/2)*((m*pi/a)^4*(n*pi/b))*(D_x*D_xy)/(1 - v_xy*v_yx) ...
    - S_x*(n*pi/b)*(((n*pi/b)^2)*(D_y/(1 - v_xy*v_yx)) + ...
    ((m*pi/a)^2)*(D_xy + (v_xy*D_y)/(1 - v_xy*v_yx))));

Z_mn = (m*pi/a)*X_mn - (n*pi/b)*Y_mn;

W_mn = -(1/(S_x*S_y))*((1/2)*D_xy*((m*pi/a)^4*(D_x/(1-v_xy*v_yx)) - ...
    (m*pi/a)^2*(n*pi/b)^2*(v_xy*D_x+v_yx*D_y)/(1-v_xy*v_yx) + ...
    (n*pi/b)^4*(D_y/(1-v_xy*v_yx))) + (m*pi/a)^2*(n*pi/b)^2*D_x*D_y/(1-v_xy*v_yx)) - ...
    ((m*pi/a)^2*D_x/(S_x*(1-v_xy*v_yx)) + (n*pi/b)^2*D_y/(S_y*(1-v_xy*v_yx)))...
    - (1/2)*D_xy*((m*pi/a)^2/S_y + (n*pi/b)^2/S_x) - 1;

% Equation for q_mn dist
q_mn_dist = 4/(a*b) * int(int(q_dist*sin((m*pi*x)/a)*sin((n*pi*y)/b),y,0,b),x,0,a);
% Compute q_mn conc
q_mn_conc = 4/(a*b) * int(int(q_conc*sin((m*pi*x)/a)*sin((n*pi*y)/b),y,ext_2_min,ext_2_max),x,ext_1_min,ext_1_max);
% Equation for q_mn
q_mn = q_mn_dist + q_mn_conc;



%% Computing deflection w

syms m n
syms x y

% Equation for deflections
w_equ = - ((W_mn*q_mn)/Z_mn) * sin(m*pi*x/a) * sin(n*pi*y/b);

% Deflection
Deflection=0;
Deflection1D = subs(w_equ,[x],[Cx1]);
Deflection2D = subs(Deflection1D,[y],[Cx2]);
Deflection2D = matlabFunction(Deflection2D);
for j=1:1:limit
for i=1:1:limit
    m = i; n = j;
    Deflection = Deflection + Deflection2D(m,n);
end
end
Deflection = double(Deflection);

colormap(figure,winter);
surf(Cx1,Cx2,Deflection)
title('Deflection')
xlabel('x [mm]'), ylabel('y [mm]'), zlabel('w [mm]')
ax = gca; 
ax.FontSize = 14; 
colorbar



%% Computing moments

syms m n
syms x y

% Equation for moments
M_x_equ = (-D_x/(1 - v_xy*v_yx)) * (W_mn * ((m*pi/a)^2 + v_yx*(n*pi/b)^2) ...
    + (m*pi/a)*(X_mn/S_x) - v_yx*(n*pi/b)*(Y_mn/S_y)) * (q_mn/Z_mn) * ...
    sin(m*pi*x/a) * sin(n*pi*y/b);
M_y_equ = (-D_y/(1 - v_xy*v_yx)) * (W_mn * ((n*pi/b)^2 + v_xy*(m*pi/a)^2) ...
    - (n*pi/b)*(Y_mn/S_y) + v_xy*(m*pi/a)*(X_mn/S_x)) * (q_mn/Z_mn) * ...
    sin(m*pi*x/a) * sin(n*pi*y/b);

% equations for tensions
sigma_fx1_equ = (M_x_equ*E_x1*E_x2*t_2*d)/(D_x*(E_x1*t_1+E_x2*t_2));
sigma_fy1_equ = (M_y_equ*E_y1*E_y2*t_2*d)/(D_y*(E_y1*t_1+E_y2*t_2));

% Tension
Tension_x=0; Tension_y=0;
Tension_x_1D = subs(sigma_fx1_equ,[x],[Cx1]);
Tension_y_1D = subs(sigma_fy1_equ,[x],[Cx1]);
Tension_x_2D = subs(Tension_x_1D,[y],[Cx2]);
Tension_y_2D = subs(Tension_y_1D,[y],[Cx2]);
Tension_x_2D = matlabFunction(Tension_x_2D);
Tension_y_2D = matlabFunction(Tension_y_2D);
for j=1:1:limit
for i=1:1:limit
    m = i; n = j;
    Tension_x = Tension_x + Tension_x_2D(m,n);
    Tension_y = Tension_y + Tension_y_2D(m,n);
end
end
Tension_x = double(Tension_x);
Tension_y = double(Tension_y);

colormap(figure,winter);
surf(Cx1,Cx2,Tension_x)
title('Tension on upper face sheet - \sigma_x')
xlabel('x [mm]'), ylabel('y [mm]'), zlabel('Tension \sigma_x [MPa]')
ax = gca; 
ax.FontSize = 14; 
colorbar
colormap(figure,winter);
surf(Cx1,Cx2,Tension_y)
title('Tension on upper face sheet - \sigma_y')
xlabel('x [mm]'), ylabel('y [mm]'), zlabel('Tension \sigma_y [MPa]')
ax = gca; 
ax.FontSize = 14; 
colorbar



%% Computing transverse forces

syms m n
syms x y

% Equation for transverse forces
T_x_equ = (X_mn*q_mn)/Z_mn * cos(m*pi*x/a) * sin(n*pi*y/b);
T_y_equ = -(Y_mn*q_mn)/Z_mn * sin(m*pi*x/a) * cos(n*pi*y/b);

tau_cxz_equ = T_x_equ/d;
tau_cyz_equ = T_y_equ/d;

% Tension
tau_cxz=0; 
tau_cyz=0; 
tau_cxz_1D = subs(tau_cxz_equ,[x],[Cx1]);
tau_cyz_1D = subs(tau_cyz_equ,[x],[Cx1]);
tau_cxz_2D = subs(tau_cxz_1D,[y],[Cx2]);
tau_cyz_2D = subs(tau_cyz_1D,[y],[Cx2]);
tau_cxz_2D = matlabFunction(tau_cxz_2D);
tau_cyz_2D = matlabFunction(tau_cyz_2D);
for j=1:1:limit
for i=1:1:limit
    m = i; n = j;
    tau_cxz = tau_cxz + tau_cxz_2D(m,n);
    tau_cyz = tau_cyz + tau_cyz_2D(m,n);
end
end
tau_cxz = double(tau_cxz);
tau_cyz = double(tau_cyz);

colormap(figure,winter);
surf(Cx1,Cx2,tau_cxz)
title('Shear in the core - \tau_{xz}')
xlabel('x [mm]'), ylabel('y [mm]'), zlabel('Shear \tau_{xz} [MPa]')
ax = gca; 
ax.FontSize = 14; 
colorbar
colormap(figure,winter)
surf(Cx1,Cx2,tau_cyz)
title('Shear in the core - \tau_{yz}')
xlabel('x [mm]'), ylabel('y [mm]'), zlabel('Shear \tau_{yz} [MPa]')
ax = gca; 
ax.FontSize = 14; 
colorbar


