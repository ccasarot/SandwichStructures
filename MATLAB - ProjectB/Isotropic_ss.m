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
%%%%%%%     First file for isotropic simply supported sandwich      %%%%%%%
%%%%%%%     panels. Might not be perfect as the focus moved         %%%%%%%
%%%%%%%     toward orthotropic case                                 %%%%%%%
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

%% Input data

% % % % % % % % % % % % PLATE PARAMETERS % % % % % % % % % % % % % % % % % 

a = 2060;        % [mm] Length for edge a
b = 650;         % [mm] Length for edge b
t_1 = 2;         % [mm] Thickness for face sheet 1
t_2 = 2;         % [mm] Thickness for face sheet 2
t_c = 26;        % [mm] Thickness of the core
d   = 28;        % [mm] distance between axis of face sheets

% % % % % % % % % % % % % % % % LOAD % % % % % % % % % % % % % % % % % % % 

q = 800*9.81/1000000;    % [MPa] Uniformely distributed load

% % % % % % % % % % % % MATERIAL PARAMETERS % % % % % % % % % % % % % % % % 

% Upper sheet (1)           % Lower sheet (2)
E_x1 = 70000;               E_x2 = 70000;   % [MPa] Young modulus along x
E_y1 = E_x1;                E_y2 = E_x2;    % [MPa] Young modulus along y
G_xy1 = E_x1/(2*(1+0.3));   % [MPa] Shear modulus 
G_xy2 = E_x2/(2*(1+0.3));   % [MPa] Shear modulus 
% (metallic materials E_x1/(2*(1+0.3)) -> you might need to change it

% Core
G_cx = 86.9;          % [MPa] Shear modulus core in x
G_cy = 86.9;          % [MPa] Shear modulus core in y
E_cx = 1.13e-2;       % [MPa] mostly used to check for weak core conditon
E_cy = 1.13e-2;       % [MPa] mostly used to check for weak core conditon
E_c  = 1.13e-2;       % [MPa] mostly used to check for weak core conditon

% % % % % % % % % % COORDINATES FOR CALCULATIONS % % % % % % % % % % % % % 

% Forces T            % Moments M           % Deflections w 
x_value_Tx = 0;       x_value_Mx = a/2;     x_value_w = a/2;
y_value_Tx = b/2;     y_value_Mx = b/2;     y_value_w = b/2;
x_value_Ty = a/2;     x_value_My = a/2;
y_value_Ty = 0;       y_value_My = b/2;

%% Calculating rigidities

% Stiffnesses
D_x = (E_x1*t_1*E_x2*t_2*d^2) / (E_x1*t_1 + E_x2*t_2);
D_y = D_x;
D   = sqrt(D_x*D_y);

S_x = (G_cx*d^2) / t_c;   % Shear stiffness
S_y = S_x;               
S   = sqrt(S_x*S_y);

% Poisson ratios of the whole plate, approximations!
% v_xy = 0.25*sqrt(D_x/D_y);
% v_yx = v_xy*(D_x/D_y);

% Poisson for isotropic
v = 0.3;            
v_xy = v;
v_yx = v;

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

%% Computing deflections w_b and w_s

syms m n
syms x y

% Equation for q_mn, common for both
q_mn = 4/(a*b) * int(int(q*sin((m*pi*x)/a)*sin((n*pi*y)/b),y,0,b),x,0,a);

% Create function to iterate
w_b_equ = (q_mn * (1-v^2) * sin(m*pi*x/a) * sin(n*pi*y/b)) / (D * ((m*pi/a)^2 + (n*pi/b)^2)^2);
w_b_loop = vpa(subs(w_b_equ, [x y],[x_value_w y_value_w]));
w_b_loop = matlabFunction(w_b_loop);
w_s_equ = (q_mn * sin(m*pi*x/a) * sin(n*pi*y/b)) / (S * ((m*pi/a)^2 + (n*pi/b)^2));
w_s_loop = vpa(subs(w_s_equ, [x y],[x_value_w y_value_w]));
w_s_loop = matlabFunction(w_s_loop);

% Loop for m and n
w_b_current = 0; w_s_current = 0;  
for j=1:1:50
for i=1:1:50
    m = i; 
    n = j;
    w_b_current = w_b_current + w_b_loop(m,n);
    w_s_current = w_s_current + w_s_loop(m,n);
end
end
w_b = w_b_current
w_s = w_s_current

% Overall deflection w
w = w_b + w_s

%% Computation for stresses, M_x, M_y, T_x, T_y - Approximate solution!

syms m n
syms x y

% Create function to iterate - using w_b and w_s
M_x_equ = - D_x/(1-v_xy*v_yx) * (diff(w_b_equ,x,2)+v_yx*diff(w_b_equ,y,2));
M_x_equ = vpa(subs(M_x_equ, [x y], [x_value_Mx y_value_Mx]));
M_x_equ = matlabFunction(M_x_equ);
M_y_equ = - D_y/(1-v_xy*v_yx) * (diff(w_b_equ,y,2)+v_xy*diff(w_b_equ,x,2));
M_y_equ = vpa(subs(M_y_equ, [x y], [x_value_My y_value_My]));
M_y_equ = matlabFunction(M_y_equ);

% Loop for m and n
M_x_current = 0; M_y_current = 0;  
for j=1:1:50
for i=1:1:50
    m = i;
    n = j;
    M_x_current = M_x_current + M_x_equ(m,n);
    M_y_current = M_y_current + M_y_equ(m,n);
end
end
M_x = M_x_current
M_y = M_y_current

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

syms m n
syms x y

% Create function to iterate
T_x_equ = S_x * diff(w_s_equ,x);
T_x_loop = vpa(subs(T_x_equ, [x y],[x_value_Tx y_value_Tx]));
T_x_loop = matlabFunction(T_x_loop);
T_y_equ = S_y * diff(w_s_equ,y);
T_y_loop = vpa(subs(T_y_equ, [x y],[x_value_Ty y_value_Ty]));
T_y_loop = matlabFunction(T_y_loop);

% Loop for m and n
T_x_current = 0; T_y_current = 0;  
for j=1:1:50
for i=1:1:50
    m = i;
    n = j;
    T_x_current = T_x_current + T_x_loop(m,n);
    T_y_current = T_y_current + T_y_loop(m,n);
end
end
T_x = T_x_current
T_y = T_y_current



%% Computation for w directly to check the result

% syms m n
% syms x y
% 
% % Preparing equations for w
% w_equ = (1 + (D / (S*(1-v^2))) * ((m*pi/a)^2 + (n*pi/b)^2)) / ...
%       (m*n*((m*pi/a)^2 + (n*pi/b)^2)^2) * sin(m*pi*x/a) * sin(n*pi*y/b);
% w_equ = vpa(subs(w_equ, [x y], [x_value_w y_value_w]));
% w_equ = matlabFunction(w_equ);
% 
% % Loop to substitute m and n
% w = 0;
% for i=1:1:50
% for j=1:1:50
%     m = i;
%     n = j;
%     w = w + w_equ(m,n);
% end
% end
% w = (16*q*(1-v^2))/(pi^2*D) * w



%% Stress calculation

% 
% 
% sigma_fx1 = (N_x*E_x1)/(E_x1*t_1+E_cx*t_c+E_x2*t_2) - (M_x*E_x1*E_x2*t_2*d)/(D_x*(E_x1*t_1+E_x2*t_2))
% sigma_fy1 = (N_y*E_y1)/(E_y1*t_1+E_cy*t_c+E_y2*t_2) - (M_y*E_y1*E_y2*t_2*d)/(D_y*(E_y1*t_1+E_y2*t_2))
% sigma_fx2 = (N_x*E_x2)/(E_x1*t_1+E_cx*t_c+E_x2*t_2) + (M_x*E_x1*E_x2*t_1*d)/(D_x*(E_x1*t_1+E_x2*t_2))
% sigma_fy2 = (N_y*E_y2)/(E_y1*t_1+E_cy*t_c+E_y2*t_2) + (M_y*E_y1*E_y2*t_1*d)/(D_y*(E_y1*t_1+E_y2*t_2))
% 
% tau_fxy = M_xy/(sqrt(t_1*t_2)*d) % sligtly modified to have one t_f
% tau_cxz = T_x/d
% tau_cyz = T_y/d


