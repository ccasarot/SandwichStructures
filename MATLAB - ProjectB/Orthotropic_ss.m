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
%%%%%%%    Final file for the orthotropic simply supported case     %%%%%%%
%%%%%%%    In the first part needs to declare material and type     %%%%%%%
%%%%%%%    of core. Failure criterias were added as well, with      %%%%%%%
%%%%%%%    the exception of the ones depending on the previous      %%%%%%%
%%%%%%%    course of compistes, which were not working. No          %%%%%%%
%%%%%%%    solution was found in time.                              %%%%%%%
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

% Collection of dimensions of plate panels from the hellicopter floor:
%                 a1   a2   b1   b2  b3  b4   c1   c2   c3
% a_dimension = [2060 2060 1630 670 980 1630 1030 1030 2320];
% b_dimension = [650  650  825  650 670 825  870  870  780];

%% Input data

% % % % % % % % % % % % PLATE PARAMETERS % % % % % % % % % % % % % % % % % 

a = 2320;          % [mm] Length for edge a
b = 780;           % [mm] Length for edge b
t_1 = 2;           % [mm] Thickness for face sheet 1
t_2 = t_1;         % [mm] Thickness for face sheet 2
t_c = 12;          % [mm] Thickness of the core
d   = t_c+t_1;     % [mm] distance between axis of face sheets

% % % % % % % % % % % % % % % % LOAD % % % % % % % % % % % % % % % % % % % 

q_dist = 35.3/1000; % 800 * 9.81 * 4.5/(10^6);      % [MPa] Uniformely distributed load
P_local = 100;              % [N] Local load
A_local = 100;              % [mm^2] Area where the local load acts
q_conc = P_local/A_local;   % [MPa] concentrated load

% % % % % % % % % % % % MATERIAL PARAMETERS % % % % % % % % % % % % % % % % 

% Upper sheet (1)        % Lower sheet (2)
E_x1  = 71000;           E_x2  = 71000;       % [MPa] Young modulus along x
E_y1  = E_x1;            E_y2  = E_x2;        % [MPa] Young modulus along y
G_xy1 = 26692;           G_xy2 = 26692;       % [MPa] Shear modulus 
% (metallic materials E_x1/(2*(1+0.3)) -> you might need to change it

v = 0.33; % If metals (if non metallic face sheets, approx. are later given)   
% it is the same in the other direction

% Core
G_cx = 155;        % [MPa] Shear modulus core in x
G_cy = 320.5;      % [MPa] Shear modulus core in y
E_cx = 3.05e-1;    % [MPa] mostly used to check for weak core conditon
E_cy = 7.395e-1;   % [MPa] mostly used to check for weak core conditon
E_c  = 552;        % [MPa] (z!) mostly used to check for weak core conditon



% % % % % % % % % % COORDINATES FOR CALCULATIONS % % % % % % % % % % % % % 

% Forces T            % Moments M           % Deflections w 
x_value_Tx = 0;       x_value_Mx = a/2;     x_value_w = a/2;
y_value_Tx = b/2;     y_value_Mx = b/2;     y_value_w = b/2;
x_value_Ty = a/2;     x_value_My = a/2;
y_value_Ty = 0;       y_value_My = b/2;
                      x_value_Mxy = 0;
                      y_value_Mxy = 0;

% Compute the integration extremes - q_mn conc
ext_1_min = a/2 - sqrt(A_local)/2; % The are where the load is applied is 100x100 mm
ext_1_max = a/2 + sqrt(A_local)/2; % and centered in the plate
ext_2_min = b/2 - sqrt(A_local)/2;
ext_2_max = b/2 + sqrt(A_local)/2;

%% TYPE OF MATERIAL - CONTROL PANEL

% MAXIMUM STRESS FOR MATERIAL
sigma_y = 270;       % [MPa] Face sheets yielding
tau_max = 1.5;       % [MPa] Max shear for the core
sigma_max_compression_core = 2.31;
sf_composite = 1.5*1.2;            % [\] Safety factor for face sheets
sf = 1.5;

sigma_max = sigma_y / sf;   % [MPa] face sheets, acount for safety factor ec
tau_max = tau_max / sf;     % core
sigma_max_compression_core = sigma_max_compression_core / sf;

% FACE SHEETS
Metallic_Face_Sheets = 'y';

% IF HONEYCOMB
Using_Honeycomb = 'y';
SquareCells = 'n';
HexagonalCells = 'y';
s = 5; % [mm] Ratio of the circle "inside the hexagon" - If hexagonal cells



%% Calculating rigidities (From previous course)

if Metallic_Face_Sheets == 'n'
    E1 = 121000;   % Mpa        % these are the data from granta
    E2 = 8600;      % Mpa
    G12 = 47000;    % Mpa
    V12 = 0.27;
    
    S_2D=[1/E1     -V12/E1  0
         -V12/E1   1/E2     0
          0         0       1/G12] ;
    
    Qloc=inv(S_2D) ;
    
    % % Layup 
    layup = [0, 45, -45, 90] ;  % layout we used. And ptoportional ply thickness     
    h1 = t_1/4;
    h2 = t_1/4;
    h3 = t_1/4;
    h4 = t_1/4;
    
    t_tot = t_1;
    
    nply = size(layup, 2) ;
    Z = [-h2/2-h3, -h2/2, h2/2, h2/2+h3, h2/2+h3+h4];     % in mm
    [T,Q] = build_TQ(layup, Qloc) ;
    [A,B,D] = build_ABD(Q, Z) ;
    ABD = [A,B;B,D] ;
    
    v_xy =  A(1, 2)/A(2, 2);
    v_yx =  A(1, 2)/A(1, 1);
    
    E_x1 = (A(1, 1) * (1-v_xy*v_yx)) / t_tot;
    E_y1 = (A(2, 2) * (1-v_xy*v_yx)) / t_tot;
    G_xy1 = A(3, 3) / t_tot;

    E_x2 = E_x1;
    E_y2 = E_y1;
    G_xy2 = G_xy1;
    E_x2
    E_y2

end



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

% % % % % % % % % % % % % % SUPERIMPOSITION % % % % % % % % % % % % % % % % 
% Superimposition of parameters, if they are already given

% % IF YOU ARE USING THE CODE FOR ISOTROPIC D_xy is different !
% D=sqrt(D_x*D_y);
% D_xy=D/(1+v)         % not sure if we need this, probably not



% % % % % % % % % % % % % % PROGRAM STARTS % % % % % % % % % % % % % % % % 



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

%% Computing transverse forces

syms m n
syms x y

% Equation for transverse forces
T_x_equ = (X_mn*q_mn)/Z_mn * cos(m*pi*x/a) * sin(n*pi*y/b);
T_x_loop = vpa(subs(T_x_equ, [x y],[x_value_Tx y_value_Tx]));
T_x_loop = matlabFunction(T_x_loop);

T_y_equ = -(Y_mn*q_mn)/Z_mn * sin(m*pi*x/a) * cos(n*pi*y/b);
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

%% Computing deflection w

syms m n
syms x y

% Equation for deflections
w_equ = - ((W_mn*q_mn)/Z_mn) * sin(m*pi*x/a) * sin(n*pi*y/b);
w_loop = vpa(subs(w_equ, [x y],[x_value_w y_value_w]));
w_loop = matlabFunction(w_loop);

w_current = 0;
for j=1:1:50
for i=1:1:50
    m = i; 
    n = j;
    w_current = w_current + w_loop(m,n);
end
end
w = w_current

%% Computing moments

syms m n
syms x y

% Equation for moments
M_x_loop = (-D_x/(1 - v_xy*v_yx)) * (W_mn * ((m*pi/a)^2 + v_yx*(n*pi/b)^2) ...
    + (m*pi/a)*(X_mn/S_x) - v_yx*(n*pi/b)*(Y_mn/S_y)) * (q_mn/Z_mn) * ...
    sin(m*pi*x/a) * sin(n*pi*y/b);
M_x_loop = vpa(subs(M_x_loop, [x y],[x_value_Mx y_value_Mx]));
M_x_loop = matlabFunction(M_x_loop);

M_y_loop = (-D_y/(1 - v_xy*v_yx)) * (W_mn * ((n*pi/b)^2 + v_xy*(m*pi/a)^2) ...
    - (n*pi/b)*(Y_mn/S_y) + v_xy*(m*pi/a)*(X_mn/S_x)) * (q_mn/Z_mn) * ...
    sin(m*pi*x/a) * sin(n*pi*y/b);
M_y_loop = vpa(subs(M_y_loop, [x y],[x_value_My y_value_My]));
M_y_loop = matlabFunction(M_y_loop);

% Formulation from simple derivation, not from the ones with terms ecc...
M_xy_loop = - D_xy/2 * (diff((diff(w_equ,y)-T_y_equ/S_y),x) + diff((diff(w_equ,x)-T_x_equ/S_x),y));
M_xy_loop = vpa(subs(M_xy_loop, [x y],[x_value_Mxy y_value_Mxy]));
M_xy_loop = matlabFunction(M_xy_loop);

% Loop for m and n
M_x_current = 0; M_y_current = 0; M_xy_current = 0;
for j=1:1:50
for i=1:1:50
    m = i; 
    n = j;
    M_x_current = M_x_current + M_x_loop(m,n);
    M_y_current = M_y_current + M_y_loop(m,n);
    M_xy_current = M_xy_current + M_xy_loop(m,n);
end
end
M_x = M_x_current
M_y = M_y_current
M_xy = M_xy_current



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%% Stress calculation
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

N_x=0;
N_y=0;
             % First term: 0 if N stresses = 0
sigma_fx1 = (N_x*E_x1)/(E_x1*t_1+E_cx*t_c+E_x2*t_2) - (M_x*E_x1*E_x2*t_2*d)/(D_x*(E_x1*t_1+E_x2*t_2))
sigma_fy1 = (N_y*E_y1)/(E_y1*t_1+E_cy*t_c+E_y2*t_2) - (M_y*E_y1*E_y2*t_2*d)/(D_y*(E_y1*t_1+E_y2*t_2))
sigma_fx2 = (N_x*E_x2)/(E_x1*t_1+E_cx*t_c+E_x2*t_2) + (M_x*E_x1*E_x2*t_1*d)/(D_x*(E_x1*t_1+E_x2*t_2))
sigma_fy2 = (N_y*E_y2)/(E_y1*t_1+E_cy*t_c+E_y2*t_2) + (M_y*E_y1*E_y2*t_1*d)/(D_y*(E_y1*t_1+E_y2*t_2))

tau_fxy1 = M_xy/(t_1*d) 
tau_fxy2 = M_xy/(t_2*d) 

tau_cxz = T_x/d
tau_cyz = T_y/d



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%% Failure criterias
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% 1 Face and core fracture 
% 2 Wrinkling (local face buckling) 
% 3 Face dimpling 
% 4 Global buckling – Beams 
% 5 Shear crimping  – Beams 
% 6 Global buckling 
% 7 Core indentation 

% % % % % % % % % % Wrinkling (local face buckling) % % % % % % % % % % % %
disp('Wrinkling')

% "In practice formulas" check for both faces and both directions
sigma_cr1x = 0.5*(E_x1*E_cx*G_cx)^(1/3)
sigma_cr1y = 0.5*(E_y1*E_cy*G_cy)^(1/3)
sigma_cr2x = 0.5*(E_x2*E_cx*G_cx)^(1/3)
sigma_cr2y = 0.5*(E_y2*E_cy*G_cy)^(1/3)
if sigma_fx1 > sigma_cr1x || sigma_fy1 > sigma_cr1y || sigma_fx2 > sigma_cr2x || sigma_fy2 > sigma_cr2y 
    disp('Wrinkling - NOT VERIFIED')
end

% % % % % % % % % % % % % Face dimpling % % % % % % % % % % % % % % % % % % 
disp('Face dimpling')

if Using_Honeycomb=='y'
    if SquareCells=='y'
    % Square honeycombs
    sigma_cr1 = 2.5 * sqrt(E_x1*E_y1) * (t_1/a)^2 
    sigma_cr2 = 2.5 * sqrt(E_x2*E_y2) * (t_2/a)^2 
        if sigma_fx1 > sigma_cr1 || sigma_fy1 > sigma_cr1 || sigma_fx2 > sigma_cr2 || sigma_fy2 > sigma_cr2
            disp('Dimpling - NOT VERIFIED')
        end
    end
    if HexagonalCells=='y'
    % Hexagonal honeycombs
    sigma_cr1 = (2*sqrt(E_x1*E_y1)/(1-(sqrt(v_xy*v_yx))^2))*(t_1/s)^2
    sigma_cr2 = (2*sqrt(E_x2*E_y2)/(1-(sqrt(v_xy*v_yx))^2))*(t_2/s)^2
        if sigma_fx1 > sigma_cr1 || sigma_fy1 > sigma_cr1 || sigma_fx2 > sigma_cr2 || sigma_fy2 > sigma_cr2
            disp('Dimpling - NOT VERIFIED')
        end
    end
end

% % % % % % % % % % % % % Shear crimping % % % % % % % % % % % % % % % % % 

% This should be for beams 

% % % % % % % % % % % % % Global buckling % % % % % % % % % % % % % % % % % 
% Remember that P loads are in force * width unit
% disp('Global buckling')
% 
% % Define shear coefficient
% theta_x = D_x/(a^2*S_x*(1-v_xy*v_yx));
% % Mind minimum K trough minimsation process
% K_vector = zeros(15,1);
%     for i=1:length(K_vector)
%         m=i;
%         K_vector(i) = ((1 + S_x/S_y * (a/(m*b))^2) * ((m*b/a)^2 + D_y/D_x * (a/(m*b))^2 + 2 * (v_yx + D_xy * (1 - v_xy*v_yx)/D_x))) / ((1 + S_x/S_y * (a/(m*b))^2) + pi^2 * theta_x * (a/b)^2 * ((m*b/a)^2 + D_y/D_x * (a/(m*b))^2 + 2 * (v_yx + D_xy * (1 - v_xy*v_yx)/D_x)));
%     end
% K = min(K_vector);
% % Find critical load
% P_cx = (pi^2*D_x*K)/(b^2*(1-v_xy*v_yx))
%     if P_mio_x > P_cx
%         disp('Global buckling, load on x direction - NOT VERIFIED')
%     end
% 
% % Define shear coefficient
% theta_y = D_y/(b^2*S_y*(1-v_xy*v_yx));
% % Mind minimum K trough minimsation process
% K_vector = zeros(15,1);
%     for i=1:length(K_vector)
%         m=i;
%         K_vector(i) = ((1 + S_y/S_x * (b/(m*a))^2) * ((m*a/b)^2 + D_x/D_y * (b/(m*a))^2 + 2 * (v_yx + D_xy * (1 - v_xy*v_yx)/D_y))) / ((1 + S_y/S_x * (b/(m*a))^2) + pi^2 * theta_y * (b/a)^2 * ((m*a/b)^2 + D_x/D_y * (b/(m*a))^2 + 2 * (v_yx + D_xy * (1 - v_xy*v_yx)/D_y)));
%     end
% K = min(K_vector);
% % Find critical load
% P_cy = (pi^2*D_y*K)/(a^2*(1-v_xy*v_yx))
%     if P_mio_y > P_cy
%         disp('Global buckling, load on y direction - NOT VERIFIED')
%     end

% This K formulation is for two sides clamped and 2 SS. Load on the SS side (see slide)
% K = ((1 + 4*S_y/(3*S_x)*(a/(m*b))^2) * ((m*b/a)^2 + 16*D_y/(3*D_x)*(a/(m*b))^2 + 8/3*(v_yx + D_xy*(1 - v_xy*v_yx)/D_x))) / ((1 + 4*S_y/(3*S_x)*(a/(m*b))^2) + pi^2*theta_x*(a/b)^2 * ((m*b/a)^2 + 16*D_y/(3*D_x)*(a/(m*b))^2 + 8/3*(v_yx + D_xy*(1 - v_xy*v_yx)/D_x)));

% % % % % % % % % % % % % Core indentation % % % % % % % % % % % % % % % % 
disp('Core indentation')

% Basically a simple verification for compression of the core
sigma_compression_core = P_local/A_local
    if sigma_compression_core > sigma_max_compression_core
        disp('Core indentation - NOT VERIFIED')
    end



% % % % % % % % % % % Face and core fracture % % % % % % % % % % % % % % %



% % From the other course. something is wrong here, but we never had the
% time to check what it was.



disp('Face and core fracture')

if Metallic_Face_Sheets == 'y'
    % Then you use Mohr circles to find the stress
    sigma_cr1 = (sigma_fx1+sigma_fy1)/2 + sqrt(((sigma_fx1-sigma_fy1)/2)^2 + tau_fxy1^2)
    sigma_cr2 = (sigma_fx2+sigma_fy2)/2 + sqrt(((sigma_fx2-sigma_fy2)/2)^2 + tau_fxy2^2)
    if abs(sigma_cr1) > sigma_max || abs(sigma_cr2) > sigma_max
        disp('Metalic face fracture - NOT VERIFIED')
    end
end

% Simple check that the core does not collapse due to shear
if tau_cxz > tau_max || tau_cyz > tau_max
    disp('Core fracture - NOT VERIFIED')
end

if Metallic_Face_Sheets == 'n'
    
    % Strength parameters for failure criterions calculations
    
    sig_1t = 2231*sf_composite;   % Mpa      old values:  1100   
    sig_1c = 1082*sf_composite;   % Mpa                   600
    sig_2t = 29*sf_composite;     % Mpa                   20
    sig_2c = 100*sf_composite;    % Mpa                   140
    
    eps_1t = sig_1t / E1 ;   % Mpa
    eps_1c = sig_1c / E1 ;   % Mpa
    eps_2t = sig_2t / E2 ;   % Mpa
    eps_2c = sig_2c / E2 ;   % Mpa
    
    tau_12_max = 60*sf_composite;        % Mpa
    gam_12_max = tau_12_max/G12 ;
    
    % I have to do an integration considering q distributed load
    N = [0; 0; 0] ;    % N
    M = [M_x; M_y; M_xy] ;   % N/mm
    incr_load_vect = [N;M]   ;
    
    strain_vect = ABD\incr_load_vect ; % incre just because its iterative.
    
    % List of stress vectors for each ply (top and bottom stresses)
    eps_top_loc = {} ;
    eps_bot_loc = {} ;
    
    
    sigma_top_loc = {} ;
    sigma_bot_loc = {} ;
    
    
    for i=1:nply
        
        % Calculate local stresses and strains 
        
        eps0 = strain_vect(1:3);
        K0 = strain_vect(4:6);
        
        eps_top = eps0 + Z(i) * K0 ;
        eps_bot = eps0 + Z(i+1) * K0 ; 
        eps_mid = eps0 + (Z(i+1)+Z(i))/2*K0 ;
        
        eps_top_loc{i} = T{i}'*eps_top ; 
        eps_bot_loc{i} = T{i}'*eps_bot ;
        eps_mid_loc{i} = T{i}'*eps_mid ;
        
        sigma_top_loc{i} = Qloc*eps_top_loc{i} ;
        sigma_bot_loc{i} = Qloc*eps_bot_loc{i} ;
        sigma_mid_loc{i} = Qloc*eps_mid_loc{i} ;
    
        % for each fibre of the ply 
        for fibre=1:2 
            if fibre==1 %top
                sig_1 = sigma_top_loc{i}(1) ;
                sig_2 = sigma_top_loc{i}(2) ;
                tau_12 = sigma_top_loc{i}(3) ;
                eps_1 = eps_top_loc{i}(1) ;
                eps_2 = eps_top_loc{i}(2) ;
                gam_12 = eps_top_loc{i}(3) ;
            elseif fibre==2 %bottom
                sig_1 = sigma_bot_loc{i}(1) ;
                sig_2 = sigma_bot_loc{i}(2) ;
                tau_12 = sigma_bot_loc{i}(3) ;
                eps_1 = eps_bot_loc{i}(1) ;
                eps_2 = eps_bot_loc{i}(2) ;
                gam_12 = eps_bot_loc{i}(3) ;
    % %         elseif fibre==3 %middle
    % %             sig_1 = sigma_mid_loc{i}(1) ;
    % %             sig_2 = sigma_mid_loc{i}(2) ;
    % %             tau_12 = sigma_mid_loc{i}(3) ;
    % %             eps_1 = eps_mid_loc{i}(1) ;
    % %             eps_2 = eps_mid_loc{i}(2) ;
    % %             gam_12 = eps_mid_loc{i}(3) ;
            end
        end
            % Check failure criterion
            i;
            HOF = hoffman(sig_1, sig_2, sig_1t, sig_1c, sig_2t, sig_2c, tau_12, tau_12_max) ; 
            if HOF == 1
               disp('Ply number'); disp(i);
               disp('Face fracture - NOT VERIFIED')
            end
        end
    end

% functions we need


%% Build Q, T, A B D matrixes
function [A,B,D] = build_ABD(Q, Z)
nply = size(Z,2) - 1  ;
A = zeros(); 
B = zeros(); 
D = zeros(); 

    for i=1:nply 
    
        A = A + Q{i} * (Z(i+1) - Z(i)) ;
        B = B + 1/2 * Q{i} * (Z(i+1)^2 - Z(i)^2); 
        D = D + 1/3 * Q{i} * (Z(i+1)^3 - Z(i)^3); 
        
    end 
end 


function [T,Q] = build_TQ(layup, Qloc)

nply = size(layup,2) ;
Q = {} ;
T = {} ;

for i=1:nply 
    teta = layup(i) ;
    c= cos(teta*pi/180) ;
    s= sin(teta*pi/180) ;
    
    T{i} = [c^2 s^2 -2*s*c
        s^2  c^2  2*s*c  
        s*c -s*c  c^2-s^2] ;
    
    
    Q{i} = T{i}*Qloc*T{i}' ;
    
end 
end 

%% Failure criterions function

function [fail] = hoffman(sig_1, sig_2, sig_1t, sig_1c, sig_2t, sig_2c, tau_12, tau_12_max)

H = (sig_1)^2 / (sig_1c*sig_1t) - (sig_1 * sig_2)/(sig_1c * sig_1t) + (sig_2^2/(sig_2c * sig_2t)) - ((sig_1t - sig_1c)/(sig_1c*sig_1t)*sig_1) - (((sig_2t - sig_2c)/(sig_2c*sig_2t)*sig_2)) + (tau_12 / tau_12_max)^2 ;
if H < 1
    fail = 0; 
else
    fail = 1; 
end
end 


