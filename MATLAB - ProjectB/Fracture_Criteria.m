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
%%%%%%%    This file collects the fracture criteria before          %%%%%%%
%%%%%%%    implementing them in the orthotropic file.               %%%%%%%
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%% Stress calculation
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

P_local = 100   % [N] Local load
A_local = 100   % [mm^2] Area where the local load acts

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

%% TYPE OF MATERIAL - CONTROL PANEL

% IF HONEYCOMB
Using_Honeycomb = 'y';
SquareCells = 'n';
HexagonalCells = 'y';
s = 5; % [mm] Ratio of the circle "inside the hexagon" - If hexagonal cells

% FACE SHEETS
Metallic_Face_Sheets = 'y';

% MAXIMUM STRESS FOR MATERIAL
sigma_y = 240
sf = 1.5 

sigma_max = sigma_y/sf  % [Mpa] face sheets, acount for safety factor ec
tau_max = 50  % core
P_mio_x = 50 % (N loads !!!) load on the long side normal to the surface... load per width unit?
P_mio_y = 200
sigma_max_compression_core = 200

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

% % % % % % % % % % % Face and core fracture % % % % % % % % % % % % % % %
disp('Face and core fracture')

if Metallic_Face_Sheets == 'y'
    % Then you use Mohr circles to find the stress
    sigma_cr1 = (sigma_fx1+sigma_fy1)/2 + sqrt(((sigma_fx1-sigma_fy1)/2)^2 + tau_fxy1^2)
    sigma_cr2 = (sigma_fx2+sigma_fy2)/2 + sqrt(((sigma_fx2-sigma_fy2)/2)^2 + tau_fxy2^2)
    if abs(sigma_cr1) > sigma_max || abs(sigma_cr2) > sigma_max
        disp('Face fracture - NOT VERIFIED')
    end
end

if Metallic_Face_Sheets == 'n'
    % Then you use Mohr circles to find the stress
    disp('Face and core fracture - NOT IMPLEMENTED YET, you should use CLT')
end

% Simple check that the core does not collapse due to shear
if tau_cxz > tau_max || tau_cyz > tau_max
    disp('Core fracture - NOT VERIFIED')
end

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

% This should be for beams as far as I understood

% % % % % % % % % % % % % Global buckling % % % % % % % % % % % % % % % % % 
% Remember that P loads are in force * width unit
disp('Global buckling')

% Define shear coefficient
theta_x = D_x/(a^2*S_x*(1-v_xy*v_yx));
% Mind minimum K trough minimsation process
K_vector = zeros(15,1);
    for i=1:length(K_vector)
        m=i;
        K_vector(i) = ((1 + S_x/S_y * (a/(m*b))^2) * ((m*b/a)^2 + D_y/D_x * (a/(m*b))^2 + 2 * (v_yx + D_xy * (1 - v_xy*v_yx)/D_x))) / ((1 + S_x/S_y * (a/(m*b))^2) + pi^2 * theta_x * (a/b)^2 * ((m*b/a)^2 + D_y/D_x * (a/(m*b))^2 + 2 * (v_yx + D_xy * (1 - v_xy*v_yx)/D_x)));
    end
K = min(K_vector);
% Find critical load
P_cx = (pi^2*D_x*K)/(b^2*(1-v_xy*v_yx))
    if P_mio_x > P_cx
        disp('Global buckling, load on x direction - NOT VERIFIED')
    end

% Define shear coefficient
theta_y = D_y/(b^2*S_y*(1-v_xy*v_yx));
% Mind minimum K trough minimsation process
K_vector = zeros(15,1);
    for i=1:length(K_vector)
        m=i;
        K_vector(i) = ((1 + S_y/S_x * (b/(m*a))^2) * ((m*a/b)^2 + D_x/D_y * (b/(m*a))^2 + 2 * (v_yx + D_xy * (1 - v_xy*v_yx)/D_y))) / ((1 + S_y/S_x * (b/(m*a))^2) + pi^2 * theta_y * (b/a)^2 * ((m*a/b)^2 + D_x/D_y * (b/(m*a))^2 + 2 * (v_yx + D_xy * (1 - v_xy*v_yx)/D_y)));
    end
K = min(K_vector);
% Find critical load
P_cy = (pi^2*D_y*K)/(a^2*(1-v_xy*v_yx))
    if P_mio_y > P_cy
        disp('Global buckling, load on y direction - NOT VERIFIED')
    end

% This K formulation is for two sides clamped and 2 SS. Load on the SS side (see slide)
% K = ((1 + 4*S_y/(3*S_x)*(a/(m*b))^2) * ((m*b/a)^2 + 16*D_y/(3*D_x)*(a/(m*b))^2 + 8/3*(v_yx + D_xy*(1 - v_xy*v_yx)/D_x))) / ((1 + 4*S_y/(3*S_x)*(a/(m*b))^2) + pi^2*theta_x*(a/b)^2 * ((m*b/a)^2 + 16*D_y/(3*D_x)*(a/(m*b))^2 + 8/3*(v_yx + D_xy*(1 - v_xy*v_yx)/D_x)));

% % % % % % % % % % % % % Core indentation % % % % % % % % % % % % % % % % 
disp('Core indentation')

% Basically a simple verification for compression of the core
sigma_compression_core = P_local/A_local
    if sigma_compression_core > sigma_max_compression_core
        disp('Core indentation - NOT VERIFIED')
    end


    