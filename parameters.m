%parameters file for mean-field model Main.m file
%loads any fixed parameters that we have in the model; other parameters are
%inputs in the function Main itself

%Erin Angelini, 5.25.18

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters specific to Main.m
r_N = 4; %spindle radius, same value used in stochastic model

b_basal = 15; %wild type aspect ratio is a:b = 25:15, so the basal value of 
%b is set at 15 when we elongate aspect ratio 
a_basal = 15; %this is the basal short axis is if we instead elongate along 
%the y-axis and keep the x-axis fixed

%for volume scaling (by factor of 1.6), AR = 1.67 (wild type)
a_WTscale = 40; b_WTscale = 24;
N_d_basal = 128; %this value is for when we want to fix the number of 
%cortical dyneins instead of scaling them with cell dimensions

arc = 8; %Let-99 push band arc length

d_alpha = pi/200; %step size for alpha values W
end_alpha = pi; %we want to go from alpha = 0 to alpha = pi
A = 0:d_alpha:end_alpha; %range of alpha values with the above step size

%parameters specific to work.m
start = 0; %starting alpah value for approximating integral
d_alpha_work = 0.005; %step size for Riemann sum

%parameters specific to torque.m
phi_pi = pi; %this is just for the case that aster spread is pi
rho_asymm = 1157/2.324; %MT angular density for asymmetric envs, value 
%from computational model; #mts/angle
rho_symm = 1000/(2*pi/3); %MT angular density for symmetric envs, value 
%from computational model; #mts/angle

%parameters specific to find_prob.m
myCutoff = 0.2*a; %this is the horizontal 60:40 mark of the cell
P_p = 1; %posterior binding probability, from computational model
P_a = 0.65; %anterior binding probability, from computational model

%parameters specific to find_pull.m
beta = 1; %this is the scaling factor of pulling forces
