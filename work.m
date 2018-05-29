function W = work(alpha, a, b, set_psi, r_N, ctr, phi_1, phi_2, tvec, AP, LD, elas)
%calculates work done at given angle alpha by approximating the integral of the
%torque from 0 to alpha
%inputs are spindle angle alpha, major/minor axes a and b of the ellipse,
%set_psi dynein locations

%for varying center and spindle envelopes: also takes spindle radius r_N,
%center ctr, and upper/lower envelope spreads phi_1 and phi_2

%AP, LD, and elas are strings that control presence of anterior to posterior 
%difference in MT binding probabilities (AP), presence of MT-length 
%dependent forces, and presence of MT buckling (respectively)


%Erin Angelini, 5.25.18

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load parameters specific to work.m
parameters

if alpha < start %if we are starting integral from an angle > alpha
    all_alphas = alpha:d_alpha_work:start-d_alpha_work;
else
    all_alphas = start:d_alpha_work:alpha-d_alpha_work;
end

%Riemann summ approx. of work integral
W = 0;
for i = 1:length(all_alphas)
    alpha_prime = all_alphas(i);
    T = torque(alpha_prime, a, b, set_psi, r_N, ctr, phi_1, phi_2, tvec, AP, LD, elas);
    W = W + T*d_alpha; 
end
if alpha >= start
    W =-W;%reverse sign
end

end
