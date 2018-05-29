function alpha_halfplane(a, r_N, ctr, phi_1, phi_2)
%for calculating the angle at which envelope endpoints go above or below
%x-axis
%runs each time you run Main so you can adjust a, r_N, ctr, phi_1, or phi_2

%Erin Angelini, 6.30.17

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

syms z
k1 = sqrt((r_N*cos(z) + ctr - a)^2 + r_N^2*sin(z)^2);
k2 = sqrt((r_N*cos(z) + ctr + a)^2 + r_N^2*sin(z)^2);

Alphaq1 = vpasolve(z == asin((k1/(a-ctr))*sin(phi_1/2)), z);
alphaq1 = double(Alphaq1); 
Alphaq2 = vpasolve(z == pi + asin((k2/(a+ctr))*sin(phi_1/2)), z);
alphaq2 = double(Alphaq2);
%if alpha >= alphaq1 and alpha < alphaq2, then q above

Alphap1 = vpasolve(z == 2*pi - asin((k1/(a-ctr))*sin(phi_1/2)), z);
alphap1 = double(Alphap1); 
Alphap2 = vpasolve(z == pi - asin((k2/(a+ctr))*sin(phi_1/2)), z);
alphap2 = double(Alphap2);
%if alpha >= alphap1 or alpha < alphap2, then p above

save('halfplane1.mat', 'alphaq1', 'alphap1', 'alphaq2', 'alphap2')

AlphaQ1 = vpasolve(z == asin((k1/(a-ctr))*sin(phi_2/2)), z);
alphaQ1 = double(AlphaQ1); 
AlphaQ2 = vpasolve(z == pi + asin((k2/(a+ctr))*sin(phi_2/2)), z);
alphaQ2 = double(AlphaQ2);


AlphaP1 = vpasolve(z == 2*pi - asin((k1/(a-ctr))*sin(phi_2/2)), z);
alphaP1 = double(AlphaP1); 
AlphaP2 = vpasolve(z == pi - asin((k2/(a+ctr))*sin(phi_2/2)), z);
alphaP2 = double(AlphaP2);

save('halfplane2.mat', 'alphaQ1', 'alphaP1', 'alphaQ2', 'alphaP2')
end