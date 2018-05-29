function tvec = let99(a, b, arc)
%finds the positions of the Let-99 bands of arclength arc localized at the 
%60:40 mark along the cell cortex

%Erin Angelini, 7.6.17

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%find x-coordinate of band center: 60% of cell diameter
k = 0.2*a;
%find positions tL1 and tL2 of upper/lower
tL1 = acos(k/a);
tL2 = 2*pi - acos(k/a);

syms t ts te;
%arc length integrand
gammaprime = sqrt((a*sin(t))^2+ (b*cos(t))^2);

%upper band starting point
tL1_s = double(vpasolve(0.5*arc-int(gammaprime,t,ts,tL1)==0,ts));
%upper band ending point
tL1_e = double(vpasolve(0.5*arc-int(gammaprime,t,tL1,te)==0,te));

%lower band starting point
tL2_s = double(vpasolve(0.5*arc-int(gammaprime,t,ts,tL2)==0,ts));
%lower band ending point
tL2_e = double(vpasolve(0.5*arc-int(gammaprime,t,tL2,te)==0,te));

tvec = [tL1_s, tL1_e, tL2_s, tL2_e];
end