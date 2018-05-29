function set_psi = equal_arcsCircle(a, N_d)
%function that places dyneins at equal arc length in the case of a circle
%with radius a, # dyneins = num_dyn

%Erin Angelini, 6.21.17

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

total_arc = 2*pi*a; %circumference
L = total_arc/N_d; %constant arc length between each dynein
theta = L/a; %central angle between each pair of adjacent dyneins:2pi/N_d
set_psi = 0; %first dynein placed at 0
for i = 2:N_d
    set_psi = [set_psi, (i-1)*theta];
end

end