function set_psi = equal_arcsHannah(a, b, N_d)
%EA: this function (written by Hannah Wayment-Steele in 2016) calculates
%the dynein locations along the cortex for a given major axis a, minor axis
%b, and number of dyneins N_d

%Starts at 0 and goes to 2*pi, giving the angles corresponding to equal arc
%lengths. Each arc length unit is determined by total_arc/num_dynein.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%"defaults":
%a = 25;
%b = 15;
%num_dyn = 128
%total_arc = 127.64;

gammaprime = @(t) sqrt(a^2*(sin(t)).^2 + b^2*(cos(t)).^2);
total_arc = integral(gammaprime, 0, 2*pi, 'AbsTol',1e-12); %arc length of ellipse

arc_integrand = @(x) a.*b.*sqrt((b.^4.*cos(x).^2+a.^4*sin(x).^2)./(b.^2.*cos(x).^2+a.^2*sin(x).^2).^3);

L=total_arc/N_d; %this is the fixed arc length between dyneins

q=0;
psi_n = 0.0005; %psi step size 
 
psi = 0; %position of first dynein
set_psi = [psi];

%now we calculate the dynein locations by calculating arc length between
%current dynein location psi to some position psi_n; stop when arc length
%is equal to L
%end total loop once we have swept out entire boundary of the ellipse
while psi < 2*pi
    while q < L
        q = integral(arc_integrand,psi,psi_n,'AbsTol',1e-12);
        psi_n = psi_n+0.0005;
    end
        
    psi = psi_n;
    set_psi = [set_psi psi];
    psi_n = 0.001;
    q = 0;
    
end

end
