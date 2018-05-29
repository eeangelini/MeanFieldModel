function myForces = find_pull(a,b,r_N,ctr,alpha,N1,env_d1,LD)
%calculates the MT-length dependent pulling forces

%Erin Angelini, 5.29.18

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%call parameters file
parameters

if strcmp(LD,'on')
    %length dependent pulling forces: f(L) = L^3
    myForces = zeros(1,N1);
    R = [r_N*cos(alpha) + ctr, r_N*sin(alpha)]; %MTOC vector
    for i = 1:N1
        t = env_d1(i); %dynein location
        D = [a*cos(t), b*sin(t)]; %MT-dynein contact point
        M = D - R; %MT vector
        L = length(M); %MT length
        myForces(i) = L^beta;
    end
elseif strcmp(LD,'off')
    %pulling force is one everywhere
    myForces = ones(1,N1);
end

end