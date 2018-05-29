function sinthetas1 = find_sintheta(alpha, a,b, r_N, ctr, t_q, t_p, N1,env_d1,set_psi,elas)
%finds sin(theta) of each MT in each envelope, where N is the number of MTs
%in the envelope, env_d is the vetor containing t_i's for dyneins, and
%ellipse parameters alpha/a/b/r_N are as usual
%t_q is lagging dynein end on envelope, t_p is leading

%elas is a string 'on' or 'off' that controls presence of MT buckling

%Erin Angelini, 5.8.18

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thetas1=[];
sinthetas1 = [];

%MTOC vector
R = [r_N*cos(alpha) + ctr, r_N*sin(alpha)];
%C=intersection of spindle axis and ellipse, ellipse(t_alpha) = C
if alpha == 0 || alpha == pi
    t_alpha = alpha;
else
    if ctr == 0
        if alpha <= pi/2 || alpha == 2*pi
            t_alpha = atan((a/b)*tan(alpha));
        elseif alpha > 3*pi/2 && alpha < 2*pi
            t_alpha = 2*pi + atan((a/b)*tan(alpha));
        else
            t_alpha = pi + atan((a/b)*tan(alpha));
        end
    else
        Ma = tan(alpha);
        Ba = -tan(alpha)*ctr;
        sols_a = [-a*(Ba*Ma*a+sqrt(Ma^2*a^2-Ba^2+b^2)*b)/(Ma^2*a^2+b^2), ...
            -a*(Ba*Ma*a-sqrt(Ma^2*a^2-Ba^2+b^2)*b)/(Ma^2*a^2+b^2)];
        for i = 1:2
            c1 = sols_a(i);
            if alpha >= 0 && alpha <= pi %C is above
                ta = acos(c1/a);
            else %C below
                ta = 2*pi - acos(c1/a);
            end
            C_check = [c1, b*sin(ta)];
            angle = acos( dot([ctr 0]-C_check, [ctr 0]-R)/(norm([ctr 0]-C_check)*norm([ctr 0] - R)) );
            if abs(angle) < 1
                t_alpha = ta;
            end
        end
    end
end
if not(exist('t_alpha', 'var'))
    alpha
end
C = [a*cos(t_alpha), b*sin(t_alpha)];

%vector b/t MTOC and ellipse
v1 = C - R;

%MT vectors
for i = 1:N1
    dynein = [a*cos(env_d1(i)), b*sin(env_d1(i))];
    v2 = dynein - R; %vector rep of MT
    %thetas1 = [thetas1 real(acos(dot(v1,v2)/(norm(v1)*norm(v2))))];
    thetas1 = [thetas1 atan2(abs(det([v1;v2])),dot(v1,v2))];
end

%shift thetas
for i=1:N1
    t_i = env_d1(i);
    if t_q>t_p && t_q > t_alpha %t_alpha in Q1, t_q in Q4
        if t_i < t_alpha || t_i > t_q %want theta negative
            if thetas1(i) > 0
                thetas1(i) = -thetas1(i);
            end
        end
    elseif t_q>t_p && t_q<t_alpha %t_alpha and t_q in Q4  
        if t_i < t_alpha && t_i > t_q 
            if thetas1(i) > 0
                thetas1(i) = -thetas1(i);
            end
        end
    else
        if t_i < t_alpha
            if thetas1(i) > 0
                thetas1(i) = -thetas1(i);
            end
        end
    end
end

%this is for including MT elasticity - perturb angle of attachment/pulling
if strcmp(elas,'on') 
    %update dynein vectors to include adjacent ones for envelope endpoints
    %determine these correctly
    if t_p < t_q %t_q is bigger i.e. envelope split across 0/2pi line
        k = 1;
        t_k = env_d1(k);
        while (t_k < t_q) && (k <= length(env_d1(k))) %get first dynein within t_q
            k = k + 1;
            t_k = env_d1(k);
        end
        oldend1 = t_k;
        for j = 1:length(env_d1)
            t_j = env_d1(j);
            if t_j <= t_p %get last dynein within t_p
                oldend2 = t_j;
            end
        end
    else
        oldend1 = env_d1(1);
        oldend2 = env_d1(length(env_d1));
    end
    k = find(set_psi == oldend1);
    j = find(set_psi == oldend2);
    if k == 1
        newend1 = set_psi(length(set_psi));
    else
        newend1 = set_psi(k-1);
    end
    if j == length(set_psi)
        newend2 = set_psi(1);
    else
        newend2 = set_psi(j+1);
    end
    envd = [newend1, env_d1, newend2];
    
    %now loop through all dyneins in envelope, perturb MT attachement
    %angles uniformly
    I = length(thetas1);
    for i = 2:I
        %vectors connecting adjacent dyneins
        d0 = envd(i); %location of dynein that the MT connects to
        D0 = [a*cos(d0), b*sin(d0)]; %actual vector
        d2 = envd(i-1);
        D2 = [a*cos(d2), b*sin(d2)] - D0; %vector to righthand dynein of d0
        d1 = envd(i+1);
        D1 = [a*cos(d1), b*sin(d1)] - D0; %vector to lefthand dynein of d0
        Phi = atan2(D1(2), D1(1)) - atan2(D2(2), D2(1)); %angle bt D1 & D2
        omega = Phi*rand(1,1); %random angle b/t 0 and Phi (uniform)
        %represents angle of buckled MT connection
        
        %unit tangent vector, use to find angel bt d2 and unit normal
        %unit normal is *RIGID* MT pulling direction, we want to update
        T = [-a*sin(d0) b*cos(d0)]./sqrt(a^2*(sin(d0))^2 + b^2*(cos(d0))^2); 
        %angle bt opp unit tangent and D2
        myAngle = atan2(D2(2), D2(1)) - atan2(T(2), T(1)); 
        if myAngle < 0
            myAngle = myAngle + 2*pi;
        end
        if myAngle > pi/2
            myAngle = 2*pi - myAngle;
        end
        psi = (pi/2) + myAngle; %angle bt unit normal MT vector and D2
        if omega <= psi %pulling vector directed to RHS of unit normal
            phi = psi - omega; %angle bt new (unit) pulling vector and unit normal
            thetas1(i) = thetas1(i) + phi; %ADD to old theta
        else %pulling vector directed to LHS of unit normal
            phi = omega - psi; %angle bt new (unit) pulling vector and unit normal
            thetas1(i) = thetas1(i) - phi; %SUBTRACT from old theta
        end
    end
end

%finally, calculate sins:
for i = 1:N1
    sinthetas1 = [sinthetas1 sin(thetas1(i))];
end

end