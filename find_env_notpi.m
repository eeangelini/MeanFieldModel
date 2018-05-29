function [t_p, t_q] = find_env_notpi(alpha, a, b, r_N, ctr, phi, alphaq1, ...
    alphaq2, alphap1, alphap2)
%function that calculates endpoints of spindle envelope with angle phi ~= pi 
%given a spindle angle alpha, major/minor axes a and b, and spindle radius r_N

%Erin Angelini, 5.14.18

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if alpha == 2*pi
    alpha = 0; %same case and code being weird with 2*pi?
end

%find endpoints A & B of segment perpendicular to spindle
[t_a, t_b] = find_env_pi(alpha, a, b, r_N,ctr);
B = [a*cos(t_b), b*sin(t_b)]; 

%C=intersection of spindle axis and ellipse, ellipse(t_alpha) = C
R = [r_N*cos(alpha) + ctr, r_N*sin(alpha)]; %MTOC
if alpha == 0
    C = [a 0];
elseif alpha == pi
    C = [-a 0];
elseif abs(pi/2 - alpha) < 10e-6 || abs(3*pi/2 - alpha) < 10e-6
    if alpha >= 0 && alpha <= pi %C is above
        myY = b*sqrt(1 - (ctr/a)^2);
    else %C below
        myY = -b*sqrt(1 - (ctr/a)^2);
    end
    C = [ctr, myY];
else
    if ctr == 0
        if alpha <= pi/2
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
            angle = acos( dot([ctr 0]-C_check, [ctr 0]-R)/(norm([ctr 0]-C_check)*norm([ctr 0]-R)) );
            if abs(angle) < 10e-5
                t_alpha = ta;
            end
        end
    end
    C = [a*cos(t_alpha), b*sin(t_alpha)];
end



%define the desired vectors
CR = R-C;
CB = B-C;
arg = dot(CR, CB)/(norm(CR)*norm(CB));
gamma = acos(arg);
beta = pi - gamma - (phi/2);
x = norm(CR)*(sin(phi/2)/sin(beta));
y = norm(CB) - x;

%find D
D = find_end(B,C,x);
d_1 = D(1); d_2 = D(2);

%find the point E: reflect D across spindle axis
if alpha == 0 || alpha == pi
    E = [d_1, -d_2];
%else, get E by vector reflection and scaling
elseif abs(pi/2 - alpha) < 10e-6 || abs(3*pi/2 - alpha) < 10e-6
    %midpoint of line ED perpendicular to spindle axis
    mid1 = ctr; 
    mid2 = d_2;
    E(1) = 2*mid1 - d_1; %midpoint formula
    E(2) = 2*mid2 - d_2;
else
    mid1 = (cot(alpha)*d_1 + tan(alpha)*ctr + d_2)/(tan(alpha)+cot(alpha));
    mid2 = tan(alpha)*(mid1-ctr);
    E(1) = 2*mid1 - d_1; %midpoint formula
    E(2) = 2*mid2 - d_2;
end

%now we know the lines defd by RE and RD
Mq = (D(2)-R(2))/(D(1)-R(1));
Bq = R(2) - R(1)*Mq;
Mp = (E(2)-R(2))/(E(1)-R(1));
Bp = R(2) - R(1)*Mp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%symbolic solns from Maple
solsQ = [-a*(Bq*Mq*a+sqrt(Mq^2*a^2-Bq^2+b^2)*b)/(Mq^2*a^2+b^2), ...
    -a*(Bq*Mq*a-sqrt(Mq^2*a^2-Bq^2+b^2)*b)/(Mq^2*a^2+b^2)];
solsP = [-a*(Bp*Mp*a+sqrt(Mp^2*a^2-Bp^2+b^2)*b)/(Mp^2*a^2+b^2),...
    -a*(Bp*Mp*a-sqrt(Mp^2*a^2-Bp^2+b^2)*b)/(Mp^2*a^2+b^2)];

for i = 1:2
    q1 = solsQ(i);
    if alpha >= alphaq1 && alpha <= alphaq2 %q above
        tq = acos(q1/a);
    else %q below
        tq = 2*pi - acos(q1/a);
    end
    for j = 1:2
        p1 = solsP(j);
        if alpha > alphap1 || alpha <= alphap2 %p above
            tp = acos(p1/a);
        else %p below
            tp = 2*pi - acos(p1/a);
        end
        %check phi
        Q = [a*cos(tq) b*sin(tq)]; P = [a*cos(tp) b*sin(tp)];
        vec1 = R - P; vec2 = R - Q;
        phi_check = acos(dot(vec2, vec1)/(norm(vec2)*norm(vec1)));
        err = abs(phi_check - phi);
        if err < 10e-4
            if alpha > 2*pi - acos(r_N/(a-ctr)) && alpha <= alphap1
                %A above but P below in right half plane
                if alpha < acos(r_N/(a-ctr)) && alpha >= alphaq1
                    %Q above but B below in right half plane
                    if (tp <= t_a +2*pi) && (tq+2*pi >= t_b)
                        t_q = tq;
                        t_p = tp;
                    end
                else
                    if (tp <= t_a+2*pi) && (tq >= t_b)
                        t_q = tq;
                        t_p = tp;
                    end
                end
            else
                if alpha < acos(r_N/(a-ctr)) && alpha >= alphaq1
                    %Q above but B below in right half plane
                    if (tp <= t_a) && (tq+2*pi >= t_b)
                        t_q = tq;
                        t_p = tp;
                    end
                else
                    if (tp <= t_a) && (tq >= t_b)
                        t_q = tq;
                        t_p = tp;
                    end
                end
            end
        end
    end
end 

%for debugging
if not(exist('t_p', 'var')) || not(exist('t_q', 'var'))
    alpha;
    %t_q = tq;
    %t_p = tp;
end

end