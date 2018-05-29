function [t_p, t_q] = find_env_pi(alpha, a, b, r_N, ctr)
%function that calculates endpoints of spindle envelope (phi = pi) given a
%spindle angle alpha, major/minor axes a and b, spindle radius r_N and
%spindle x-coordinate ctr

%Erin Angelini, 6.19.17

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = -cot(alpha);
B = r_N/sin(alpha) + cot(alpha)*ctr;
R = [r_N*cos(alpha) + ctr, r_N*sin(alpha)]; %MTOC

alphaq1 = acos(r_N/(a-ctr)); alphaq2 = acos(r_N/(a+ctr)) +pi;
alphap1 = 2*pi - acos(r_N/(a-ctr)); alphap2 = pi - acos(r_N/(a+ctr));
    
if alpha == 0 || alpha == 2*pi %envelope line is x = r_N + ctr
    t_p = acos((r_N+ctr)/a);
    t_q = 2*pi - acos((r_N+ctr)/a);
elseif abs(alpha - pi/2) < 10e-6 %line is y = r_N
    t_p = pi - asin(r_N/b);
    t_q = asin(r_N/b);
elseif alpha == pi %line is x = -r_N + ctr
    t_p = 2*pi - acos((-r_N+ctr)/a);
    t_q = acos((-r_N+ctr)/a);
elseif abs(alpha - 3*pi/2) < 10e-6 %line is y = -r_N
    t_p = 2*pi - asin(r_N/b);
    t_q = pi+asin(r_N/b);
else
    sols = [-a*(B*M*a+sqrt(M^2*a^2-B^2+b^2)*b)/(M^2*a^2+b^2),...
        -a*(B*M*a-sqrt(M^2*a^2-B^2+b^2)*b)/(M^2*a^2+b^2)];
    for i = 1:2
        q1 = sols(i);
        if alpha >= alphaq1 && alpha <= alphaq2 %q above above x-axis
            tq = acos(q1/a);
        else %q below x-axis
            tq = 2*pi - acos(q1/a);
        end
        for j = 1:2
            p1 = sols(j);
            if alpha > alphap1 || alpha <= alphap2 %p above
                tp = acos(p1/a);
            else
                tp = 2*pi - acos(p1/a);
            end
            %check that angle spread is pi
            Q = [a*cos(tq) b*sin(tq)]; P= [a*cos(tp) b*sin(tp)];
            angle = acos(dot(R-Q,R-P)/(norm(R-Q)*norm(R-P)));
            err = abs(angle - pi);
            if err < 10e-4
                t_q = tq;
                t_p = tp;
            end
        end
    end
        
    %want to check that t_p and t_q are in correct order
    if (alpha < alphaq1 || alpha > alphaq2) && ...
            (alpha > alphap1 || alpha <= alphap2) %q below, p above
        %should have t_q > t_p
        if t_q < t_p
            temp_t = t_q;
            t_q = t_p;
            t_p = temp_t;
        end
    else %should have t_q < t_p
        if t_q > t_p
            temp_t = t_q;
            t_q = t_p;
            t_p = temp_t;
        end
    end
end

end