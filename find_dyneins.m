function [env_d, let] = find_dyneins(b, t_p, t_q, set_psi, tvec)
%takes regular ellipse inputs plus parameter's t_p&t_q for leading (P) and 
%lagging (Q) endpoints of the spindle envelope, as well as the vector of 
%dynein t_k's set_psi, and returns the vector of t_i's for the dyneins 
%countained in the envelope
%also identifies which locations overlap with the Let-99 band(s), sorts
%them out and into the let vector

%Erin Angelini, 6.21.17

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

env_d = []; let =[];
%looking for cases P and Q NOT divided by line {y=0, x>0} i.e. alpha=0=2*pi
%happens when Q is below and P is above

q2 = b*sin(t_q);
p2 = b*sin(t_p);

if not(isempty(tvec)) %include Let99 push band
    tL1_s = tvec(1); tL1_e = tvec(2); tL2_s = tvec(3); tL2_e = tvec(4);
    if p2 > 0 && q2 < 0
        for k = 1:length(set_psi)
            t_k = set_psi(k);
            if t_k <= t_p && t_k >=0
                if (t_k >= tL1_s && t_k <= tL1_e) || ...
                        (t_k >= tL2_s && t_k <= tL2_e) %location in a Let-99 band
                    let = [let, t_k];
                else
                    env_d = [env_d t_k];
                end
            elseif t_k >= t_q && t_k <= 2*pi
                if (t_k >= tL1_s && t_k <= tL1_e) || ...
                        (t_k >= tL2_s && t_k <= tL2_e) %location in a Let-99 band
                    let = [let, t_k];
                else
                    env_d = [env_d t_k];
                end
            end
        end
    else 
        for k = 1:length(set_psi)
            t_k = set_psi(k);
            if t_k <= t_p && t_k >= t_q
                if (t_k >= tL1_s && t_k <= tL1_e) || ...
                        (t_k >= tL2_s && t_k <= tL2_e) %location in a Let-99 band
                    let = [let, t_k];
                else
                    env_d = [env_d t_k];
                end
            end
        end
    end
else
    if p2 > 0 && q2 < 0
        for k = 1:length(set_psi)
            t_k = set_psi(k);
            if t_k <= t_p && t_k >=0
                env_d = [env_d t_k];
            elseif t_k >= t_q && t_k <= 2*pi
                env_d = [env_d t_k];
            end
        end
    else 
        for k = 1:length(set_psi)
            t_k = set_psi(k);
            if t_k <= t_p && t_k >= t_q
                env_d = [env_d t_k];
            end
        end
    end
end

end
