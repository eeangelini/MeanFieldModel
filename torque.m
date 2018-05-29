function T = torque(alpha, a, b, set_psi, r_N, ctr, phi_1, phi_2, tvec, AP, LD, elas)
%function that plots the spindle setup with given spindle angle of
%orientation alpha, major/minor axes a and b of ellipse, set_psi positions
%(in terms of t) of dyneins

%for varying center and spindle envelopes: also takes spindle radius r_N,
%center ctr, and upper/lower envelope spreads phi_1 and phi_2

%tvec has Let-99 band locations - see workplot.m for function call

%AP, LD, and elas are strings that control presence of anterior to posterior 
%difference in MT binding probabilities (AP), presence of MT-length 
%dependent forces, and presence of MT buckling (respectively)

%Erin Angelini, 5.25.18

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load parameters specific to torque.m
parameters

%calibration for input angle
if alpha > 2*pi
    alpha = alpha - 2*pi;
elseif alpha < 0
    alpha = alpha + 2*pi;
end

%find coordinates of upper envelope endpoints
%load variables from previous call of alpha_halfplane.m:
load('halfplane1.mat')
if phi_1 == phi_pi
    [t_p, t_q] = find_env_pi(alpha, a, b, r_N,ctr);
else
    [t_p, t_q] = find_env_notpi(alpha,a,b,r_N, ctr, phi_1,...
        alphaq1, alphaq2, alphap1, alphap2);
end

%identify Let-99 band and t_k's of dyneins within envelope
[env_d1, let1] = find_dyneins(b, t_p, t_q, set_psi, tvec);
N1 = length(env_d1);
Nlet1 = length(let1);

%now do same for lower envelope
%find coordinates of lower envelope endpoints; v is leading, u is lagging
%load variables from previous call of alpha_halfplane.m:
load('halfplane2.mat')
if phi_2 == phi_pi
    if alpha <= pi
        [t_v, t_u] = find_env_pi(alpha + pi, a, b, r_N,ctr);
    else
        [t_v, t_u] = find_env_pi(alpha - pi, a, b, r_N,ctr);
    end
else
    if alpha <= pi
        [t_v, t_u] = find_env_notpi(alpha + pi, a, b, r_N, ctr,phi_2,... 
            alphaQ1, alphaQ2, alphaP1, alphaP2);
    else
        [t_v, t_u] = find_env_notpi(alpha - pi, a, b, r_N, ctr,phi_2,... 
            alphaQ1, alphaQ2, alphaP1, alphaP2);
    end
end

%find dyneins in lower envelope
[env_d2, let2] = find_dyneins(b, t_v, t_u, set_psi, tvec);
N2 = length(env_d2);
Nlet2 = length(let2);
if N2 == 0 %check point for error in t_v, t_u
    [env_d2, let2] = find_dyneins(b, t_u, t_v, set_psi, tvec);
    N2 = length(env_d2);
    Nlet2 = length(let2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%now use t_i's to calculate sin(theta_i)
sinthetas1 = find_sintheta(alpha, a,b, r_N, ctr, t_q, t_p, N1, env_d1,set_psi,elas);

if alpha < pi
    sinthetas2 = find_sintheta(alpha+pi, a,b, r_N, ctr, t_u, t_v, N2, env_d2, set_psi, elas);
else
    sinthetas2 = find_sintheta(alpha-pi, a,b, r_N, ctr, t_u, t_v, N2, env_d2, set_psi, elas);
end

%separately for Let-99 MTs
if not(isempty(let1))
    tL1_s = tvec(1); tL1_e = tvec(2);
    sinthetaslet1 = find_sintheta(alpha, a,b, r_N, ctr, tL1_s, tL1_e, Nlet1, let1, set_psi, elas);
end

if not(isempty(let2))
    if alpha <= pi
        tL2_s = tvec(3); tL2_e = tvec(4);
        sinthetaslet2 = find_sintheta(alpha+pi, a,b, r_N, ctr, tL2_s, tL2_e, Nlet2, let2, set_psi, elas);
    else
        tL2_s = tvec(3); tL2_e = tvec(4);
        sinthetaslet2 = find_sintheta(alpha-pi, a,b, r_N, ctr, tL2_s, tL2_e, Nlet2, let2, set_psi, elas);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate MT binding probabilities
pvec1 = find_prob(a,N1,env_d1,AP); %upper envelope

pvec2 = find_prob(a,N2,env_d2,AP); %lower envelope

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calculate MT-dynein pulling forces
Fpull1 = find_pull(a,b,r_N,ctr,alpha,N1,env_d1,LD); %upper envelope

Fpull2 = find_pull(a,b,r_N,ctr,alpha,N2,env_d2,LD); %lower envelope

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%finally, calculate torque

%MT density: #MTs per angle in radians; values from paramters.m
if phi_1 > phi_2 %asymmetric envelope, top aster larger
    rho1 = rho_asymm;
else %symmetric envelopes
    rho1 = rho_symm;
end

if phi_2 > phi_1 %asymmetric envelope, bottom aster larger
    rho2 = rho_asymm;
else %symmetric envelopes
    rho2 = rho_symm;
end


%sum sin(theta_i)*(p_i.*Fpull_i) <- torque term for dynein_i
if isempty(pvec1) || isempty(Fpull1)
    alpha
else
    sum1 = dot(sinthetas1,(pvec1.*Fpull1)); 
end
if isempty(pvec2) || isempty(Fpull2)
    alpha
else
    sum2 = dot(sinthetas2,(pvec2.*Fpull2)); 
end


T1 = r_N*rho1*sum1; 
T2 = r_N*rho2*sum2;
Fpush = 1; %for now, pushing forces=1
if exist('sinthetaslet1', 'var')
    sum3 = sum(sinthetaslet1);
    T3 = r_N*Fpush*rho1*sum3;
else
    T3 = 0;
end
if exist('sinthetaslet2', 'var')
    sum4 = sum(sinthetaslet2);
    T4 = r_N*Fpush*rho2*sum4;
else
    T4 = 0;
end
T=T1+T2+T3+T4; %global torque on PNC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOTTING CODE

% %plot ellipse
% ellipse(0,0,a,b)
% hold on
% 
% %plot dyneins
% figure(1)
% plot(a*cos(set_psi),b*sin(set_psi),'o', 'MarkerSize', 7,...
%     'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b')
% hold on
% % 
% 
% %plot Let-99 bands
% if not(isempty(tvec))
%     tband1 = tL1_s:0.1:tL1_e; tband2 = tL2_s:0.1:tL2_e;
%     plot(a*cos(tband1), b*sin(tband1), 'r-', 'LineWidth', 10)
%     hold on
%     plot(a*cos(tband2), b*sin(tband2), 'r-', 'LineWidth', 10)
% end
% 
% %envelope boundaries
% mtoc1 = [r_N*cos(alpha)+ctr r_N*sin(alpha)];
% mtoc2 = [-r_N*cos(alpha)+ctr -r_N*sin(alpha)];
% 
% plot([a*cos(t_p) mtoc1(1)], [b*sin(t_p) mtoc1(2)],'r--')
% hold on
% plot([a*cos(t_q) mtoc1(1)], [b*sin(t_q) mtoc1(2)],'r--')
% hold on
% plot([a*cos(t_v) mtoc2(1)], [b*sin(t_v) mtoc2(2)],'r--')
% hold on
% plot([a*cos(t_u) mtoc2(1)], [b*sin(t_u) mtoc2(2)],'r--')
% hold on
% 
% %vectors from MTOC to each dynein
% for i = 1:N1
%    dynein1 = [a*cos(env_d1(i)) b*sin(env_d1(i))];
%    plot([mtoc1(1) dynein1(1)], [mtoc1(2) dynein1(2)], 'g-', 'LineWidth', 2)
%    hold on
% end
% for i = 1:N2
%    dynein2 = [a*cos(env_d2(i)) b*sin(env_d2(i))];
%    plot([mtoc2(1) dynein2(1)], [mtoc2(2) dynein2(2)], 'g-', 'LineWidth', 2)
%    hold on
% end
% 
% %and to Let-99 bands
% if not(isempty(tvec))
% for i = 1:Nlet1
%    mt1 = [a*cos(let1(i)) b*sin(let1(i))];
%    plot([mtoc1(1) mt1(1)], [mtoc1(2) mt1(2)], 'g-', 'LineWidth', 1.5)
%    hold on
% end
% for i = 1:Nlet2
%    mt2 = [a*cos(let2(i)) b*sin(let2(i))];
%    plot([mtoc2(1) mt2(1)], [mtoc2(2) mt2(2)], 'g-', 'LineWidth', 1.5)
%    hold on
% end
% end
% 
% %plot spindle axis
% plot([ctr r_N*cos(alpha) + ctr], [0 r_N*sin(alpha)], 'Color', [.95,.90,0], 'LineWidth', 3)
% hold on
% plot([ctr -r_N*cos(alpha) + ctr], [0 -r_N*sin(alpha)], 'Color', [.95,.90,0], 'LineWidth', 3)
% hold on
% 
% %plot spindle center and MTOC
% plot(ctr, 0, 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 5)
% plot(r_N*cos(alpha)+ctr, r_N*sin(alpha), 'o','MarkerFaceColor',...
%     [.95,.90,0], 'MarkerEdgeColor', [.95,.90,0], 'MarkerSize', 10)
% plot(-r_N*cos(alpha)+ctr, -r_N*sin(alpha), 'o','MarkerFaceColor',...
%     [.95,.90,0], 'MarkerEdgeColor', [.95,.90,0], 'MarkerSize', 10)
% 
% pbaspect([a/b 1 1])

end
