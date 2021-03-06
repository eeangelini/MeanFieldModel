function W = MainWithFilter(a, b, ctr, phi_1, phi_2, push, AP, LD, elas)
%for generating alpha vs. work plots efficiently
%THIS IS THE VERSION WITH A HIGH PASS FILTER - NOTE: THIS FILTER ONLY WORKS
%FOR A SYMMETRIC WELL

%inputs are ellipse major & minor axes a & b, spindle radius r_N, spindle
%center (x-coord) ctr, and upper & lower envelope spreads phi_1 & phi_2
%push is a string 'on' or 'off' that controls presence of the Let-99 band
%AP is a string 'on' or 'off' that controls presence of anterior to 
%posterior difference in MT binding probabilities
%LD is a string 'on' or 'off' that controls presence of MT-length dependent
%forces
%elas is a string 'on' or 'off' that controls presence of MT buckling


%Erin Angelini, 5.25.18

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load fixed parameters from parameters.m
parameters

%call alpha_halfplane for envelope endpoints above/below
alpha_halfplane(a, r_N, ctr, phi_1, phi_2)
%r_N is spindle radius from parameters.m

%get dynein locations
if a==b %circle
    gammaprime = @(t) sqrt(a^2*(sin(t)).^2 + a^2*(cos(t)).^2);
    total_arc = integral(gammaprime, 0, 2*pi, 'AbsTol',1e-12);
    N_d = round(total_arc); %want to keep interval of arc length ~ 1
    set_psi = equal_arcsCircle(a, N_d); %dynein locations
else %ellipse
    %get precalculated dynein locations
    if (a>=16 && a<=30) && b == b_basal %elongating x-axis, y-axis fixed
        cd 'mat files'
        load('dyneins16to30.mat')
        cd '../'
        set_psiCell = psiCell(a);
        set_psi = cell2mat(set_psiCell);
    elseif a==a_basal && (b>=16 && b<=25) %elongating y-axis, x-axis fixed
        cd 'mat files'
        load('dyneins_b16to35.mat')
        cd '../'
        set_psiCell = psibCell(b);
        set_psi = cell2mat(set_psiCell);
    elseif a==a_WTscale && b==b_WTscale %for volume scaling (wild type)
        cd 'mat files'
        load('dyneinsAR1p67VolScaleBy1p6')
        cd '../'
        %set_psi = scaleDynWithVol;
        %this is if we FIX N_d, comment out above line
        set_psi = fixDynAt128;
    else %or calculate them
        gammaprime = @(t) sqrt(a^2*(sin(t)).^2 + b^2*(cos(t)).^2);
        total_arc = integral(gammaprime, 0, 2*pi, 'AbsTol',1e-12);
        N_d = round(total_arc); %want to keep interval of arc length ~ 1
        set_psi = equal_arcsHannah(a, b, N_d); %calculate dyneins for ellipse
    end
end

%check for push bands
if strcmp(push,'on')
    %get Let-99 bands's locations
    tvec = let99(a, b, arc); %arc = Let-99 band arc length, in parameters file
elseif strcmp(push,'off')
    tvec = [];
end

%calculate energy landscape
W = [];
for i = 1:length(A) %A is vector of alpha values from parameters.m
    alpha = A(i);
    W = [W, work(alpha, a, b, set_psi, r_N, ctr, phi_1, phi_2, tvec, AP, LD, elas)];
end

%this is if we want to extract the max and min of W (for BVP solve_mfpt_new.m)
Wmax = max(W);
Wmin = min(W);

%plot the energy landscape
figure(1)
plot(A, W,'r-','LineWidth',4)
xlim([0 pi])
xticks([0 pi/4 pi/2 3*pi/4 pi])
xticklabels({'0','\pi/4','\pi/2','3\pi/4', '\pi'})
set(gca,'FontSize',30)
xlabel('\alpha')
ylabel('W(\alpha)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%do high-pass filter on energy landscape

% signal = W';
% %if a==b && ctr ==0
% %    lower_freq = 0.5;
% %    upper_freq = 0;
% %else
%     lower_freq = 0;
%     upper_freq = 0.5;
% %end
% 
% dt = 0.008;
% Fs = 1/dt;
% N = size(signal,1);
% dF = Fs/N;
% f = (-Fs/2:dF:Fs/2-dF)';
% 
% BPF = ((lower_freq < abs(f)) & (abs(f) < upper_freq));
% 
% signal=signal-mean(signal);
% 
% spektrum = fftshift(fft(signal))/N;
% 
% spektrum = real(BPF.*spektrum); 
% %returned as complex double but all 0 anyway?
% 
% signal=real(ifft(ifftshift(spektrum))); %inverse ifft
% %again complex part = 0 so????
% 
% if a~=b %only need to do for actual ellipse
%     signalShift = signal-(min(signal));
% else
%     signalShift = signal;
% end
% 
% figure(1)
% plot(A, signalShift, 'r-', 'LineWidth', 3)
% xlim([0 pi])
% if max(signalShift) ~= 0
%     ylim([min(signalShift) max(signalShift)+0.25*max(signalShift)])
% else
%     ylim([-0.5 0.5])
% end
% set(gca,'fontsize',15)
% xlabel('\alpha', 'FontSize', 18)
% ylabel('W(\alpha)', 'FontSize', 18)