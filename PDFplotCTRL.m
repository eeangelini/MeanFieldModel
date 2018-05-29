function P = PDFplotCTRL(a)
%plots pdfs of spindle orientation at steady state from mean-field model;
%for the case of the on center, symmetric spindle
%takes in values of a from 16 to 30, b fixed at 15 (AR = a/b)

%Erin Angelini, 5.25.18

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get preset parameters
parameters

%load vector of k=k_1*W_max values
cd 'mat files'
load('kValsCtrlAR1p06to2.mat')
cd '../'
%now extract the value corresponding to AR = a/15
k = kvec(a-16);

W = 0.5*(1-cos(2.*A)); %unscaled functional form for control case
%A is range of alpha values from paramters.m

%now calclulate normalization constant
P_int = 0;
for i=2:length(A) %approximate integral over 0 to pi
    P_int = P_int + exp(-k.*W(i))*d_alpha;
end
N = 1/P_int;


P = N*exp(-k.*W);
plot(A, P, 'LineWidth', 2.5)
xlim([0 pi])
xticks([0 pi/4 pi/2 3*pi/4 pi])
xticklabels({'0','\pi/4','\pi/2','3\pi/4', '\pi'})
xlabel('\alpha', 'FontSize', 30)
ylabel('p(\alpha)', 'FontSize', 30)
end