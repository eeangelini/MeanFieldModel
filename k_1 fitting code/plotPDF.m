function resnorm = plotPDF(psi,Wmax)
%takes final data from simulations and plots normalized histogram & PDF of
%the data against angle of orientation on x-axis
%uses result Wmax(AR) from mean-field model (max of the energy landscape
%from Main.m) to fit the data to the closed form of the energy landscape
%vector of final data comes from .csv file

%Erin Angelini, 5.29.18

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%histogram version of PDF
figure(1)
[t,N,X] = nhist(psi,'f', 20, 'smooth', 'xmin', 0, 'xmax', 2*pi, 'xlabel', 'Final Orientation');
xlim([0,pi])

%normalized line plot version of PDF
myGuess = round(length(X)/2); %bin containing pi near middle
k = [myGuess-2, myGuess-1, myGuess, myGuess+1, myGuess+2];
stopIndex = 0;
for i = 1:length(k)
    if X(k(i)) <= pi
        stopIndex = k(i);
    end
end
Xnew = X(1:stopIndex);

P = [];
for i = 1:stopIndex
    P = [P, N(i)/max(N(1:stopIndex))]; %for normalizing
    %P = [P, N(i)];
end
plot(Xnew, P, 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b')
xlim([0,pi])
xticks([0 pi/4 pi/2 3*pi/4 pi])
xticklabels({'0','\pi/4','\pi/2','3\pi/4', '\pi'})
set(gca,'FontSize',30)
xlabel('Final angle \alpha')


%for fitting parameter
x0=0; %arbitrary start point
%NOTE: this is energy landscape for *symmetric wells*, change for asymm
W = @(Xnew) Wmax*0.5*(1-cos(2.*(Xnew))); 
F = @(k1,Xnew) exp(-k1*(W(Xnew))); 
[x,resnorm] = lsqcurvefit(F,x0,Xnew,P);
hold on
plot(Xnew, F(x,Xnew), 'm-', 'LineWidth', 2)
legend('Simulation Data', 'Matlab Fit')
mystr1 = ['k_1 = ', num2str(x)];
mystr2 = ['MSE = ', num2str(resnorm)];
text(2, 0.85, mystr1,'FontSize',12)
text(2, 0.8, mystr2,'FontSize',12)
ylabel('P(\alpha)')
end
