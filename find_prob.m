function myProbs = find_prob(a,N1,env_d1,AP)
%calculates the binding probability based on dynein location

%Erin Angelini, 5.29.18

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%call parameters file
parameters

%now calculate binding probabilities
if strcmp(AP,'on')
    %anterior and posterior have different binding probabilities
    %anterior = 1, posterior = 0.65
    %AP-cuttoff is 60:40 P:A
    myProbs = zeros(1,N1);
    for i = 1:N1
        if a*cos(env_d1(i)) >= myCutoff %dynein in the anterior region
            myProbs(i) = P_p;
        else
            myProbs(i) = P_a;
        end
    end
elseif strcmp(AP,'off')
    %probability is one everywhere
    myProbs = ones(1,N1);
end

end