function D = find_end(B,C,x)
%takes points B and C (1x2 vectors) as well as (scalar) length x and finds 
%the point D = [d_1, d_2] with distance x from C along the line defined by 
%C & B

%Erin Angelini, 5.28.18

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find vector xHat with magnitude x in direction of CB
CB = B-C;
normCB = CB./norm(CB);
xHat = x.*normCB;
%translate w.r.t the point C (gives us the point D, D - C = xHat)
D = xHat + C;

end