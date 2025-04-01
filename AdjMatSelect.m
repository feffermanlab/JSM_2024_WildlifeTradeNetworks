function [W, N, r]= AdjMatSelect(optW, optr,l)
%ADJMATSELECT Auxilary function to select the adjacency matrix
%   optW: the adjacency matrix W identified by an integer 1-6
%   optR: the vector of intrinsic benefites identified by an integer 1-3
%   l: constant multiplier for intrinsic weight default to 1

%An auxilary function which keeps all the example adjacency matrices
%easily accessible

if optW == 1
    W = [0,0,0,0,0,0,0,0;
            0,0,0,0,0,0,0,0;
            0.5,0.5,0,0,0,0,0,0;
            0.5,0.5,0,0,0,0,0,0;
            0,0,1,0,0,0,0,0;
            0,0,0.5,0.5,0,0,0,0;
            0,0,0.5,0.5,0,0,0,0;
            0,0,0,1,0,0,0,0];
elseif optW == 2
    W = [0,0,0,0,0,0,0,0;
            0,0,0,0,0,0,0,0;
            0.5,0.5,0,0,0,0,0,0;
            0.5,0.5,0,0,0,0,0,0;
            0,0,1,0,0,0,0,0;
            0,0,1,0,0,0,0,0;
            0,0,1,0,0,0,0,0;
            0,0,0.5,0.5,0,0,0,0];
elseif optW == 3
    W = [0,0,0,0,0,0,0,0;
            0,0,0,0,0,0,0,0;
            1,0,0,0,0,0,0,0;
            0.5,0.5,0,0,0,0,0,0;
            0,0,1,0,0,0,0,0;
            0,0,0.5,0.5,0,0,0,0;
            0,0,0.5,0.5,0,0,0,0;
            0,0,0,1,0,0,0,0];
elseif optW==4
    W = [0,0,0,0,0,0,0,0;
            0,0,0,0,0,0,0,0;
            1,0,0,0,0,0,0,0;
            0.5,0.5,0,0,0,0,0,0;
            0,0,1,0,0,0,0,0;
            0,0,1,0,0,0,0,0;
            0,0,1,0,0,0,0,0;
            0,0,0.5,0.5,0,0,0,0];
elseif optW==5
    W=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
        1/2,1/2,0,0,0,0,0,0,0,0,0,0,0,0,0;
        1/2,1/2,0,0,0,0,0,0,0,0,0,0,0,0,0;
        1/2,1/2,0,0,0,0,0,0,0,0,0,0,0,0,0;
        0,0,1/3,1/3,1/3,0,0,0,0,0,0,0,0,0,0;
        0,0,1/3,1/3,1/3,0,0,0,0,0,0,0,0,0,0;
        0,0,1/3,1/3,1/3,0,0,0,0,0,0,0,0,0,0;
        0,0,1/3,1/3,1/3,0,0,0,0,0,0,0,0,0,0;
        0,0,1/3,1/3,1/3,0,0,0,0,0,0,0,0,0,0;
        0,0,1/3,1/3,1/3,0,0,0,0,0,0,0,0,0,0;
        0,0,1/3,1/3,1/3,0,0,0,0,0,0,0,0,0,0;
        0,0,1/3,1/3,1/3,0,0,0,0,0,0,0,0,0,0;
        0,0,1/3,1/3,1/3,0,0,0,0,0,0,0,0,0,0;
        0,0,1/3,1/3,1/3,0,0,0,0,0,0,0,0,0,0;
        ];
elseif optW==6
     W=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
        1/2,1/2,0,0,0,0,0,0,0,0,0,0,0,0,0;
        1/2,1/2,0,0,0,0,0,0,0,0,0,0,0,0,0;
        1/2,1/2,0,0,0,0,0,0,0,0,0,0,0,0,0;
        0,0,1/2,1/2,0,0,0,0,0,0,0,0,0,0,0;
        0,0,1/2,1/2,0,0,0,0,0,0,0,0,0,0,0;
        0,0,0,1,0,0,0,0,0,0,0,0,0,0,0;
        0,0,0,1,0,0,0,0,0,0,0,0,0,0,0;
        0,0,0,1,0,0,0,0,0,0,0,0,0,0,0;
        0,0,0,1,0,0,0,0,0,0,0,0,0,0,0;
        0,0,0,1,0,0,0,0,0,0,0,0,0,0,0;
        0,0,0,1,0,0,0,0,0,0,0,0,0,0,0;
        0,0,0,1/2,1/2,0,0,0,0,0,0,0,0,0,0;
        0,0,0,1/2,1/2,0,0,0,0,0,0,0,0,0,0;
        ];
else
    W=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
        1/2,1/2,0,0,0,0,0,0,0,0,0,0,0,0,0;
        1/2,1/2,0,0,0,0,0,0,0,0,0,0,0,0,0;
        1/2,1/2,0,0,0,0,0,0,0,0,0,0,0,0,0;
        0,0,0,1,0,0,0,0,0,0,0,0,0,0,0;
        0,0,0,1,0,0,0,0,0,0,0,0,0,0,0;
        0,0,0,1,0,0,0,0,0,0,0,0,0,0,0;
        0,0,0,1,0,0,0,0,0,0,0,0,0,0,0;
        0,0,0,1,0,0,0,0,0,0,0,0,0,0,0;
        0,0,0,1,0,0,0,0,0,0,0,0,0,0,0;
        0,0,0,1,0,0,0,0,0,0,0,0,0,0,0;
        0,0,0,1,0,0,0,0,0,0,0,0,0,0,0;
        0,0,0,1,0,0,0,0,0,0,0,0,0,0,0;
        0,0,0,1,0,0,0,0,0,0,0,0,0,0,0;
        ];

end
N = length(W(1,:));
mult = 1;
if nargin == 3
    mult = l;
end
if optr == 1
    r = mult*[1,1,2,2,4,4,4,4];
elseif optr==2
    r = mult*[1,1,1,1,1,1,1,1];
elseif optr==3
    r = mult*[1,1,2,2,2,4,4,4,4,4,4,4,4,4,4];
else
    r = mult*[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
end
end