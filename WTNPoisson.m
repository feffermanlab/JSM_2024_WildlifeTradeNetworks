%%This script runs a numerical example for the downstream contagion with
%   with infection determined by a random variable rather than the binary
%   process of the original manuscipt. 

%Select an adjacency matrix
%[W,N,r]=AdjMatSelect(2,1);
W= [0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0;
    0.5,0.5,0,0,0,0,0,0;
    0.5,0.5,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0;
    0,0,1,0,0,0,0,0;
    0,0,1,0,0,0,0,0;
    0,0,0.5,0.5,0,0,0,0];

%Save number of players.
N=length(W(1,:));

%risk of contamination by environment
el=0.2;

%initial strategy profile
x=rand(N,1);
%initial thresholds
%taus = 0.09*ones(N,1);%...................................Constant
taus = [0.2,0.2,0.1,0.1,0.05,0.05,0.05,0.05]';%............Increasing


disp("infection load")
pstar(x,N,W,el)
disp("payoff")
payoffs(x,N,W,taus,el)
disp('best responses')
MBR(x,N,W,taus,el)
disp('Nash Equilibrium')
MBRseqence(x,N,W, taus,el,10)


function eq = MBRseqence(x,N,W,taus,el, iter)
%MBRsequence simulates a myopic best response sequence which, under certain
%conditions is guarenteed to converge
%   x:      initial strategy profile
%   N:      Number of Players
%   W:      Adjacency Matrix
%   taus:   Threshold Values 
%   el:     environmental introduction probability
%   iter:   Iteration limit on MBR

%This function simply forms the sequence of myopic best responses of size
%iter starting from the intitial strategy profile x. In the case where MBR
%has a unique fixed point, this function is sure to converge to the unique
%Nash equilibrium for sufficiently large iter. 

    for i =1:iter
        x = MBR(x,N,W,taus,el);
    end
    eq=x;
end


function xprime = MBR(x,N,W,taus,el)
%MBR computes the best response for every player against a given strategy
%profile
%   x:      Given strategy profile
%   N:      Number of Players
%   W:      Adjacency Matrix
%   taus:   Threshold Values 
%   el:     environmental introduction probability

%For each player, this function finds the strategy which maximizes their
%payoff assuming that all other players do not change their strategy. When
%this mapping from strategy profiles to strategy profiles satisfied the
%conditions of Brouwer's Fixed Point Theorem, then the mapping has a unique
%fixed point
    xprime= zeros(N,1);
    for idx =1:N
        fun=@(xi)-1*focalPayoff(xi,x,N,W,taus,el,idx);
        br= fminbnd(fun,0,1);
        xprime(idx)=br;
    end
end

function wi = focalPayoff(focalx,x,N,W,taus,el,focalidx)
%Focal Payoff computes the payoff for a focal individual with a specified
%strategy when all other players play according to the strategy profile x
%   focalx: Strategy of focal individual
%   x:      Strategy profile of all other players
%   N:      Number of Players
%   W:      Adjacency Matrix
%   taus:   Threshold Values 
%   el:     environmental introduction probability
%   focalidx:   The index of the focal individual.

% This helper function replaces the focal individuals strategy from x with 
% a preselected strategy, focalx, then computes the payoff for the focal 
% player 
    x(focalidx)=focalx;
    ws = payoffs(x,N,W,taus,el);
    wi=ws(focalidx);
end


function ws = payoffs(x,N,W,taus,el)
%Payoffs computes the payoff for every player given the specified strategy
%profils x
%   x:      Strategy profile
%   N:      Number of Players
%   W:      Adjacency Matrix
%   taus:   Threshold Values 
%   el:     environmental introduction probability

%This function computes the payoff for every player using the formulation
%which involves a payout from a successful sale, a loss from an unsuccesful
%sale, and a cost associated to the strategy x
   
    %Cost function
    C=0.8*-x.^2;
    %Loss from unsuccessful sale
    L=(0.5-x).^2;
    %payout from successful sale
    R=ones(N,1);
    
    %Compute stationary loads
    lambdas = pstar(x,N,W,el);
    %Compute probabilites for succesful sales
    probs = poisscdf(taus,lambdas);
    %Calculate payoffs
    ws= (L+R).*probs+C-L;
end


function p = pstar(x,N,W,el)
%PSTAR computes the stationary probability of infection for a given network
%and strategy profile
%   x: a strategy profile
%   N: number of players
%   W: Adjacency matrix
%   el: Probability of infection from the environment

    %The pass through rate associated to each pairwise interaction
    %Here the pass through rate is modeled in the following way but it can
    %be changed:
    %   A transaction from i to j may happen with probability W_{i,j}
    %   That transaction probability is modulated by i and j's strategies
    %       If j uses a high strategy, his decision to purchase from i will
    %       be affected linearly by i's strategy, if not, his diecision
    %       will not depend on i's strategy
    %   If a transaction happens, the probability of transmission decreases
    %   with i's strategy and with j's strategy
    A=W.*(x*x'+(1-x*ones(1,N))).*((1-x)*(1-x)');

    %The environmental contamination rate is simply given by a constant
    %multiple of 1-x
    B=el.*(1-x);

    %The stationary load is calculated a described in the manuscript. 
    p=(eye(N)-A)\B;
end