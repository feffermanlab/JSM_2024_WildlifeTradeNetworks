%%This script runs a numerical example for the downstream contagion with
%   with infection determined by a random variable rather than the binary
%   process of the original manuscipt. 

%example() is a function that demonstrates what the functions in this
%script can do
%example()


%Numerical Experiment 1
%NumericalExperiment1()

%Numerical Experiment 2
%NumericalExperiment2()

function ret = NumericalExperiment1()
    %subsection4.1
    
    %environmental introduction rate
    el=0.2;
    %initial thresholds
    %taus = 0.09*ones(N,1);%...................................Constant
    taus = [0.2,0.2,0.1,0.1,0.05,0.05,0.05,0.05]';%............Increasing
    

    %get a symmetric and ansymmetric adjacency matrix and find equilibrium
    %strategies through best response simulation
    [W1,N1, r1]=AdjMatSelect(1,1,2);
    [W2,N2, r2]=AdjMatSelect(2,1,2);
    res1=simulate(N1,W1,taus,el,100,100);
    res2=simulate(N2,W2,taus,el,100,100);
 
    %display results for the symmetric case 
    disp("symmetric")
    %mean and stdtdev for strategies for all repetitions
    disp("mean strategy")
    res1(:,1)
    disp("stddev")
    res1(:,2)
    %mean and stddev for infection probabilites for all repetitions
    disp("mean infection load")
    res1(:,3)
    disp("stddev")
    res1(:,4)
    %When we have appropriate adjusted risk measurments they will go here

    %display results for the symmetric case 
    disp("asymmetric")
    %mean and stdtdev for strategies for all repetitions
    disp("mean strategy")
    res2(:,1)
    disp("stddev")
    res2(:,2)
    %mean and stddev for infection probabilites for all repetitions
    disp("mean infection load")
    res2(:,3)
    disp("stddev")
    res2(:,4)
    %When we have appropriate adjusted risk measurments they will go here
end

function ret = NumericalExperiment2()

    %environmental introduction rate
    el=0.2;
    %initial thresholds
    %taus = 0.09*ones(N,1);%...................................Constant
    taus = [0.2,0.2,0.1,0.1,0.1,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05]';%............Increasing

    %Subsection4.2
    %select adjacency matrices for a competetive, mildly monopolistic and
    %completely monopolistic trade network and find nash equilibria through
    %best response simulation
    [W3,N3,r3]=AdjMatSelect(5,3,1);
    [W4,N4,r4]=AdjMatSelect(6,3,1);
    [W5,N5,r5]=AdjMatSelect(7,3,1);

    res3=simulate(N3,W3,taus,el,100,100);
    res4=simulate(N4,W4,taus,el,100,100);
    res5=simulate(N5,W5,taus,el,100,100);

    %Display the mean and stddev of strategy, mean and stddev of infection
    %probability and both measures of risk for each network. 
    disp("competative")
    %mean and stdtdev for strategies for all repetitions
    disp("mean strategy")
    res3(:,1)
    disp("stddev")
    res3(:,2)
    %mean and stddev for infection probabilites for all repetitions
    disp("mean infection load")
    res3(:,3)
    disp("stddev")
    res3(:,4)

    disp("mild monopoly")
    %mean and stdtdev for strategies for all repetitions
    disp("mean strategy")
    res4(:,1)
    disp("stddev")
    res4(:,2)
    %mean and stddev for infection probabilites for all repetitions
    disp("mean infection load")
    res4(:,3)
    disp("stddev")
    res4(:,4)
    
    disp("Total Monopoly")
    %mean and stdtdev for strategies for all repetitions
    disp("mean strategy")
    res5(:,1)
    disp("stddev")
    res5(:,2)
    %mean and stddev for infection probabilites for all repetitions
    disp("mean infection load")
    res5(:,3)
    disp("stddev")
    res5(:,4)
end

function ret = example()
%%Example is a function which demonstrates what the other basic functions
%%in this script can do

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
    MBRsequence(x,N,W, taus,el,10)
    disp('Mean and standard deviation with repetition')
    simulate(N,W,taus,el,10,100)
    
end

function ret = simulate(N,W,taus, el, iter, reps)
%Simulate finds the Nash equilibrium by the MBR process for many initial
%data and records the variance in the case that there is bistability
%   N:      Number of Players
%   W:      Adjacency Matrix
%   taus:   Threshold Values 
%   el:     environmental introduction probability
%   iter:   Iteration limit on MBR
%   reps:   Number of repetitions with random intitial data

    
    res = zeros(N,reps);
    pres = zeros(N,reps);
    for i=1:reps
        x=rand(N,1);
        sol = MBRsequence(x,N,W,taus,el,iter);
        res(:,i)= sol;
        pres(:,1)= pstar(sol,N,W,el);
    end
    ret=[mean(res,2),std(res,0,2),mean(pres,2),std(pres,0,2)];
end

function eq = MBRsequence(x,N,W,taus,el, iter)
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
    i=1;
    finished = false;
    while ~finished
        oldx = x;
        x = MBR(oldx,N,W,taus,el);
        finished= or(~any(abs(x-oldx)>10^-5), i>=iter);
        i=1+1;
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