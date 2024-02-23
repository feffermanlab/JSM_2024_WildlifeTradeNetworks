%This script takes a populations seeded by WTNSimulateinit
%and computes the payoff assuming a long run infections state. 

WTNsimulateInit;


%sigma
sigmaB = rand(NB,1,"double");
sigmaD = rand(ND,1,"double");
sigmaR = rand(NR,1,"double");
sigmaO = rand(NO,1,"double");

%l
lB = rand(NB,1,'double');
lD = rand(ND,1,'double');
lR = rand(NR,1,'double');
lO = rand(NO,1,'double');

%alpha
alphas = rand(NB,1,'double');
alphabar = mean(alphas);

%beta
betas = rand(ND,1,'double');
betabar = mean(betas);

%gamma
gammas = rand(NR,1,'double');
gammabar = mean(gammas);


%compute equilibrium infection rates

%construct the matrix $B$ and $D$

%solve the nonlinear problem with fsolve


%compute payoffs

%Transactions
   WBD = ABD.*(1-lD*(IB.')).';
   WDR = ADR.*(1-lR*(ID.')).';
   WRO = ARO.*(1-lO*(IR.')).';
   
   piB = piB +(KB.*sigmaB) + alphabar.*(WBD*MD);
   piD = piD +((KD./MD).*sigmaD/alphabar) + (betabar/alphabar).*(WDR*MR./MD) - (WBD.'*ones(NB,1));
   piR = piR +((KR./MR).*sigmaR/betabar) + (gammabar/betabar).*(WRO*MO./MR) - (WDR.'*ones(ND,1));
   piO = piO +((KO./MO).*sigmaO/gammabar) - (WRO.'*ones(NR,1))+ UO -JO.*IO;



