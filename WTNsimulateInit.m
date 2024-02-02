%This script seeds the immutable parameters for a simulation

rng("default")

%Number of Breeders
NB=10;

%Number of Distributors
ND = 3;

%Number of Retailers
NR = 6;

%Number of Hobbyists 
NO = 20;

%Random Infectivity
eps = 0.01;

%Costs of sigma
CB = 1;
CD = 1;
CR = 1;
CO = 1;

%Governement Regulation
GB=1;
GD=1;
GR=1;

%Personal Beliefs 
pB = rand(NB,1,"double");
pD = rand(ND,1,"double");
pR = rand(NR,1,"double");
pO = rand(NO,1,"double");

%Propensity to spend 
MD = rand(ND,1,"double");
MR = rand(NR,1,"double");
MO = rand(NO,1,"double");

%Owner Utility 
UO = rand(NO,1,"double");

%Owner Utility Loss
JO = rand(NO,1,"double");

%Simplified Coefficient;
KB = (-CB+GB+pB);
KD = (-CD+GD+pD)./MD;
KR = (-CR+GR+pR)./MR;
KO = (-CO +pO)./MO;

%InfectionState;
IB = zeros(NB,1);
ID = zeros(ND,1);
IR = zeros(NR,1);
IO = zeros(NO,1);

%Seed Adjacency matrix 
ABD = zeros(NB,ND);
for i = 1:1:ND
    ABD(:,i)= rand(NB,1,"double");
    ABD(:,i)=ABD(:,i)/sum(ABD(:,i));
end

ADR = zeros(ND,NR);
for i = 1:1:NR
    ADR(:,i)= rand(ND,1,"double");
    ADR(:,i)=ADR(:,i)/sum(ADR(:,i));
end

ARO = zeros(NR,NO);
for i = 1:1:NO
    ARO(:,i)= rand(NR,1,"double");
    ARO(:,i)=ARO(:,i)/sum(ARO(:,i));
end

