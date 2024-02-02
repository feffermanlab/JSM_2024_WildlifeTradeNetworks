%This Script takes a population seeded by WTNSimulateInit
%and runs it for T time steps and calculates final fitness

WTNsimulateInit;

T=100000;

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

%Total payoffs
piB = zeros(NB,1,'double');
piD = zeros(ND,1,'double');
piR = zeros(NR,1,'double');
piO = zeros(NO,1,'double');
 

for t = 1:T
    %Transactions
    WBD = ABD.*(1-lD*(IB.')).';
    WDR = ADR.*(1-lR*(ID.')).';
    WRO = ARO.*(1-lO*(IR.')).';
    
    piB = piB +(KB.*sigmaB) + (WBD*MD);
    piD = piD +((KD./MD).*sigmaD) + (WDR*MR./MD) - (WBD.'*ones(NB,1));
    piR = piR +((KR./MR).*sigmaR) + (WRO*MO./MR) - (WDR.'*ones(ND,1));
    piO = piO +((KO./MO).*sigmaO) - (WRO.'*ones(NR,1))+ UO -JO.*IO;

    %Transactional Infection 
    % Each Stakeholder has a porbability of infection based on upstream
    % interaction
    InfectD =((WBD.'*IB)-sigmaD>0);
    InfectR =((WDR.'*ID)-sigmaR>0);
    InfectO =((WRO.'*IR)-sigmaO>0);

    ID= or(ID, InfectD);
    IR= or(IR, InfectR);
    IO= or(IO, InfectO);

    %Random Infection
    %Each stakeholder has a probability of infection based on sigma
    RInfectB = (rand(NB,1)<eps);
    RInfectD = (rand(ND,1)<eps);
    RInfectR = (rand(NR,1)<eps);
    RInfectO = (rand(NO,1)<eps);
    
    IB= or(IB, RInfectB);
    ID= or(ID, RInfectD);
    IR= or(IR, RInfectR);
    IO= or(IO, RInfectO);

    %Recovery
    RecoverB = (rand(NB,1)-0.1.*sigmaB<0);
    RecoverD = (rand(ND,1)-0.1.*sigmaD<0);
    RecoverR = (rand(NR,1)-0.1.*sigmaR<0);
    RecoverO = (rand(NO,1)-0.1.*sigmaO<0);

    IB = and(IB, not(RecoverB));
    ID = and(ID, not(RecoverD));
    IR = and(IR, not(RecoverR));
    IO = and(IO, not(RecoverO));
        
end

piB
piD
piR
piO

IB
ID
IR
IO

