%%This script runs a numerical example for the downstream contagion game 

%max time steps
T=10;

%risk of contamination by environment
el=0.1;
%resolution with which to search for maxima in the strategy space
res = 1000;

%subsection4.1
%get a symmetric and ansymmetric adjacency matrix and find equilibrium
%strategies through best response simulation
[W1,N1, r1]=AdjMatSelect(1,1,2);
[W2,N2, r2]=AdjMatSelect(2,1,2);
[xbar1,xdev1,xres1,pbar1,pdev1]=simulate(T,el,res,W1,N1,r1,10);
[xbar2,xdev2,xres2,pbar2,pdev2]=simulate(T,el,res,W2,N2,r2,10);
 
%display results for the symmetric case 
disp("symmetric")
%mean and stdtdev for strategies for all repetitions
xbar1
xdev1
%mean and stddev for infection probabilites for all repetitions
pbar1
pdev1
%Naive risk
disp('Naive risk')
1-prod(1-pbar1)
%weighted risk
disp('weighted risk')
weights1=weights(W1,xbar1,N1);
weights1'*pbar1 

%display results for the symmetric case
disp("asymmetric")
%mean and stdtdev for strategies for all repetitions
xbar2
xdev2
%mean and stddev for infection probabilites for all repetitions
pbar2
pdev2
%naive risk
disp('Naive risk')
1-prod(1-pbar2)
%weighted risk
disp('weighted risk')
weights2=weights(W2,xbar2,N2);
weights2'*pbar2 

%Subsection4.2
%select adjacency matrices for a competetive, mildly monopolistic and
%completely monopolistic trade network and find nash equilibria through
%best response simulation
[W3,N3,r3]=AdjMatSelect(5,3,1);
[W4,N4,r4]=AdjMatSelect(6,3,1);
[W5,N5,r5]=AdjMatSelect(7,3,1);

[xbar3,xdev3,xres3,pbar3,pdev3]=simulate(T,el,res,W3,N3,r3,10);
[xbar4,xdev4,xres4,pbar4,pdev4]=simulate(T,el,res,W4,N4,r4,10);
[xbar5,xdev5,xres5,pbar5,pdev5]=simulate(T,el,res,W5,N5,r5,10);

%Display the mean and stddev of strategy, mean and stddev of infection
%probability and both measures of risk for each network. 
disp("competative")
xbar3
xdev3
pbar3
pdev3
disp('Naive risk')
1-prod(1-pbar3)
disp('weighted risk')
weights3=weights(W3,xbar3,N3)
weights3'*pbar3 
disp("mild monopoly")
xbar4
xdev4
pbar4
pdev4
disp('Naive risk')
1-prod(1-pbar4)
disp('weighted risk')
weights4=weights(W4,xbar4,N4)
weights4'*pbar4 
disp("total monopoly")
xbar5
xdev5
pbar5
pdev5
disp('Naive risk')
1-prod(1-pbar5)
disp('weighted risk')
weights5=weights(W5,xbar5,N5)
weights5'*pbar5 

%section 4.3
% for a symmetric and asymmetric network, change the weight of intrinsic
% benefits for owners and observe visually the change in equlibrium
% strategy
%line weight for the graph so I didn't have to change it in 7 places
lw=1.5;
%symmetric network
figure()
%select symmetric adjacency matrix and simulate on a strip in the parameter
%space
[W1,N1,r1]=AdjMatSelect(1,2,1);
[sol, band]=parameterscan(T,el,res,W1,N1,100,0,1.5,80,r1,[5,6,7,8]);
%plot the results
h=plot(sol(:,1),sol(:,7),"-",sol(:,1),sol(:,3),"--",sol(:,1),sol(:,4),":");
set(h(3),'Color',"#0000a4","LineWidth",lw)
set(h(2),'Color','#bc272d',"LineWidth",lw)
set(h(1),'Color','#e9c716',"LineWidth",lw)
title('Parameter Cascades in Symmetric Network')
ylim([0,0.55])
xlabel('Intrinsic Benefit Weight')
ylabel('Equilibrium Strategy')
legend('Consumer','Producer','Distributor','location','southeast')


%asymmetric network network
figure()
%select asymmetric adjacency matrix and simulate on a strip in the
%parameter space
[W1,N1,r1]=AdjMatSelect(2,2,1);
[sol, band]=parameterscan(T,el,res,W1,N1,100,0,1.5,80,r1,[5,6,7,8]);
%plot the results
h=plot(sol(:,1),sol(:,7),"-", sol(:,1),sol(:,3),"--",sol(:,1),sol(:,4),":",sol(:,1),sol(:,5),"-.");
set(h(3),'Color',"#0000a4","LineWidth",lw)
set(h(2),'Color','#bc272d',"LineWidth",lw)
set(h(1),'Color','#e9c716',"LineWidth",lw)
set(h(4),'Color','#50ad9f',"LineWidth",lw)
title('Parameter Cascades in Asymmetric Network')
ylim([0,0.55])
xlabel('Intrinsic Benefit Weight')
ylabel('Equilibrium Strategy')
legend('Consumer','Producer','Distributor (a1)','Distributor (a2)', 'location','southeast')

%subsection 4.4
%Select the symmetric adjacency matrix and find Nash equilibria through myopic
%best response in the case that everyone is rational, in the case that a
%selected player defects to 0 and in the case the selected player defects
%to 1
[W1,N1, r1]=AdjMatSelect(1,1,2);
[xbar1,xdev1,xres1,pbar1,pdev1]=simulate(T,el,res,W1,N1,r1,10);
[xbar2,xdev2,xres2,pbar2,pdev2]=defector(T,el,res,W1,N1,r1,10,2,0);
[xbar3,xdev3,xres3,pbar3,pdev3]=defector(T,el,res,W1,N1,r1,10,2,1);

%display the same statistics as in subsection 4.1 and 4.2 for each of the
%cases simulate3d above
disp("Rational")
xbar1
xdev1
pbar1
pdev1
disp('Naive risk')
1-prod(1-pbar1)
disp('weighted risk')
weights1=weights(W1,xbar1,N1);
weights1'*pbar1 

disp("Consumer defects to 0")
xbar2
xdev2
pbar2
pdev2
disp('Naive risk')
1-prod(1-pbar2)
disp('weighted risk')
weights2=weights(W1,xbar2,N1);
weights2'*pbar2 

disp("Consumer defects to 1")
xbar3
xdev3
pbar3
pdev3
disp('Naive risk')
1-prod(1-pbar3)
disp('weighted risk')
weights3=weights(W1,xbar3,N1);
weights3'*pbar3 


function [curve, band] = parameterscan(T,el,res, W,N,reps, minr, maxr, rres, rdefault, rfoci)
%PARAMETERSCAN finds Nash equlibria through myopic best response simulation
%in a band of the parameter space
%   T: Max time
%   el: risk of contamination from the environment
%   res: number of time steps in a single time unit
%   W: adjacency matrix
%   N: size of network
%   reps: number of repetitions 
%   minr: minimum of the parameter to be changed
%   maxr: maximum of the parameter to be changed
%   rres: number of parameter test values in parameter interval
%   rdefault: parametrer value for players not being changed
%   rfoci: players being changed

%This function tests a range of different r values for a subset of players
%(rfoci) for each test point in the range of r values, the function finds a
%Nash equilibrium through simulatino of myopic best response and saves the
%average equilibrium strategy and infection probability as well as the weighted risk
%in curve and saves the stddevs in band 
    curve = zeros(rres, 2*N+2);
    band = zeros(rres,2*N+1);
    rspace = linspace(minr,maxr,rres);
    for i = 1:rres
        disp(i)
        rs = rdefault;
        rs(rfoci)=rspace(i);
        [xbar,xdev,xres,pbar,pdev,pres] = simulate(T,el,res,W,N,rs,reps);
        curve(i,1)=rspace(i);
        band(i,1)=rspace(i);
        curve(i,2:N+1)=xbar;
        band(i,2:N+1)=xdev;
        curve(i,N+2:2*N+1)=pbar;
        band(i,N+2:2*N+1)=pdev;
        curve(i,2*N+2)=weightedrisk(W,xbar,N,pbar);
    end
end



function [xbar, xdev, xres, pbar, pdev, pres] =simulate(T, el, res, W, N, r, reps)
%SIMULATE finds Nash equilibrium through simulated myopic best response
%   T: Max time
%   el: risk of contamination from the environment
%   res: number of time steps in a single time unit
%   W: adjacency matrix
%   N: size of network
%   reps: number of repetitions 

%For a single point in the parameter space and for a single network, this
%function finds the Nash equilibrium through repeated myopic best response.
%It is not rigorously proven that this process will always converge to the
%same equilibrium (i.e. it is not clear that this process in monostable) so
%we repeat the process many times (reps) and measure the mean and stddev

xbar = zeros(N,1);
xdev = zeros(N,1);
pbar = zeros(N,1);
pdev = zeros(N,1);
xres = zeros(N,reps);
pres = zeros(N,reps);
for rep = 1:reps    
    %initialize strategies 
    x= rand(N,1);
    sol =zeros(N,T+1);
    sol(:,1)=x;

    %start simulation
    complete = false;
    t=1;
    while(t<T && ~complete) 
        x=sol(:,t);
    for i =1:N
        %for focal individual i, calculate payoff everywhere in strategy space
        f = @(xi)pii(i,xi,x,N,W,el,r(i));
        strats = linspace(0,1,res);
        y= zeros(res,1);
        for j =1:res
            y(j)=f(strats(j));
        end
        %find strategy which maximizes individual payoff
        [M,I]=max(y);
        %add this to solution 
        sol(i,t+1)=strats(I);
    end
    %check if the strategy has changed
    complete = (~all(sol(:,t+1)-sol(:,t)));
    t= t+1;
    end
    %print out final time, final strategy profile, and final infection
    %probability
    fin = sol(:,t);
    p=pstar(fin,N,W,el);
    xres(:,rep)=fin;
    pres(:,rep)=p;
end
xbar = mean(xres,2);
xdev = std(xres,0,2);
pbar = mean(pres,2);
pdev = std(pres,0,2);
end


function [w] = pii(focus, xfocus,x,N,W,el,rfocus)
%PII is the payoff function for a player using a prescribed strategy
%   focus: the focal player
%   xfocus: the strategy being played by the focal player
%   x: the strategy being played by all other players
%   N: the number of players
%   W: Adjacency matrix for the network
%   el: the propbability of infection from the environment
%   rfocus: Weight of intrinsic benefit for the focal player

%This function depends heavily on the form of the toy model. With different
%functional forms this entire function would need to be rewritten. 

    %compute p*
    p = zeros(N,1);
    for i = 1:N
        d1 = el*(1-xfocus);
        d2 = (1-xfocus)^2.*W(i,:)*p;
        d=d1+d2;
        p(i)= d/(xfocus+d);
    end 
    %individual payoff 
    d=rfocus*xfocus*(1-xfocus)*(1-p(focus));
    
    %downstream interactions
    F=(1-xfocus).*W(:,focus).*(1-p(focus)*x);
    b=sum(F);
   
    %upstream interactions
    G = (x'-1).*W(focus,:).*(1-xfocus*p');
    c = sum(G);

    w= d+b+c;
end


function p = pstar(x,N,W,el)
%PSTAR computes the stationary probability of infection for a given network
%and strategy profile
%   x: a strategy profile
%   N: number of players
%   W: ASdjacency matrix
%   el: Probability of infection from the environment

%This function depends heavilty on the functional forms of the toy model.
%The same function would not work with different functional forms. We are
%using the results outlined in Theorem 1 of the manuscript. 
    p = zeros(N,1);
    for i = 1:N
        d1 = el*(1-x(i));
        d2 = (1-x(i))^2*W(i,:)*p;
        %d2 = sum(W(i,:)*(1-x(i))^2.*p);
        d=d1+d2;
        p(i)= d/(x(i)+d);
    end 
end

function z = weights(W,X,N)
%WEIGHTS computes the weights needed form computing weighted risk
%   W: adjacency matrix
%   X: Strategy profile
%   N: Number of players

%This function computes the number of paths of infection that each
%individual is accountable for in a network
    z = zeros(N,1);
    %Mat = eye(N);
    Mat=(repmat((1-X).^2,1,size(W,2))'.*W);
    while nnz(Mat)>0
        z=z+sum(Mat)';
        Mat = Mat*(repmat((1-X).^2,1,size(W,2))'.*W);
    end
end

function Rg = weightedrisk(W,X,N,pstar)
%WEIGHTEDRISK computes the weighted risk of a strategy profile on a network
%   W: adjacency matrix
%   X: stategy profile
%   N: Number of players
%   pstar: the stationary probability of infection

%This function computes weighted risk by multiplying the number of routs of
%infection each individual is potent to start by the likelihood that they
%may be infected themselves. 
    w =weights(W,X,N);
    Rg=w'*pstar;
end

function [xbar, xdev, xres, pbar, pdev, pres] =defector(T, el, res, W, N, r, reps, defi, defx)
%DEFECTOR Finds Nash equilibria through simulated myopic best response in
%the case that a player is not playing rationally and selects a single
%strategy
%   T: Max time
%   el: risk of contamination from the environment
%   res: number of time steps in a single time unit
%   W: adjacency matrix
%   N: size of network
%   reps: number of repetitions 
%   defi: player who will defect
%   defx: strategy to which they will defect

%This function computes through myopic best response the Nash equilibrium
%of the game. However, one individual (defi) does not update their strategy
%through MBR but rather always plays strategy defx
xbar = zeros(N,1);
xdev = zeros(N,1);
pbar = zeros(N,1);
pdev = zeros(N,1);
xres = zeros(N,reps);
pres = zeros(N,reps);
for rep = 1:reps    
    %initialize strategies 
    x= rand(N,1);
    x(defi)=defx;
    sol =zeros(N,T+1);
    sol(:,1)=x;

    %start simulation
    complete = false;
    t=1;
    while(t<T && ~complete) 
        x=sol(:,t);
    for i =setdiff(1:N,defi)
        %for focal individual i, calculate payoff everywhere in strategy space
        f = @(xi)pii(i,xi,x,N,W,el,r(i));
        strats = linspace(0,1,res);
        y= zeros(res,1);
        for j =1:res
            y(j)=f(strats(j));
        end
        %find strategy which maximizes individual payoff
        [M,I]=max(y);
        %add this to solution 
        sol(i,t+1)=strats(I);
    end
    sol(defi,t+1)=defx;
    %check if the strategy has changed
    complete = (~all(sol(:,t+1)-sol(:,t)));
    t= t+1;
    end
    %print out final time, final strategy profile, and final infection
    %probability
    fin = sol(:,t);
    p=pstar(fin,N,W,el);
    xres(:,rep)=fin;
    pres(:,rep)=p;
end
xbar = mean(xres,2);
xdev = std(xres,0,2);
pbar = mean(pres,2);
pdev = std(pres,0,2);

end


