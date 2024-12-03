%%This script runs a numerical example for the downstream contagion game 

%max time steps
T=10;

%risk of contamination by environment
el=0.1;
%resolution with which to search for maxima in the strategy space
res = 1000;

%subsection4.1
%[W1,N1, r1]=AdjMatSelect(1,1,2);
%[W2,N2, r2]=AdjMatSelect(2,1,2);
% [xbar1,xdev1,xres1,pbar1,pdev1]=simulate(T,el,res,W1,N1,r1,10);
% [xbar2,xdev2,xres2,pbar2,pdev2]=simulate(T,el,res,W2,N2,r2,10);
% 
% 
% disp("symmetric")
% xbar1
% xdev1
% pbar1
% pdev1
% disp('Naive risk')
% 1-prod(1-pbar1)
% disp('weighted risk')
% weights1=weights(W1,xbar1,N1)
% weights1'*pbar1 
% disp("asymmetric")
% xbar2
% xdev2
% pbar2
% pdev2
% disp('Naive risk')
% 1-prod(1-pbar2)
% disp('weighted risk')
% weights2=weights(W2,xbar2,N2)
% weights2'*pbar2 

%Subsection4.2
% [W3,N3,r3]=AdjMatSelect(5,3,1);
% [W4,N4,r4]=AdjMatSelect(6,3,1);
% [W5,N5,r5]=AdjMatSelect(7,3,1);
% 
% [xbar3,xdev3,xres3,pbar3,pdev3]=simulate(T,el,res,W3,N3,r3,10);
% [xbar4,xdev4,xres4,pbar4,pdev4]=simulate(T,el,res,W4,N4,r4,10);
% [xbar5,xdev5,xres5,pbar5,pdev5]=simulate(T,el,res,W5,N5,r5,10);
% 
% 
% disp("competative")
% xbar3
% xdev3
% pbar3
% pdev3
% disp('Naive risk')
% 1-prod(1-pbar3)
% disp('weighted risk')
% weights3=weights(W3,xbar3,N3)
% weights3'*pbar3 
% disp("mild monopoly")
% xbar4
% xdev4
% pbar4
% pdev4
% disp('Naive risk')
% 1-prod(1-pbar4)
% disp('weighted risk')
% weights4=weights(W4,xbar4,N4)
% weights4'*pbar4 
% disp("total monopoly")
% xbar5
% xdev5
% pbar5
% pdev5
% disp('Naive risk')
% 1-prod(1-pbar5)
% disp('weighted risk')
% weights5=weights(W5,xbar5,N5)
% weights5'*pbar5 


%robust network
% figure()
% [W1,N1,r1]=AdjMatSelect(1,2,1);
% [sol, band]=parameterscan(T,el,res,W1,N1,100,0,1.5,80,r1,[5,6,7,8]);
% plot(sol(:,1),sol(:,7),sol(:,1),sol(:,3),sol(:,1),sol(:,4),sol(:,1),sol(:,6))
% title('Parameter Cascades in Symmetric Network')
% ylim([0,0.55])
% xlabel('Intrinsic Benefit Weight')
% ylabel('Equilibrium Strategy')
% legend('Consumer 1','Producer','Distributor','Consumer 2')


%fragile network
% figure(r1)
% [W1,N1,r1]=AdjMatSelect(2,2,1);
% [sol, band]=parameterscan(T,el,res,W1,N1,100,0,1.5,80,r1,[5,6,7,8]);
% plot(sol(:,1),sol(:,7),sol(:,1),sol(:,3),sol(:,1),sol(:,4),sol(:,1),sol(:,6))
% title('Parameter Cascades in Asymmetric Network')
% ylim([0,0.55])
% xlabel('Intrinsic Benefit Weight')
% ylabel('Equilibrium Strategy')
% legend('Consumer 1','Producer','Distributor','Consumer 2')

%defector logic
[W1,N1, r1]=AdjMatSelect(1,1,2);
[xbar1,xdev1,xres1,pbar1,pdev1]=simulate(T,el,res,W1,N1,r1,10);
[xbar2,xdev2,xres2,pbar2,pdev2]=defector(T,el,res,W1,N1,r1,10,2,0);
[xbar3,xdev3,xres3,pbar3,pdev3]=defector(T,el,res,W1,N1,r1,10,2,1);
disp("symmetric")
xbar1
xdev1
pbar1
pdev1
disp('Naive risk')
1-prod(1-pbar1)
disp('weighted risk')
weights1=weights(W1,xbar1,N1)
weights1'*pbar1 
disp("Consumer defects to 0")
xbar2
xdev2
pbar2
pdev2
disp('Naive risk')
1-prod(1-pbar2)
disp('weighted risk')
weights2=weights(W1,xbar2,N1)
weights2'*pbar2 

disp("Consumer defects to 1")
xbar3
xdev3
pbar3
pdev3
disp('Naive risk')
1-prod(1-pbar3)
disp('weighted risk')
weights3=weights(W1,xbar3,N1)
weights3'*pbar3 


function [curve, band] = parameterscan(T,el,res, W,N,reps, minr, maxr, rres, rdefault, rfoci)
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


%payoff function for pi 
function [w] = pii(focus, xfocus,x,N,W,el,rfocus)
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

%Stationary probability of infection
function p = pstar(x,N,W,el)
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
    z = zeros(N,1);
    %Mat = eye(N);
    Mat=(repmat((1-X).^2,1,size(W,2))'.*W);
    while nnz(Mat)>0
        z=z+sum(Mat)';
        Mat = Mat*(repmat((1-X).^2,1,size(W,2))'.*W);
    end
end

function Rg = weightedrisk(W,X,N,pstar)
    w =weights(W,X,N);
    Rg=w'*pstar;
end

function [xbar, xdev, xres, pbar, pdev, pres] =defector(T, el, res, W, N, r, reps, defi, defx)
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


