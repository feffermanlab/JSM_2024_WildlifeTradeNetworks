%This script was used as a numerical investigation of the monostability of
%of the propbability of infection difference equation. The rigorous result
%is found as Theorem 1 in the manuscript

%Example weighted adjacency matrix
B = [0,0,0,0;0.01,0,0,0;0.01,0.8,0,0;0.01,0,0.8,0];
%example recovery complements
D = diag([1,0.9,0.9,0.9]);
%Minor of B
Bminor = [0,0,0;0.8,0,0;0,0.8,0];
%Minor of D
Dminor = diag([0.9,0.9,0.9]);
%epsilon
eps=0.01;
%variable for the results
res = zeros(1000,6);

%solve the difference equation and find where it reaches an apparent
%equilibrium
for i = 0:9
    for j = 0:9
        for k = 0:9
            pnew = [i/10;j/10;k/10];
            p = [10;10;10];
            while norm(p-pnew,2)>10^(-5)
                p=pnew;
                pnew = T(p,Bminor,Dminor,eps);
            end
            res(i*100+j*10+k+1,:)=[i,j,k,pnew'];
        end
    end
    
end
%print results
res

%Functions for each formulation of T
F= @(p)T(p,Bminor,Dminor,eps)-p;
G= @(p)T2(p,B,D)-p;
H =@(p)T3(p,Bminor,Dminor,eps)-p;

%F([1;0;0])
%G([1;0;0])

%Proof of concept that each of the formulations for T are the same
[p,fval]=fsolve(F,[0.5;0.5;0.5])
[p,fval]=fsolve(G,[0.5;0.5;0.5])
[p,fval]=fsolve(H,[0.5;0.5;0.5])

%Proof of concept that the fixed point is directly computable
p = zeros(4,1);
p(1)=1;
for i = 2:4
    s = B(i,:)*p;
    p(i)=s/(1-D(i,i)+s);
end

p
%We suspect that each of the results of the above for loop as well as the
%resultf of T T2 and T3 should be equivalent


function pplus = T(p,B,D,eps)
%T next generation operator for probability of infection
    pplus =B*p+D*p-(B*p+eps).*p+eps;
end
    

function pplus = T2(p,B,D)
%T2 another formulation of the next genereation operator for probability of
%   infection
    l =[1;p];
    pplus = (B+D)*l-((B*l).*l);
    pplus = pplus(2:end,:);
end

function pplus = T3(p,B,D,eps)
%T3 yet another formulation of the next generation operator for probability
%   of infection
    pplus = (B*p+eps).*(1-p)+D*p;
end



