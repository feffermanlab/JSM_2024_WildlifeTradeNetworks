

B = [0,0,0,0;0.01,0,0,0;0.01,0.8,0,0;0.01,0,0.8,0];
D = diag([1,0.9,0.9,0.9]);

Bminor = [0,0,0;0.8,0,0;0,0.8,0];
Dminor = diag([0.9,0.9,0.9]);
eps=0.01;
res = zeros(1000,6);

%T([0.5;0.5;0.5],Bminor,Dminor,eps)
%T2([0.5;0.5;0.5],B,D)

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

res

F= @(p)T(p,Bminor,Dminor,eps)-p;
G= @(p)T2(p,B,D)-p;
H =@(p)T3(p,Bminor,Dminor,eps)-p;

%F([1;0;0])
%G([1;0;0])

[p,fval]=fsolve(F,[0.5;0.5;0.5])
[p,fval]=fsolve(G,[0.5;0.5;0.5])
[p,fval]=fsolve(H,[0.5;0.5;0.5])

p = zeros(4,1);
p(1)=1;
for i = 2:4
    s = B(i,:)*p;
    p(i)=s/(1-D(i,i)+s);
end

p


function pplus = T(p,B,D,eps)
    pplus =B*p+D*p-(B*p+eps).*p+eps;
end
    

function pplus = T2(p,B,D)
    l =[1;p];
    pplus = (B+D)*l-((B*l).*l);
    pplus = pplus(2:end,:);
end

function pplus = T3(p,B,D,eps)
    pplus = (B*p+eps).*(1-p)+D*p;
end



