function [free_net,free_xch,all_free] = fluxinit(model,ineq,eq,n,ref,seed)
%% initialize random fluxes that fall within constraints
% input: model and ineq that specify independent (free) metabolic reactions
% ; inequality and equality constratins; n, the number of flux 
% distributions; reference reaction to which all fluxex are normalized; and
% seed (default 0).
% output: randomized set of free net and exchange reactions that fit
% constraints

if nargin<5 || isempty(ref)
    ref=1;
    seed=0;
elseif nargin<6 || isempty(seed)
    seed=0;
end

% Aeq=[model.kernel_net(ref,:) model.kernel_xch(ref,:)];
% beq=1;
rng(seed);

net_l=size(model.kernel_net,2);
xch_l=size(model.kernel_xch,2);
net_min=zeros(net_l,1);
xch_min=zeros(xch_l,1);
all_min=[net_min; xch_min];
all_max=all_min;
options=optimoptions('linprog','Display','off');

for i=1:length(all_min)
    f=zeros(1,length(all_min));
    f(i)=1;
    x=linprog(f,ineq.A,ineq.b,eq.A,eq.b,[],[],options);
    if isempty(x)
        all_min(i)=-10;
    else
        all_min(i)=x(i);
    end
    f(i)=-1;
    x=linprog(f,ineq.A,ineq.b,eq.A,eq.b,[],[],options);
    if isempty(x)
        all_max(i)=10;
    else
        all_max(i)=x(i);
    end
end

free_net=zeros(net_l,n);
free_xch=zeros(xch_l,n);
m=1;
for i=1:20000000*n
    net=all_min(1:net_l)+(all_max(1:net_l)-all_min(1:net_l)).*rand(net_l,1);
    xch=all_max(net_l+1:end).*rand(xch_l,1);
    refflux=model.kernel_net(ref,:)*net;
    net=net/refflux;
    xch=xch/refflux;
    if ineq.A*[net;xch]<=ineq.b
        free_net(:,m)=net;
        free_xch(:,m)=xch;
        m=m+1;
    end
    if m>n
        break
    end
end

all_free=[free_net; free_xch];
remove_col=find(all(all_free==0));
free_net(:,remove_col)=[];
free_xch(:,remove_col)=[];
all_free(:,remove_col)=[];
if n==1 && isempty(free_net)
    net=all_min(1:net_l)+(all_max(1:net_l)-all_min(1:net_l)).*rand(net_l,1);
    xch=all_max(net_l+1:end).*rand(xch_l,1);
    refflux=model.kernel_net(ref,:)*net;
    net=net/refflux;
    xch=xch/refflux;
    free_net(:,n)=net;
    free_xch(:,n)=xch;
    warning('an infeasible initial flux distribution (seed %i)',seed)
end
end