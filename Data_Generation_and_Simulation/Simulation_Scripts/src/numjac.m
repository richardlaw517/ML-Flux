function J = numjac(simulate,net,xch,input,mea,fmea,var,isMMID,A,DeMa)
%% calculate numerical jacobian matrix
%persistent F net_l xch_l RL 20221102

% RL 20221102
% if isempty(net_l) || net_l~=length(net) || xch_l~=length(xch)
%     net_l=length(net);
%     xch_l=length(xch);
%     if isMMID
%         F=@(x)varssr_mmid(simulate,x(1:net_l),x(net_l+1:end),input,mea,fmea,var,model,A,DeMa);
%     else
%         F=@(x)varssr(simulate,x(1:net_l),x(net_l+1:end),input,mea,fmea,var);
%     end
% 
% end

net_l=length(net);
xch_l=length(xch);
if isMMID
    F=@(x)varssr_mmid(simulate,x(1:net_l),x(net_l+1:end),input,mea,fmea,var,model,A,DeMa);
else
    F=@(x)varssr(simulate,x(1:net_l),x(net_l+1:end),input,mea,fmea,var);
end
[~,~,~,~,~,~,J] = lsqnonlin(F,[net; xch],[],[],optimset('MaxFunEval',0,'display','off')); % lsqnonlin from optimization toolbox
%Jalt=jacobianest(F,[net; xch]); % jacobianest from DERVESTsuite toolbox (https://www.mathworks.com/matlabcentral/fileexchange/13490-adaptive-robust-numerical-differentiation)
%J=jacobianest(F,[net; xch]); % slow/read the license info under MATLAB/DERIVESTsuite
J=full(J);
%Jalt = full(Jalt);
end

