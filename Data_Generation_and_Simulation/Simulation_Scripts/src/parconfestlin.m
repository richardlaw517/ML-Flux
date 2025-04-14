function [sm,Hess,inverH,varFlux] = parconfestlin(simulate,net_opt,xch_opt,input,mea,fmea,covar,model,isMMID,A,DeMa)
% Calculate sensitivity matrix (Jacobian)
sm = numjac(simulate,net_opt,xch_opt,input,mea,fmea,covar,isMMID,A,DeMa);
% Calculate jacobian (jacobian vector) and hessian
inverCovar = pinv(covar);
Hess = sm'*inverCovar*sm;
inverH = pinv(Hess);
nullbase = blkdiag(model.kernel_net,model.kernel_xch);
varFlux = nullbase*inverH*nullbase';
end