function [net_opt,xch_opt,info] = fluxest2(simulate,model,net,xch,ineq,eq,input,mea,fmea,var,temp,iteration,prevOpt,CIter,keeplast,useParallel,parDataQueue)
%% gradient-descent method for the cumomer model to estimate optimal fluxes
% Dependencies: xlsreadineq.m, varssr.m
% input: simulate (function handle), model (MATLAB structure containing
% metabolic model), net and xch (free flux vectors), input (input isotope
% substrate labeling patterns), mea and var (average measured labeling
% patterns and variances), fmea (structure of flux measurements),
% iter (maximum number of iterations)
% output: net_opt and xch_opt (optimized free net and exchange fluxes),
% info (information regarding fit)

% Reshape net and xch matrices into column vectors
net=net(:);
xch=xch(:);

app1=MFEA_GUI_Hist_v3_;
net_l=length(net);
net_ll=net_l+1;
F=@(x)objfunc(simulate,x(1:net_l),x(net_ll:end),input,mea,fmea,var);

%% optimization steps
% When no inequality constraints are specified
if isempty(ineq)
    [ineq]=xlsreadineq([],model,net,xch);
end
% for flux estimation
fluxest=0;
if isempty(eq)
    eq.A=zeros(1,net_l+length(xch));
    eq.b=0;
    fluxest=1;
end
% initialize
init_score=F([net;xch]);
init_fea=[min(ineq.A*[net;xch]<=ineq.b) min(eq.A*[net;xch]==eq.b)];
disp(['optimization begins with score = ', num2str(init_score), '; ineq = ', num2str(init_fea(1)), '; eq = ', num2str(init_fea(2))])

if nargin<14 || isempty(prevOpt)
    prevOpt=init_score-1;
end
if nargin<15 || isempty(CIter)
    CIter=900;
end
if nargin<16 || isempty(keeplast)
    keeplast=0;
end

% if init_score is less than the opt score, shorten the optimization
if sum(init_fea)==2 && init_score<=prevOpt
    CIter=50;
    disp(['MaxIter = ', num2str(CIter)])
end

if fluxest
    % for flux estimation
    options=optimset('OutputFcn',@(x,optimValue,state) outfun(x,optimValue,state,parDataQueue),'Algorithm','interior-point' ...
        ,'TolX',1e-10,'TolCon',1e-30,'TolFun',1e-7 ...
        ,'MaxFunEvals',110000,'MaxIter',1100,'PlotFcn','optimplotfval');
else
    % for confidence interval estimation
    options=optimset('OutputFcn',@(x,optimValue,state) outfun(x,optimValue,state,parDataQueue),'Algorithm','interior-point' ...
        ,'Display','notify-detailed' ...
        ,'TolX',1e-10,'TolCon',1e-30,'TolFun',1e-7 ...
        ,'MaxFunEvals',100*CIter,'MaxIter',CIter,'PlotFcn','optimplotfval');
end

% MultiStart, ga, GlobalSearch
while isnan(init_score)
    net=-2.5+5*rand(size(net));
    xch=5*rand(size(xch));
    init_score=F([net;xch]);
    init_fea=[min(ineq.A*[net;xch]<=ineq.b) min(eq.A*[net;xch]==eq.b)];
    disp(['optimization begins with score = ', num2str(init_score), '; ineq = ', num2str(init_fea(1)), '; eq = ', num2str(init_fea(2))])
end
% options=optimset('Display','iter','MaxFunEvals',5000,'MaxIter',5000);
% [info.x,info.fval,info.exitflag,info.output]=fminsearch(F,[net;xch],options)

% Turn off unnecessary warnings
warning('off','MATLAB:singularMatrix')
warning('off','MATLAB:nearlySingularMatrix')
warning('off','MATLAB:illConditionedMatrix')
try
    [info.x,info.fval,info.exitflag,info.output]=fmincon(F,[net;xch],ineq.A,ineq.b,eq.A,eq.b,[],[],[],options);
catch
    warning('Error in fmincon, trying a new flux set with fewer iterations')
    % if there's error, try a new flux set with fewer iterations
    options=optimset('OutputFcn',@(x,optimValue,state) outfun(x,optimValue,state,parDataQueue),'Algorithm','interior-point' ...
        ,'MaxFunEvals',30000,'MaxIter',300,'PlotFcn','optimplotfval');
    net=-1.5+3*rand(size(net));
    xch=3*rand(size(xch));
    [info.x,info.fval,info.exitflag,info.output]=fmincon(F,[net;xch],ineq.A,ineq.b,eq.A,eq.b,[],[],[],options);
end

opt_fea=[min(ineq.A*info.x<=ineq.b) min(eq.A*info.x==eq.b)];
if sum(opt_fea)<2
    warning(['infeasible opt: init = ' num2str(init_score), '; opt = ' num2str(info.fval), '; ineq = ', num2str(opt_fea(1)), '; eq = ', num2str(opt_fea(2))])
    disp(['eq.A*[net;xch] = ', num2str(eq.A*info.x), '; eq.b = ', num2str(eq.b)])
end

if sum(init_fea)==2 && init_score<info.fval && keeplast==0
    disp(['better init: init = ' num2str(init_score), '; opt = ' num2str(info.fval), '; ineq = ', num2str(opt_fea(1)), '; eq = ', num2str(opt_fea(2))])
    net_opt=net;
    xch_opt=xch;
    info.fval=init_score;
else
    net_opt=info.x(1:net_l);
    xch_opt=info.x(net_ll:end);
end
previousScore = [];

app1.FinalVSSREditField.Value = info.fval;
h = figure('visible','off');
uit = app1.MeasuredMID_Table;
d = get(uit,'data');
mea_names = fieldnames(mea);
uit.ColumnName = mea_names;
k = 0;
for i = 1:numel(mea_names)
    a = length(getfield(mea,mea_names{i}));
    if k <= a
        k = a;
    end
end
%         Above for loop calculates maximum length of mea structure
%         fields.
for i = 1:numel(mea_names)
    l = getfield(mea,mea_names{i})';
    l(end+1:k) = 0;
    d(:,i) = l;
end
set(uit, 'data', d);
%         Above for loop insert measured labeling fraction data with adding
%         zero until the maximum length of the fields (i.e. k).
uit2 = app1.SimulatedMID_Table;
d2 = get(uit,'data');
uit2.ColumnName = mea_names;
sim=simulate(net_opt,xch_opt,input);
for i = 1:numel(mea_names)
    l = getfield(sim,mea_names{i})';
    l(end+1:k) = 0;
    d2(:,i) = l;
end
set(uit2, 'data', d2);
% Label bar graph
% Making string of legend component with for loop
% Finding maximum labeling index(eg. M+9) for each metabolites
maxLabIndex1 = [];
for j = 1:size(d,2)
    for i = size(d,1):-1:1
        if d(i,j) ~= 0
            maxLabIndex1(j) = i;
            break
        end
    end
end
totMaxLab1 = max(maxLabIndex1);
legendGen1 = {};
for i = 1:totMaxLab1
    legendGen1{i} = ['M+',num2str(i-1)];
end
maxLabIndex2 = [];
for j = 1:size(d2,2)
    for i = size(d2,1):-1:1
        if d2(i,j) ~= 0
            maxLabIndex2(j) = i;
            break
        end
    end
end
totMaxLab2 = max(maxLabIndex2);
legendGen2 = {};
for i = 1:totMaxLab2
    legendGen2{i} = ['M+',num2str(i-1)];
end
% Truncate bar handles to max labeling

uit.RowName = legendGen1;
uit2.RowName = legendGen2;

a = model.rxns';
uit3 = app1.FluxFreeEnergy_Table;
uit3.RowName = a;
net=model.kernel_net*net_opt;
xch=model.kernel_xch*xch_opt;
f=xch+max(0,net);
b=xch-min(0,net);
fe=8.314*temp/1000*log(b./f);
set(uit3,'data', [net xch f b fe]);
hold on
d=[];d2=[];
sim=simulate(net_opt,xch_opt,input);
mea_names = fieldnames(mea);
k = 0;
for i = 1:numel(mea_names)
    a = length(getfield(mea,mea_names{i}));
    if k <= a
        k = a;
    end
end
for i = 1:numel(mea_names)
    l = getfield(mea,mea_names{i})';
    l(end+1:k) = 0;
    d(:,i) = l;
end
for i = 1:numel(mea_names)
    l = getfield(sim,mea_names{i})';
    l(end+1:k) = 0;
    d2(:,i) = l;
end
dupmea = {};
for i = 1:2*size(mea_names,1)
    if i<=size(mea_names,1)
        dupmea{i} = mea_names{i};
    else
        dupmea{i} = mea_names{i-size(mea_names,1)};
    end
end
subplot(2,1,1)
b = bar([1:numel(mea_names)],d,'stacked');
ylim([0,1]);
[row1,col1] = find(d); [row2,col2] = find(d2); rowmax1 = max(row1); rowmax2 = max(row2); legendMaxInt = max(rowmax1, rowmax2);
legendName = {};
for i = 1:legendMaxInt
    a = num2str(i-1);
    legendName{i} = strcat('M+',a);
end
% Truncate bar handle vector
b_mod = b(1:legendMaxInt);
xticks([1:numel(mea_names)])
xticklabels(dupmea)
xtickangle(90);
titlePiece = 'Measured MID';
title(titlePiece);
subplot(2,1,2)
b = bar([1:numel(mea_names)],d2,'stacked');
ylim([0,1]);
legendName = {};
for i = 1:legendMaxInt
    a = num2str(i-1);
    legendName{i} = strcat('M+',a);
end
% Truncate bar handle vector
b_mod = b(1:legendMaxInt);
xticks([1:numel(mea_names)])
xticklabels(dupmea)
legend(fliplr(b_mod),fliplr(legendName),'Location','bestoutside');
xtickangle(90);
titlePiece = 'Simulated MID';
title(titlePiece);
saveas(gcf,'mea_sim','jpg');
imload = imread('mea_sim.jpg');
app1.OptimizedMIDImage.ImageSource = imload;
hold on;

%% Objective (score) function
function score = objfunc(simulate,net,xch,input,mea,fmea,var)
    [~,score]=varssr(simulate,net,xch,input,mea,fmea,var);
end

%% Out function
function stop = outfun(x,optimValues,state,parDataQueue)
    stop = false;
    h = figure('visible','off');
    switch state
        case 'init'
            hold on 
        case 'iter'
            % Saving live score data (only if parallel computing)
            if useParallel
                send(parDataQueue,[iteration,optimValues.iteration+1,optimValues.fval]);
                disp(['Parallel loop ',num2str(iteration),' iteration ',num2str(optimValues.iteration+1),' fval ',num2str(optimValues.fval)]);
            end
            
            % Running plot
            A = optimValues.iteration;
            net_opt=x(1:net_l);
            xch_opt=x(net_ll:end);
            currentScore = F([net_opt;xch_opt]);
            previousScore(A+2) = currentScore;
            if previousScore(A+2) <= previousScore(A+1)*0.9 || (norm(optimValues.stepsize) < 10^(-10) && optimValues.constrviolation < 10^(-30) )
                sim=simulate(net_opt,xch_opt,input);
                mea_names = fieldnames(mea);
                k = 0;
                for i = 1:numel(mea_names)
                    a = length(getfield(mea,mea_names{i}));
                    if k <= a
                        k = a;
                    end
                end
                for i = 1:numel(mea_names)
                    l = getfield(mea,mea_names{i})';
                    l(end+1:k) = 0;
                    d(:,i) = l;
                end
                for i = 1:numel(mea_names)
                    l = getfield(sim,mea_names{i})';
                    l(end+1:k) = 0;
                    d2(:,i) = l;
                end
                dupmea = {};
                for i = 1:2*size(mea_names,1)
                    if i<=size(mea_names,1)
                        dupmea{i} = mea_names{i};
                    else
                        dupmea{i} = mea_names{i-size(mea_names,1)};
                    end
                end
                subplot(1,2,1)
                b = bar([1:numel(mea_names)],d,'stacked');
                ylim([0,1]);
                [row1,col1] = find(d); [row2,col2] = find(d2); rowmax1 = max(row1); rowmax2 = max(row2); legendMaxInt = max(rowmax1, rowmax2);
                legendName = {};
                for i = 1:legendMaxInt
                    a = num2str(i-1);
                    legendName{i} = strcat('M+',a);
                end
                % Truncate bar handle vector
                b_mod = b(1:legendMaxInt);
                xticks([1:numel(mea_names)])
                xticklabels(dupmea)
                xtickangle(90);
                titlePiece = 'Measured MID';
                title(titlePiece);
                subplot(1,2,2)
                b = bar([1:numel(mea_names)],d2,'stacked');
                ylim([0,1]);
                legendName = {};
                for i = 1:legendMaxInt
                    a = num2str(i-1);
                    legendName{i} = strcat('M+',a);
                end
                % Truncate bar handle vector
                b_mod = b(1:legendMaxInt);
                xticks([1:numel(mea_names)])
                xticklabels(dupmea)
                legend(fliplr(b_mod),fliplr(legendName),'Location','bestoutside');
                xtickangle(90);
                scoreChar = num2str(currentScore); iterNum = num2str(A);
                titlePiece = 'Simulated MID';
                realTitle = strcat(iterNum,'th Iteration, ',' Current VSSR:', scoreChar);
                title({[titlePiece]
                    [realTitle]});
                saveas(gcf,'mea_sim','jpg');
                if app1.ShowImagesCheckBox.Value == 1
                    imload = imread('mea_sim.jpg');
                    app1.RealTimeOptImage.ImageSource = imload;
                hold on;
                end
            end                
    end
end

end
