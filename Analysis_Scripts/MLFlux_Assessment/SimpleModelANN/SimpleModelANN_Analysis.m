%% Analytical solution (Fig.1d)

[B,D] = meshgrid(0.5:0.01:1,0:0.01:1);
Btemp = B;
Dtemp = D;
D(Dtemp>(Btemp+0.005)) = NaN;
B(Dtemp>(Btemp+0.005)) = NaN;
clear Btemp Dtemp

v3 = (1-D)./(2.*B-D);
v2r = (1-B)./(B-D);
v2f = 1-2.*v3+v2r;
v2_xch = min(v2r,v2f);
 
figure(1)
clf
hold on
contourf(B,D,v3,0:0.05:1)
title("Analytical v3/v1")
xlabel("B1/B")
ylabel("D1/D")
zlabel("v3/v1")
set(gca,'ColorScale','linear')
set(gca,'TickDir','out')
clim([0,1])
colorbar

figure(2)
clf
hold on
title("Analytical v2r")

contourf(B,D,v2r,[0:0.1:1,2:1:5])
xlabel("B1/B")
ylabel("D1/D")
zlabel("v3/v1")
set(gca,'ColorScale','log')
set(gca,'TickDir','out')
clim([1E-4,10])
colorbar('Ticks', [0:0.1:1,2:1:5],'TickLabels', ...
    {'<0.1','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9', ...
    '1','2','3','4','5','6','7','8','9','>=10'});

%% ANN Solution (Fig.1g, SuppFig.1c-d)
% To change between flux sampling configurations, 'labeldir' in the the
% function 'performNN' must be changed. logUniXch is the final one used in
% Fig. 1g

kernel = load("KernelNet_Simple.txt");

[B_ANN,D_ANN] = meshgrid(0.5:0.01:1,0:0.01:1);
Btemp = B_ANN;
Dtemp = D_ANN;
D_ANN(Dtemp>Btemp) = NaN;
B_ANN(Dtemp>Btemp) = NaN;
clear Btemp Dtemp

v3_ANN = zeros(size(B_ANN));
v2x_ANN = zeros(size(B_ANN));
v2r_ANN = zeros(size(B_ANN));
v2n_ANN = zeros(size(B_ANN));
for i = 1:size(v3_ANN,1)
    for j = 1:size(v3_ANN,2)
        freeFluxes = performNN(getLabels(B_ANN(i,j),D_ANN(i,j)));
        fullFluxes = (kernel*freeFluxes(1:2,:));

        v3_ANN(i,j) = fullFluxes(3)./fullFluxes(1);
        v2x_ANN(i,j) = freeFluxes(3)./fullFluxes(1);
        v2n_ANN(i,j) = 1-2*v3_ANN(i,j);
        if v2n_ANN(i,j)>0
            v2r_ANN(i,j) = v2x_ANN(i,j);
        else
            v2r_ANN(i,j) = v2x_ANN(i,j)-v2n_ANN(i,j);
        end

    end
end

figure(3)
clf
contourf(B_ANN,D_ANN,v3_ANN,0:0.05:1)
title("ANN v3/v1")
xlabel("B1/B")
ylabel("D1/D")
zlabel("v3/v1")
set(gca,'ColorScale','linear')
set(gca,'TickDir','out')
clim([0,1])
colorbar

figure(4)
clf
hold on
title("ANN v2r")

contourf(B_ANN,D_ANN,v2r_ANN,[0:0.1:1,2:1:5])
xlabel("B1/B") 
ylabel("D1/D")
zlabel("v3/v1")
set(gca,'ColorScale','log')
set(gca,'TickDir','out')
clim([1E-4,10])
colorbar('Ticks', [0:0.1:1,2:1:5],'TickLabels', ...
    {'<0.1','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9', ...
    '1','2','3','4','5','6','7','8','9','>=10'});

%% Functions
function fluxes = performNN(labels)
    %labeldir = linUniXch,logUniXch,logNormXch,DistanceWeighted,Uniform
    labeldir = "linUniXch";
    W_In_1 = load(labeldir+"/weights_Dense_1.txt");
    W_1_2 = load(labeldir+"/weights_Dense_2.txt");
    W_2_3 = load(labeldir+"/weights_Dense_3.txt");
    W_3_4 = load(labeldir+"/weights_Dense_4.txt");
    W_4_5 = load(labeldir+"/weights_Dense_5.txt");
    W_5_Out = load(labeldir+"/weights_output.txt");
    B_In_1 = load(labeldir+"/biases_Dense_1.txt");
    B_1_2 = load(labeldir+"/biases_Dense_2.txt");
    B_2_3 = load(labeldir+"/biases_Dense_3.txt");
    B_3_4 = load(labeldir+"/biases_Dense_4.txt");
    B_4_5 = load(labeldir+"/biases_Dense_5.txt");
    B_5_Out = load(labeldir+"/biases_output.txt");

    neurons1 = max(0,W_In_1*labels+B_In_1);
    neurons2 = max(0,W_1_2*neurons1+B_1_2);
    neurons3 = max(0,W_2_3*neurons2+B_2_3);
    neurons4 = max(0,W_3_4*neurons3+B_3_4);
    neurons5 = max(0,W_4_5*neurons4+B_4_5);
    allfluxes = W_5_Out*neurons5+B_5_Out;

    netfluxes = allfluxes(1:2,:);

    xchfluxes = allfluxes(3,:);
    xchfluxes(xchfluxes<0) = 0;

    fluxes = [netfluxes;xchfluxes];
end

function allLabels = getLabels(B1,D1)
    allLabels = zeros(20,1);
    allLabels(1) = 1-B1;
    allLabels(2) = B1;
    allLabels(5) = 1;
    allLabels(9) = 1-D1;
    allLabels(10) = D1;
    allLabels(13) = 1-2/3*B1;
    allLabels(14) = 2/3*B1;
    allLabels(17) = allLabels(9);
    allLabels(18) = allLabels(10);
end