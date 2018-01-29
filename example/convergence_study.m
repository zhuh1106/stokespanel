function convergence_study()
close all
nc = 1;
k = 3;
npsc = [4;5;5;8;5;5;4;5;5;8;5;5];   % number of panels along each side
npsi = [8;3;8;3];
[u0,sc,si] = panel_main_conv(nc, k, npsc, npsi);
numOfPanel = [0;0;0;0;0;0];
maxError = [0;0;0;0;0;0];


% with npsc
err = {};
for k = 1:2
    [u,~,~] = panel_main_conv(nc, k, npsc, npsi);
    err{k}.err = abs(u - u0);
    err{k}.np = k*(sum(npsc)+sum(npsi));
    figure()
    imagesc(1:80,1:40,log10(abs(err{k}.err))),colorbar;
end
numOfPanel(1) = err{1}.np;
numOfPanel(2) = err{2}.np;
maxError(1) = max(max(err{1}.err));
maxError(2) = max(max(err{2}.err));

% with npsc2
nc = 1;
npsc2 = [3;3;3;3;3;3;3;3;3;3;3;3];
npsi2 = [3;3;3;3];
err2 = {};
for k = 1:4
    [u,~,~] = panel_main_conv(nc, k, npsc2, npsi2);
    err2{k}.err = abs(u - u0);
    err2{k}.np = k*(sum(npsc2)+sum(npsi2));
    figure()
    imagesc(1:80,1:40,log10(abs(err2{k}.err))),colorbar;
end
numOfPanel(3) = err2{1}.np;
numOfPanel(4) = err2{2}.np;
numOfPanel(5) = err2{3}.np;
numOfPanel(6) = err2{4}.np;

maxError(3) = max(max(err2{1}.err));
maxError(4) = max(max(err2{2}.err));
maxError(5) = max(max(err2{3}.err));
maxError(6) = max(max(err2{4}.err));

maxError;
numOfPanel;

