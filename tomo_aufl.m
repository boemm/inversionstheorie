%function [] = tomo_aufl(obsfile, dx, lbar)
% 1D tomography
% parameter:
% obsfile             path to file with data
% dx                  interval on bar (cm)
% lbar                length of bar


dx=0.1;
lbar=20;
obsfile='data/data.in';
%% read
num_nodes = lbar/dx;                            % number of model points
xm = linspace(0,lbar,num_nodes);                % discretize bar (model points xm)
obsdata=dlmread(obsfile);                       % read data
num_data = length(obsdata(:,1));                %number data                    
xs = obsdata(:,2);                              % get source points (cm)
xr = obsdata(:,3);                              % get receiver points (cm)
tt = obsdata(:,4);                              %travel time

background=0.01;
peak=0.02;
peakpos=5;
peakintervall=1;
peakstart=peakpos-1/2*peakintervall;
peakend=peakpos+1/2*peakintervall;

for i=1:num_nodes
    
    if i*dx < peakstart 
        slowness(i)=background;
    elseif i*dx > peakend
         slowness(i)=background;
    else
        slowness(i)=peak;
    end
end
%%
g = zeros(num_nodes,num_data);  % calculate representers
for j=1:num_data
    g(:,j) = xm >= min(xs(j),xr(j)) & xm < max(xs(j),xr(j));
end
%--------------------------------------------------------------------------
synfile='data/data_syn.in';

ttsyn = slowness*g*dx/100.;                      %predicted travel times by test model
fid = fopen(synfile, 'w');
%fprintf(fid, '%s  %s  %s\n', 'source','resiever','syn. traveltime');
fprintf(fid,'%8.3f\n',ttsyn);
fclose(fid);

% plot synthetic data
f1 = figure('Name','Slowness Modell');
plot(xm,slowness);
xlabel ('Abtastinterval [cm]'); ylabel('Slowness [s/km]');
%% Annahme exakte Daten und L2-Norm
syndata=dlmread(synfile);                  % read data
ttsyn1 = syndata(:,1);

G = g'*g*dx;                               %calculate  Gram-matrix
alpha = linsolve(G,tt);
p=100.*g*alpha;

G1 = g'*g*dx;                               %calculate  Gram-matrix
alpha1 = linsolve(G1,ttsyn1);				%solve linear system  G*alpha=d
%calculate slowness function
p1 = 100.*g*alpha1;						%slowness function (s/km)

%plot Modelldaten vs. echte Daten mit L2-Norm
f2=figure('Name', 'Modelldaten vs. echte Daten');
subplot(2,1,1); plot(xm,p1);
title ('1D-Tomographie Modellaten und L2-Norm')
xlabel ('Abtastinterval [cm]'); ylabel('Slowness [s/km]');
figure2=subplot(2,1,2); plot(xm,p);
title ('1D-Tomographie echte Daten und L2-Norm')
xlabel ('Abtastinterval [cm]'); ylabel('Slowness [s/km]');

%% Ableitungsnorm


