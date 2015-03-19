%---------------------------------------------------------------------
%%1D-Tomographie - Übungsblatt 4
%

function [] = tomo1D(inputdata,delta_x,lbar)

%%Parameter
%inputdata					%Tomo. data (number, source, resiever,travel time)
%delta_x					%sampling interval
%lbar						%bare lenght [cm]
%---------------------------------------------------------------------

num_nodes = lbar/delta_x;				%number modle points
inputdata = dlmread(inputdata);			%reading data file
num_data = length(inputdata(:,1));		%number data
xs = inputdata(:,2);					%source
xr = inputdata(:,3);					%resiever
tt = inputdata(:,4);					%travel time

%-----------------------------------------------------------------------
%
xm = linspace(0,lbar,num_nodes);		% discretise -> modle points xm
gk = zeros(num_nodes,num_data);			%matrix of representatives
%calculate matrix gk 
for j=1:num_data
	gk(:,j) = xm >= min(xs(j),xr(j)) & xm < max(xs(j),xr(j));
end 

G = gk'*gk*delta_x;						%calculate  Gram-matrix
alpha = linsolve(G,tt);					%solve linear system  G*alpha=d
%calculate slowness function
p = 100.*gk*alpha;						%slowness function (s/km)

%Plot
plot(xm,p);
title ('1D-Tomographie')
xlabel ('Abtastinterval [cm]'); ylabel('Slowness [s/km]');
