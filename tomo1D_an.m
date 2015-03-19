%---------------------------------------------------------------------
%%1D-Tomographie mit Ableitungsnorm - Übungsblatt 5
%

function [] = tomo1D_an(inputdata, delta_x,lbar)

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
hk = 0.5*(xs-xr);                       %midpoint between source and resiver
Dk = abs(xr-xs);                        %distance source resiever

%-----------------------------------------------------------------------
%
xm = linspace(0,lbar,num_nodes);		% discretise -> modle points xm
gk = zeros(num_nodes,num_data);			%matrix of representatives
%calculate matrix gk 
for j=1:num_data
    for i= 1:num_nodes
        if i*delta_x <= xs(j)
            gk(i,j)=Dk(j)*(2*lbar-hk(j));
        elseif i*delta_x > xs(j) && i*delta_x <= xr(j)
            gk(i,j)=Dk(j)*(2*lbar-hk(j))-0.5*(i*delta_x-hk(j)+0.5*Dk(j)).^2;
        elseif i*delta_x > xr(j)
            gk(i,j)=Dk(j)*(2*lbar-i*delta_x);
        end 
    end
end
Rk = zeros(num_nodes,num_data);         %needed for the new definitaion of the "Innenprodukt", integral of the box function
for j=1:num_data
    for i= 1:num_nodes
        if i*delta_x <= xs(j)
            Rk(i,j)=0;
        elseif i*delta_x > xs(j) && i*delta_x <= xr(j)
            Rk(i,j)=-i*delta_x+hk(j)-0.5*Dk(j);
        elseif i*delta_x > xr(j)
            Rk(i,j)=-Dk(j);
        end 
    end
end

G = Rk'*Rk*delta_x;                     %calculate  Gram-matrix
for j=num_data
    for i=num_data
        G(i,j)=G(i,j)+lbar*Dk(j)*Dk(i); %Gram matrix calculatet with new defined "Innenprodukt"
    end
end
alpha = linsolve(G,tt);					%solve system of equations  G*alpha=d
%calculate slowness function
p = 100.*gk*alpha;						%slowness function (s/km)

%Plot
plot(xm,p);
title ('1D-Tomographie')
xlabel ('Abtastinterval [cm]'); ylabel('Slowness [s/km]');
