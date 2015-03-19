%% ---------------------------------------------------------------------
% 1D-Tomographie mit L2-Norm und Messfehlern - Übungsblatt 6

function [] = tomo1D_misfit(inputdata,delta_x,lbar)

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
er = inputdata(:,5);                    %errors in travel time
%% -----------------------------------------------------------------------

xm = linspace(0,lbar,num_nodes);		% discretise -> modle points xm
gk = zeros(num_nodes,num_data);			%matrix of representatives
%calculate matrix gk 
for j=1:num_data
	gk(:,j) = xm >= min(xs(j),xr(j)) & xm < max(xs(j),xr(j));
end 
G = gk'*gk*delta_x;						%calculate Gram-matrix
%% -----------------------------------------------------------------------
%loop over Lagrange parametr nu
nu_l=-10;                               %lower limit of nu
nu_u=10;                                %upper limit of nu
nnu=200;                                %number of nu values
Ner=6;
nu = linspace(nu_l,nu_u,nnu);
mf = zeros(nnu,1);                      %misfit
mN = zeros(nnu,1);                      %modle norm
E2 = er.^2;
for i=1:nnu
    alpha=linsolve(G+realpow(10.,-nu(i))*diag(E2),tt);   %solwing linear system
    mf(i)= norm((tt-G*alpha)./er)^2;
    mN(i)= alpha'*G*alpha;
end
f1 = figure('Name','Misift vs. Lagrange parameter');
plot(nu,log10(mf));
xlabel('Lagrange parameter'); ylabel('Log10 Misfit');
f2 = figure('Name','Misift vs. model norm');
plot(mN,mf);
xlabel('Model norm'); ylabel('Misfit');

%% bisection
f3=figure('Name','1D-Tomographie');
for j = 1:Ner
    er_new = er*j^2;
    fprintf('Error: %10.4f\n',er_new(1));
    nu1 = nu_l; nu2 = nu_u;
    E2 = er_new.^2;
    while (nu2-nu1) > 0.001
        nu_m=0.5*(nu2+nu1);
        alfa = linsolve(G+realpow(10.,-nu_m)*diag(E2),tt);
        res = tt-G*alfa;
        mf_new = dot(res./er_new,res./er_new);
        if (mf_new > num_data)
            nu1 = nu;
        else
            nu2 = nu;
        end
        fprintf('%10.4f %10.4f\n',nu_m,mf_new);
    end
    p = 100.*gk*alpha;
    plot(xm,p);
    xlabel ('Abtastinterval [cm]'); ylabel('Slowness [s/km]');
end

%% calculate slowness function
						%slowness function (s/km)


