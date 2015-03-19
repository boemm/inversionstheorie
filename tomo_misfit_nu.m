%---------------------------------------------------------------------
% 1D tomography L2-norm with errors, misifit versus Lagrange parameter
%
obsfile='data/data.in';             % file with data
lbar = 20;                          % length of bar (20 cm)
nm = 2000;                          % number of model points
nnu = 24;                           % number of nu-values
lgnul = -8;                         % lower limit of lgnu
lgnur = +3;                         % upper limit of lgnu
%-----------------------------------------------------------------
fid = fopen(obsfile);                   % open file
obsdata = fscanf(fid,'%d %f %f %f %f'); % read in file
fclose(fid);                            % close file
xs = obsdata(2:5:length(obsdata));      % get source points (cm)
xr = obsdata(3:5:length(obsdata));      % get receiver points (cm)
tt = obsdata(4:5:length(obsdata));      % get travel times (ms)
er = obsdata(5:5:length(obsdata));      % get travel time errors (ms)
nobs = length(xr);                      % number of data
xm = linspace(0,lbar,nm);               % discretize bar (model points xm)
dx = lbar/(nm-1);                       % interval on bar (cm)
g = zeros(nm,nobs);                     % calculate representers
for j=1:nobs
    g(:,j) = xm >= min(xs(j),xr(j)) & xm < max(xs(j),xr(j));
end
gm = g'*g*dx;                           % Gram matrix
%-------------------------------------------------------------------
% Loop over Lagrange parameter nu, use log nu
lgnu = linspace(lgnul,lgnur,nnu);
er2 = er.^2;
misfit = zeros(nnu,1); modnorm = zeros(nnu,1);
for j=1:nnu
    alfa = linsolve(gm+realpow(10.,-lgnu(j))*diag(er2),tt);  % solve linear system
    misfit(j) = norm((tt-gm*alfa)./er)^2;
    modnorm(j)=alfa'*gm*alfa;
    fprintf('%6.2f %10.2e %10.2e\n',lgnu(j),misfit(j),modnorm(j));
end
f1 = figure('Name','Misift versus Lagrange parameter');
plot(lgnu,log10(misfit));
xlabel('Lagrange parameter'); ylabel('Log10 Misfit');
f2 = figure('Name','Misift versus model norm');
plot(modnorm,misfit);
xlabel('Model norm'); ylabel('Misfit');

