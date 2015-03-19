%---------------------------------------------------------------------
% 1D tomography L2-norm with errors, bisection algorithm for correct nu
%
obsfile='data/messorig.in';                  % file with data
lbar = 15;                          % length of bar (15 cm)
nm = 151;                           % number of model points
nnu = 24;                           % number of nu-values
lgnul = -8;                         % lower limit of lgnu
lgnur = +6;                         % upper limit of lgnu
nerr = 6;                           % number of different errors
%-----------------------------------------------------------------
fid = fopen(obsfile);               % open file
obsdata = fscanf(fid,'%d %f %f %f %f'); % read in file
fclose(fid);                        % close file
xs = obsdata(2:5:length(obsdata));  % get source points (cm)
xr = obsdata(3:5:length(obsdata));  % get receiver points (cm)
tt = obsdata(4:5:length(obsdata));  % get travel times (ms)
er = obsdata(5:5:length(obsdata));  % get travel time errors (ms)
nobs = length(xr);                  % number of data
xm = linspace(0,lbar,nm);           % discretize bar (model points xm)
dx = lbar/(nm-1);                   % interval on bar (cm)
g = zeros(nm,nobs);                 % calculate representers
for j=1:nobs
    g(:,j) = xm >= min(xs(j),xr(j)) & xm < max(xs(j),xr(j));
end
gm = g'*g*dx;                       % Gram matrix
%-----------------------------------------------------------------
% start bisection
f1 = figure('Name','Slowness model');
for j = 1:nerr
    ernew = er*j^2;
    fprintf('Error: %10.4f\n',ernew(1));
    lgnu1 = lgnul; lgnu2 = lgnur;
    er2 = ernew.^2;
    while (lgnu2-lgnu1) > 0.001
        lgnum=0.5*(lgnu2+lgnu1);
        alfa = linsolve(gm+realpow(10.,-lgnum)*diag(er2),tt);
        res = tt-gm*alfa;
        mf = dot(res./ernew,res./ernew);
        if (mf > nobs)
            lgnu1 = lgnum;
        else
            lgnu2 = lgnum;
        end
        fprintf('%10.4f %10.4f\n',lgnum,mf);
    end
    p = 100.*g*alfa;                    % evaluate slowness function (s/km)
    pl = plot(xm,p);                    % plot slowness
    xlabel('Offset [cm]'); ylabel('Slowness [s/km]'); title('Slowness model')
    hold on
    set(pl,'LineWidth',j);
    ttsyn = gm*alfa;                    % predicted travel times by model (ms)
    for k = 1:nobs;                     % compare
        fprintf('%8.3f %8.3f\n',tt(k),ttsyn(k));
    end
end

