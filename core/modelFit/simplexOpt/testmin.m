function [] = testmin()
% test minimization results


for ii = 1:20; fprintf('\n'); end

warning('off','optim:fmincon:SwitchingToMediumScale')
dbg = 0;

% pause time
pse = 5;
np = 1001;



%% test unconstrained
fprintf('Unconstrained test - https://en.wikipedia.org/wiki/Test_functions_for_optimization\n\n\n')

% sphere
fn = 'Sphere';
fun = @(x)(sum(x.^2));
exct = [0 0];
ic = [-2 2];
% straight simplex
[xs, fs, cs] = simplexmin(fun,ic,[],[],[],[],[],[],[],dbg);
% simplex function min
[x, val, cnt] = minfun(fun,ic,[],[],[],[],[],[],[],dbg);
% fmin function
[xm, fm, ~, out] = fminsearch(fun,ic);
% fm = feval(fun,xm);
fprintf('\n');
fprintf('%s	(%s)	',fn,func2str(fun));
fprintf('Exact soln:	'); fprintf('%g	',exct); fprintf('\n');
fprintf('ic:	'); fprintf('%g	',ic); fprintf('\n');
fprintf('fmin:	%g	%g evals	@	',fm,out.funcCount); fprintf('%g	',xm); 
fprintf('error:	'); fprintf('%g	',exct-xm); fprintf('\n');
fprintf('simplex:	%g	%g evals	@	',fs,cs); fprintf('%g,	',xs); 
fprintf('error:	'); fprintf('%g	',exct-xs'); fprintf('\n');
fprintf('minfun:	%g	%g evals	@	',val,cnt); fprintf('%g,	',x); 
fprintf('error:	'); fprintf('%g	',exct-x');
fprintf('\n');
f = figure; set(f,'color','w'); hold on; grid on;
xlabel('X1'); ylabel('X2'); colorbar; title(char('sphere',char(fun)))
xv = linspace(min([ic(1) x(1) xs(1) xm(1)]),max([ic(1) x(1) xs(1) xm(1)]),np)';
yv = linspace(min([ic(2) x(2) xs(2) xm(2)]),max([ic(2) x(2) xs(2) xm(2)]),np);
[xx, yy] = meshgrid(xv,yv);
zz = xx;
for ii = 1:np^2
	zz(ii) = feval(fun,[xx(ii) yy(ii)]);
end
imagesc(xv,yv,zz); axis image
temp = plot(xs(1),xs(2),'k*'); set(temp,'markersize',10)
temp = plot(x(1),x(2),'ko'); set(temp,'markersize',10)
temp = plot(xm(1),xm(2),'kd'); set(temp,'markersize',10)
legend('Simplex','Min','Fminsearch')
pause(pse); fprintf('\n\n');


% sphere
fn = 'Sphere';
fun = @(x)(sum(x.^2));
exct = [0 0 0 0];
ic = [-2 2 9 4];
% straight simplex
[xs, fs, cs] = simplexmin(fun,ic,[],[],[],[],[],[],[],dbg);
% simplex function min
[x, val, cnt] = minfun(fun,ic,[],[],[],[],[],[],[],dbg);
% fmin function
[xm, fm, ~, out] = fminsearch(fun,ic);
% fm = feval(fun,xm);
fprintf('\n');
fprintf('%s	(%s)	',fn,func2str(fun));
fprintf('Exact soln:	'); fprintf('%g	',exct); fprintf('\n');
fprintf('ic:	'); fprintf('%g	',ic); fprintf('\n');
fprintf('fmin:	%g	%g evals	@	',fm,out.funcCount); fprintf('%g	',xm); 
fprintf('error:	'); fprintf('%g	',exct-xm); fprintf('\n');
fprintf('simplex:	%g	%g evals	@	',fs,cs); fprintf('%g,	',xs); 
fprintf('error:	'); fprintf('%g	',exct-xs'); fprintf('\n');
fprintf('minfun:	%g	%g evals	@	',val,cnt); fprintf('%g,	',x); 
fprintf('error:	'); fprintf('%g	',exct-x');
fprintf('\n');
f = figure; set(f,'color','w'); hold on; grid on;
xlabel('X1'); ylabel('X2'); colorbar; title(char('sphere',char(fun)))
xv = linspace(min([ic(1) x(1) xs(1) xm(1)]),max([ic(1) x(1) xs(1) xm(1)]),np)';
yv = linspace(min([ic(2) x(2) xs(2) xm(2)]),max([ic(2) x(2) xs(2) xm(2)]),np);
[xx, yy] = meshgrid(xv,yv);
zz = xx;
for ii = 1:np^2
	zz(ii) = feval(fun,[xx(ii) yy(ii)]);
end
imagesc(xv,yv,zz); axis image
temp = plot(xs(1),xs(2),'k*'); set(temp,'markersize',10)
temp = plot(x(1),x(2),'ko'); set(temp,'markersize',10)
temp = plot(xm(1),xm(2),'kd'); set(temp,'markersize',10)
legend('Simplex','Min','Fminsearch')
pause(pse); fprintf('\n\n');



% sphere
fn = 'Offset Sphere';
fun = @(x)((x(1)+4)^2 + (x(2)-5)^2 + (x(3)+2)^2 + (x(4)-9)^2);
exct = [-4 5 -2 9];
ic = [-8 2 -7 4];
% straight simplex
[xs, fs, cs] = simplexmin(fun,ic,[],[],[],[],[],[],[],dbg);
% simplex function min
[x, val, cnt] = minfun(fun,ic,[],[],[],[],[],[],[],dbg);
% fmin function
[xm, fm, ~, out] = fminsearch(fun,ic);
% fm = feval(fun,xm);
fprintf('\n');
fprintf('%s	(%s)	',fn,func2str(fun));
fprintf('Exact soln:	'); fprintf('%g	',exct); fprintf('\n');
fprintf('ic:	'); fprintf('%g	',ic); fprintf('\n');
fprintf('fmin:	%g	%g evals	@	',fm,out.funcCount); fprintf('%g	',xm); 
fprintf('error:	'); fprintf('%g	',exct-xm); fprintf('\n');
fprintf('simplex:	%g	%g evals	@	',fs,cs); fprintf('%g,	',xs); 
fprintf('error:	'); fprintf('%g	',exct-xs'); fprintf('\n');
fprintf('minfun:	%g	%g evals	@	',val,cnt); fprintf('%g	',x); 
fprintf('error:	'); fprintf('%g	',exct-x');
fprintf('\n');
pause(pse); fprintf('\n\n');



% Rastrigin
fn = 'Rastrigin';
fun = @(x)(20 + sum(x.^2 - 10*cos(2 * pi * x)));
exct = [0 0];
ic = [-8 4];
% straight simplex
[xs, fs, cs] = simplexmin(fun,ic,[],[],[],[],[],[],[],dbg);
% simplex function min
[x, val, cnt] = minfun(fun,ic,[],[],[],[],[],[],[],dbg);
% fmin function
[xm, fm, ~, out] = fminsearch(fun,ic);
% fm = feval(fun,xm);
fprintf('\n');
fprintf('%s	(%s)	',fn,func2str(fun));
fprintf('Exact soln:	'); fprintf('%g	',exct); fprintf('\n');
fprintf('ic:	'); fprintf('%g	',ic); fprintf('\n');
fprintf('fmin:	%g	%g evals	@	',fm,out.funcCount); fprintf('%g	',xm); 
fprintf('error:	'); fprintf('%g	',exct-xm); fprintf('\n');
fprintf('simplex:	%g	%g evals	@	',fs,cs); fprintf('%g,	',xs); 
fprintf('error:	'); fprintf('%g	',exct-xs'); fprintf('\n');
fprintf('minfun:	%g	%g evals	@	',val,cnt); fprintf('%g	',x); 
fprintf('error:	'); fprintf('%g	',exct-x');
fprintf('\n');
f = figure; set(f,'color','w'); hold on; grid on;
xlabel('X1'); ylabel('X2'); colorbar; title(char('Rastrigin',char(fun)))
xv = linspace(min([ic(1) x(1) xs(1) xm(1)]),max([ic(1) x(1) xs(1) xm(1)]),np)';
yv = linspace(min([ic(2) x(2) xs(2) xm(2)]),max([ic(2) x(2) xs(2) xm(2)]),np);
[xx, yy] = meshgrid(xv,yv);
zz = xx;
for ii = 1:np^2
	zz(ii) = feval(fun,[xx(ii) yy(ii)]);
end
imagesc(xv,yv,zz); axis image
temp = plot(xs(1),xs(2),'k*'); set(temp,'markersize',10)
temp = plot(x(1),x(2),'ko'); set(temp,'markersize',10)
temp = plot(xm(1),xm(2),'kd'); set(temp,'markersize',10)
legend('Simplex','Min','Fminsearch')
pause(pse); fprintf('\n\n');


% Rastrigin
fn = 'Rastrigin';
fun = @(x)(40 + sum(x.^2 - 10*cos(2 * pi * x)));
exct = [0 0 0 0];
ic = [-8 4 1 -4];
% straight simplex
[xs, fs, cs] = simplexmin(fun,ic,[],[],[],[],[],[],[],dbg);
% simplex function min
[x, val, cnt] = minfun(fun,ic,[],[],[],[],[],[],[],dbg);
% fmin function
[xm, fm, ~, out] = fminsearch(fun,ic);
% fm = feval(fun,xm);
fprintf('\n');
fprintf('%s	(%s)	',fn,func2str(fun));
fprintf('Exact soln:	'); fprintf('%g	',exct); fprintf('\n');
fprintf('ic:	'); fprintf('%g	',ic); fprintf('\n');
fprintf('fmin:	%g	%g evals	@	',fm,out.funcCount); fprintf('%g	',xm); 
fprintf('error:	'); fprintf('%g	',exct-xm); fprintf('\n');
fprintf('simplex:	%g	%g evals	@	',fs,cs); fprintf('%g,	',xs); 
fprintf('error:	'); fprintf('%g	',exct-xs'); fprintf('\n');
fprintf('minfun:	%g	%g evals	@	',val,cnt); fprintf('%g	',x); 
fprintf('error:	'); fprintf('%g	',exct-x');
fprintf('\n');
f = figure; set(f,'color','w'); hold on; grid on;
xlabel('X1'); ylabel('X2'); colorbar; title(char('Rastrigin',char(fun)))
xv = linspace(min([ic(1) x(1) xs(1) xm(1)]),max([ic(1) x(1) xs(1) xm(1)]),np)';
yv = linspace(min([ic(2) x(2) xs(2) xm(2)]),max([ic(2) x(2) xs(2) xm(2)]),np);
[xx, yy] = meshgrid(xv,yv);
zz = xx;
for ii = 1:np^2
	zz(ii) = feval(fun,[xx(ii) yy(ii)]);
end
imagesc(xv,yv,zz); axis image
temp = plot(xs(1),xs(2),'k*'); set(temp,'markersize',10)
temp = plot(x(1),x(2),'ko'); set(temp,'markersize',10)
temp = plot(xm(1),xm(2),'kd'); set(temp,'markersize',10)
legend('Simplex','Min','Fminsearch')
pause(pse); fprintf('\n\n');



% Ackley
fn = 'Ackley';
fun = @(x)(-20*exp(-0.2 * sqrt(0.5 * (x(1)^2 + x(2)^2))) - ...
	exp(0.5 * (cos(2*pi*x(1)) + cos(2*pi*x(2)))) + exp(1) + 20);
exct = [0 0];
ic = [3 -3];
% straight simplex
[xs, fs, cs] = simplexmin(fun,ic,[],[],[],[],[],[],[],dbg);
% simplex function min
[x, val, cnt] = minfun(fun,ic,[],[],[],[],[],[],[],dbg);
% fmin function
[xm, fm, ~, out] = fminsearch(fun,ic);
% fm = feval(fun,xm);
fprintf('\n');
fprintf('%s	(%s)	',fn,func2str(fun));
fprintf('Exact soln:	'); fprintf('%g	',exct); fprintf('\n');
fprintf('ic:	'); fprintf('%g	',ic); fprintf('\n');
fprintf('fmin:	%g	%g evals	@	',fm,out.funcCount); fprintf('%g	',xm); 
fprintf('error:	'); fprintf('%g	',exct-xm); fprintf('\n');
fprintf('simplex:	%g	%g evals	@	',fs,cs); fprintf('%g,	',xs); 
fprintf('error:	'); fprintf('%g	',exct-xs'); fprintf('\n');
fprintf('minfun:	%g	%g evals	@	',val,cnt); fprintf('%g	',x); 
fprintf('error:	'); fprintf('%g	',exct-x');
fprintf('\n');
f = figure; set(f,'color','w'); hold on; grid on;
xlabel('X1'); ylabel('X2'); colorbar; title(char('Ackley',char(fun)))
xv = linspace(min([ic(1) x(1) xs(1) xm(1)]),max([ic(1) x(1) xs(1) xm(1)]),np)';
yv = linspace(min([ic(2) x(2) xs(2) xm(2)]),max([ic(2) x(2) xs(2) xm(2)]),np);
[xx, yy] = meshgrid(xv,yv);
zz = xx;
for ii = 1:np^2
	zz(ii) = feval(fun,[xx(ii) yy(ii)]);
end
imagesc(xv,yv,zz); axis image
temp = plot(xs(1),xs(2),'k*'); set(temp,'markersize',10)
temp = plot(x(1),x(2),'ko'); set(temp,'markersize',10)
temp = plot(xm(1),xm(2),'kd'); set(temp,'markersize',10)
legend('Simplex','Min','Fminsearch')
pause(pse); fprintf('\n\n');


% rosenbrock
fn = 'rosenbrock';
fun = @(x)(sum(100*((x(2:end) - (x(1:end-1)).^2).^2) + (1 - x(1:end-1)).^2));
exct = [1 1];
ic = [3 -3];
% straight simplex
[xs, fs, cs] = simplexmin(fun,ic,[],[],[],[],[],[],[],dbg);
% simplex function min
[x, val, cnt] = minfun(fun,ic,[],[],[],[],[],[],[],dbg);
% fmin function
[xm, fm, ~, out] = fminsearch(fun,ic);
% fm = feval(fun,xm);
fprintf('\n');
fprintf('%s	(%s)	',fn,func2str(fun));
fprintf('Exact soln:	'); fprintf('%g	',exct); fprintf('\n');
fprintf('ic:	'); fprintf('%g	',ic); fprintf('\n');
fprintf('fmin:	%g	%g evals	@	',fm,out.funcCount); fprintf('%g	',xm); 
fprintf('error:	'); fprintf('%g	',exct-xm); fprintf('\n');
fprintf('simplex:	%g	%g evals	@	',fs,cs); fprintf('%g,	',xs); 
fprintf('error:	'); fprintf('%g	',exct-xs'); fprintf('\n');
fprintf('minfun:	%g	%g evals	@	',val,cnt); fprintf('%g	',x); 
fprintf('error:	'); fprintf('%g	',exct-x');
fprintf('\n');
f = figure; set(f,'color','w'); hold on; grid on;
xlabel('X1'); ylabel('X2'); colorbar; title(char('rosenbrok',char(fun)))
xv = linspace(min([ic(1) x(1) xs(1) xm(1)]),max([ic(1) x(1) xs(1) xm(1)]),np)';
yv = linspace(min([ic(2) x(2) xs(2) xm(2)]),max([ic(2) x(2) xs(2) xm(2)]),np);
[xx, yy] = meshgrid(xv,yv);
zz = xx;
for ii = 1:np^2
	zz(ii) = feval(fun,[xx(ii) yy(ii)]);
end
imagesc(xv,yv,zz); axis image
temp = plot(xs(1),xs(2),'k*'); set(temp,'markersize',10)
temp = plot(x(1),x(2),'ko'); set(temp,'markersize',10)
temp = plot(xm(1),xm(2),'kd'); set(temp,'markersize',10)
legend('Simplex','Min','Fminsearch')
pause(pse); fprintf('\n\n');


% beals
fn = 'beale';
fun = @(x)((1.5 - x(1) + x(1)*x(2))^2 + (2.25 - x(1) + x(1)*x(2)^2)^2 + ...
	(2.625 - x(1) + x(1)*x(2)^3)^2);
exct = [3 0.5];
ic = [2 2];
% straight simplex
[xs, fs, cs] = simplexmin(fun,ic,[],[],[],[],[],[],[],dbg);
% simplex function min
[x, val, cnt] = minfun(fun,ic,[],[],[],[],[],[],[],dbg);
% fmin function
[xm, fm, ~, out] = fminsearch(fun,ic);
% fm = feval(fun,xm);
fprintf('\n');
fprintf('%s	(%s)	',fn,func2str(fun));
fprintf('Exact soln:	'); fprintf('%g	',exct); fprintf('\n');
fprintf('ic:	'); fprintf('%g	',ic); fprintf('\n');
fprintf('fmin:	%g	%g evals	@	',fm,out.funcCount); fprintf('%g	',xm); 
fprintf('error:	'); fprintf('%g	',exct-xm); fprintf('\n');
fprintf('simplex:	%g	%g evals	@	',fs,cs); fprintf('%g,	',xs); 
fprintf('error:	'); fprintf('%g	',exct-xs'); fprintf('\n');
fprintf('minfun:	%g	%g evals	@	',val,cnt); fprintf('%g	',x); 
fprintf('error:	'); fprintf('%g	',exct-x');
fprintf('\n');
f = figure; set(f,'color','w'); hold on; grid on;
xlabel('X1'); ylabel('X2'); colorbar; title(char('Beale',char(fun)))
xv = linspace(min([ic(1) x(1) xs(1) xm(1)]),max([ic(1) x(1) xs(1) xm(1)]),np)';
yv = linspace(min([ic(2) x(2) xs(2) xm(2)]),max([ic(2) x(2) xs(2) xm(2)]),np);
[xx, yy] = meshgrid(xv,yv);
zz = xx;
for ii = 1:np^2
	zz(ii) = feval(fun,[xx(ii) yy(ii)]);
end
imagesc(xv,yv,zz); axis image
temp = plot(xs(1),xs(2),'k*'); set(temp,'markersize',10)
temp = plot(x(1),x(2),'ko'); set(temp,'markersize',10)
temp = plot(xm(1),xm(2),'kd'); set(temp,'markersize',10)
legend('Simplex','Min','Fminsearch')
pause(pse); fprintf('\n\n');


% booth
fn = 'booth';
fun = @(x)((x(1)+2*x(2)-7)^2 + (2*x(1) + x(2) - 5)^2);
exct = [1 3];
ic = [-5 2];
% straight simplex
[xs, fs, cs] = simplexmin(fun,ic,[],[],[],[],[],[],[],dbg);
% simplex function min
[x, val, cnt] = minfun(fun,ic,[],[],[],[],[],[],[],dbg);
% fmin function
[xm, fm, ~, out] = fminsearch(fun,ic);
% fm = feval(fun,xm);
fprintf('\n');
fprintf('%s	(%s)	',fn,func2str(fun));
fprintf('Exact soln:	'); fprintf('%g	',exct); fprintf('\n');
fprintf('ic:	'); fprintf('%g	',ic); fprintf('\n');
fprintf('fmin:	%g	%g evals	@	',fm,out.funcCount); fprintf('%g	',xm); 
fprintf('error:	'); fprintf('%g	',exct-xm); fprintf('\n');
fprintf('simplex:	%g	%g evals	@	',fs,cs); fprintf('%g,	',xs); 
fprintf('error:	'); fprintf('%g	',exct-xs'); fprintf('\n');
fprintf('minfun:	%g	%g evals	@	',val,cnt); fprintf('%g	',x); 
fprintf('error:	'); fprintf('%g	',exct-x');
fprintf('\n');
f = figure; set(f,'color','w'); hold on; grid on;
xlabel('X1'); ylabel('X2'); colorbar; title(char('Booth',char(fun)))
xv = linspace(min([ic(1) x(1) xs(1) xm(1)]),max([ic(1) x(1) xs(1) xm(1)]),np)';
yv = linspace(min([ic(2) x(2) xs(2) xm(2)]),max([ic(2) x(2) xs(2) xm(2)]),np);
[xx, yy] = meshgrid(xv,yv);
zz = xx;
for ii = 1:np^2
	zz(ii) = feval(fun,[xx(ii) yy(ii)]);
end
imagesc(xv,yv,zz); axis image
temp = plot(xs(1),xs(2),'k*'); set(temp,'markersize',10)
temp = plot(x(1),x(2),'ko'); set(temp,'markersize',10)
temp = plot(xm(1),xm(2),'kd'); set(temp,'markersize',10)
legend('Simplex','Min','Fminsearch')
pause(pse); fprintf('\n\n');


% bukin
fn = 'bukin';
fun = @(x)(100*sqrt(abs(x(2)-0.01*x(1)^2)) + 0.01*abs(x(1) + 10));
exct = [-10 1];
ic = [-6 -5];
% straight simplex
[xs, fs, cs] = simplexmin(fun,ic,[],[],[],[],[],[],[],dbg);
% simplex function min
[x, val, cnt] = minfun(fun,ic,[],[],[],[],[],[],[],dbg);
% fmin function
[xm, fm, ~, out] = fminsearch(fun,ic);
% fm = feval(fun,xm);
fprintf('\n');
fprintf('%s	(%s)	',fn,func2str(fun));
fprintf('Exact soln:	'); fprintf('%g	',exct); fprintf('\n');
fprintf('ic:	'); fprintf('%g	',ic); fprintf('\n');
fprintf('fmin:	%g	%g evals	@	',fm,out.funcCount); fprintf('%g	',xm); 
fprintf('error:	'); fprintf('%g	',exct-xm); fprintf('\n');
fprintf('simplex:	%g	%g evals	@	',fs,cs); fprintf('%g,	',xs); 
fprintf('error:	'); fprintf('%g	',exct-xs'); fprintf('\n');
fprintf('minfun:	%g	%g evals	@	',val,cnt); fprintf('%g	',x); 
fprintf('error:	'); fprintf('%g	',exct-x');
fprintf('\n');
f = figure; set(f,'color','w'); hold on; grid on;
xlabel('X1'); ylabel('X2'); colorbar; title(char('Bukin',char(fun)))
xv = linspace(min([ic(1) x(1) xs(1) xm(1)]),max([ic(1) x(1) xs(1) xm(1)]),np)';
yv = linspace(min([ic(2) x(2) xs(2) xm(2)]),max([ic(2) x(2) xs(2) xm(2)]),np);
[xx, yy] = meshgrid(xv,yv);
zz = xx;
for ii = 1:np^2
	zz(ii) = feval(fun,[xx(ii) yy(ii)]);
end
imagesc(xv,yv,zz); axis image
temp = plot(xs(1),xs(2),'k*'); set(temp,'markersize',10)
temp = plot(x(1),x(2),'ko'); set(temp,'markersize',10)
temp = plot(xm(1),xm(2),'kd'); set(temp,'markersize',10)
legend('Simplex','Min','Fminsearch')
pause(pse); fprintf('\n\n');



% 'matyas'
fn = 'matyas';
fun = @(x)(0.26*(x(1)^2 + x(2)^2) - 0.48*x(1)*x(2));
exct = [0 0];
ic = [-3 5];
% straight simplex
[xs, fs, cs] = simplexmin(fun,ic,[],[],[],[],[],[],[],dbg);
% simplex function min
[x, val, cnt] = minfun(fun,ic,[],[],[],[],[],[],[],dbg);
% fmin function
[xm, fm, ~, out] = fminsearch(fun,ic);
% fm = feval(fun,xm);
fprintf('\n');
fprintf('%s	(%s)	',fn,func2str(fun));
fprintf('Exact soln:	'); fprintf('%g	',exct); fprintf('\n');
fprintf('ic:	'); fprintf('%g	',ic); fprintf('\n');
fprintf('fmin:	%g	%g evals	@	',fm,out.funcCount); fprintf('%g	',xm); 
fprintf('error:	'); fprintf('%g	',exct-xm); fprintf('\n');
fprintf('simplex:	%g	%g evals	@	',fs,cs); fprintf('%g,	',xs); 
fprintf('error:	'); fprintf('%g	',exct-xs'); fprintf('\n');
fprintf('minfun:	%g	%g evals	@	',val,cnt); fprintf('%g	',x); 
fprintf('error:	'); fprintf('%g	',exct-x');
fprintf('\n');
f = figure; set(f,'color','w'); hold on; grid on;
xlabel('X1'); ylabel('X2'); colorbar; title(char('Matyas',char(fun)))
xv = linspace(min([ic(1) x(1) xs(1) xm(1)]),max([ic(1) x(1) xs(1) xm(1)]),np)';
yv = linspace(min([ic(2) x(2) xs(2) xm(2)]),max([ic(2) x(2) xs(2) xm(2)]),np);
[xx, yy] = meshgrid(xv,yv);
zz = xx;
for ii = 1:np^2
	zz(ii) = feval(fun,[xx(ii) yy(ii)]);
end
imagesc(xv,yv,zz); axis image
temp = plot(xs(1),xs(2),'k*'); set(temp,'markersize',10)
temp = plot(x(1),x(2),'ko'); set(temp,'markersize',10)
temp = plot(xm(1),xm(2),'kd'); set(temp,'markersize',10)
legend('Simplex','Min','Fminsearch')
pause(pse); fprintf('\n\n');


% easom
fn = 'easom';
fun = @(x)(-cos(x(1))*cos(x(2))*exp(-((x(1)-pi)^2 +(x(2)-pi)^2)));
exct = [pi pi];
ic = [1 1];
% straight simplex
[xs, fs, cs] = simplexmin(fun,ic,[],[],[],[],[],[],[],dbg);
% simplex function min
[x, val, cnt] = minfun(fun,ic,[],[],[],[],[],[],[],dbg);
% fmin function
[xm, fm, ~, out] = fminsearch(fun,ic);
% fm = feval(fun,xm);
fprintf('\n');
fprintf('%s	(%s)	',fn,func2str(fun));
fprintf('Exact soln:	'); fprintf('%g	',exct); fprintf('\n');
fprintf('ic:	'); fprintf('%g	',ic); fprintf('\n');
fprintf('fmin:	%g	%g evals	@	',fm,out.funcCount); fprintf('%g	',xm); 
fprintf('error:	'); fprintf('%g	',exct-xm); fprintf('\n');
fprintf('simplex:	%g	%g evals	@	',fs,cs); fprintf('%g,	',xs); 
fprintf('error:	'); fprintf('%g	',exct-xs'); fprintf('\n');
fprintf('minfun:	%g	%g evals	@	',val,cnt); fprintf('%g	',x); 
fprintf('error:	'); fprintf('%g	',exct-x');
fprintf('\n');
f = figure; set(f,'color','w'); hold on; grid on;
xlabel('X1'); ylabel('X2'); colorbar; title(char('Easom',char(fun)))
xv = linspace(min([ic(1) x(1) xs(1) xm(1)]),max([ic(1) x(1) xs(1) xm(1)]),np)';
yv = linspace(min([ic(2) x(2) xs(2) xm(2)]),max([ic(2) x(2) xs(2) xm(2)]),np);
[xx, yy] = meshgrid(xv,yv);
zz = xx;
for ii = 1:np^2
	zz(ii) = feval(fun,[xx(ii) yy(ii)]);
end
imagesc(xv,yv,zz); axis image
temp = plot(xs(1),xs(2),'k*'); set(temp,'markersize',10)
temp = plot(x(1),x(2),'ko'); set(temp,'markersize',10)
temp = plot(xm(1),xm(2),'kd'); set(temp,'markersize',10)
legend('Simplex','Min','Fminsearch')
pause(pse); fprintf('\n\n');



% mccormick
fn = 'mccormick';
fun = @(x)(sin(x(1)+x(2)) + (x(1)-x(2))^2 - 1.5*x(1) + 2.5*x(2) + 2.9133);
exct = [-0.54719 -1.54719];
ic = [1 -1];
% straight simplex
[xs, fs, cs] = simplexmin(fun,ic,[],[],[],[],[],[],[],dbg);
% simplex function min
[x, val, cnt] = minfun(fun,ic,[],[],[],[],[],[],[],dbg);
% fmin function
[xm, fm, ~, out] = fminsearch(fun,ic);
% fm = feval(fun,xm);
fprintf('\n');
fprintf('%s	(%s)	',fn,func2str(fun));
fprintf('Exact soln:	'); fprintf('%g	',exct); fprintf('\n');
fprintf('ic:	'); fprintf('%g	',ic); fprintf('\n');
fprintf('fmin:	%g	%g evals	@	',fm,out.funcCount); fprintf('%g	',xm); 
fprintf('error:	'); fprintf('%g	',exct-xm); fprintf('\n');
fprintf('simplex:	%g	%g evals	@	',fs,cs); fprintf('%g,	',xs); 
fprintf('error:	'); fprintf('%g	',exct-xs'); fprintf('\n');
fprintf('minfun:	%g	%g evals	@	',val,cnt); fprintf('%g	',x); 
fprintf('error:	'); fprintf('%g	',exct-x');
fprintf('\n');
f = figure; set(f,'color','w'); hold on; grid on;
xlabel('X1'); ylabel('X2'); colorbar; title(char('Mccormick',char(fun)))
xv = linspace(min([ic(1) x(1) xs(1) xm(1)]),max([ic(1) x(1) xs(1) xm(1)]),np)';
yv = linspace(min([ic(2) x(2) xs(2) xm(2)]),max([ic(2) x(2) xs(2) xm(2)]),np);
[xx, yy] = meshgrid(xv,yv);
zz = xx;
for ii = 1:np^2
	zz(ii) = feval(fun,[xx(ii) yy(ii)]);
end
imagesc(xv,yv,zz); axis image
temp = plot(xs(1),xs(2),'k*'); set(temp,'markersize',10)
temp = plot(x(1),x(2),'ko'); set(temp,'markersize',10)
temp = plot(xm(1),xm(2),'kd'); set(temp,'markersize',10)
legend('Simplex','Min','Fminsearch')
pause(pse); fprintf('\n\n');



% schaffer
fn = 'Schaffer';
fun = @(x)(0.5 + ((sin(x(1)^2 - x(2)^2))^2 - 0.5) / ((1 + 0.001*(x(1)^2 + x(2)^2))^2));
exct = [0 0];
ic = [1 -2];
% straight simplex
[xs, fs, cs] = simplexmin(fun,ic,[],[],[],[],[],[],[],dbg);
% simplex function min
[x, val, cnt] = minfun(fun,ic,[],[],[],[],[],[],[],dbg);
% fmin function
[xm, fm, ~, out] = fminsearch(fun,ic);
% fm = feval(fun,xm);
fprintf('\n');
fprintf('%s	(%s)	',fn,func2str(fun));
fprintf('Exact soln:	'); fprintf('%g	',exct); fprintf('\n');
fprintf('ic:	'); fprintf('%g	',ic); fprintf('\n');
fprintf('fmin:	%g	%g evals	@	',fm,out.funcCount); fprintf('%g	',xm); 
fprintf('error:	'); fprintf('%g	',exct-xm); fprintf('\n');
fprintf('simplex:	%g	%g evals	@	',fs,cs); fprintf('%g,	',xs); 
fprintf('error:	'); fprintf('%g	',exct-xs'); fprintf('\n');
fprintf('minfun:	%g	%g evals	@	',val,cnt); fprintf('%g	',x); 
fprintf('error:	'); fprintf('%g	',exct-x');
fprintf('\n');
f = figure; set(f,'color','w'); hold on; grid on;
xlabel('X1'); ylabel('X2'); colorbar; title(char('Schaffer',char(fun)))
xv = linspace(min([ic(1) x(1) xs(1) xm(1)]),max([ic(1) x(1) xs(1) xm(1)]),np)';
yv = linspace(min([ic(2) x(2) xs(2) xm(2)]),max([ic(2) x(2) xs(2) xm(2)]),np);
[xx, yy] = meshgrid(xv,yv);
zz = xx;
for ii = 1:np^2
	zz(ii) = feval(fun,[xx(ii) yy(ii)]);
end
imagesc(xv,yv,zz); axis image
temp = plot(xs(1),xs(2),'k*'); set(temp,'markersize',10)
temp = plot(x(1),x(2),'ko'); set(temp,'markersize',10)
temp = plot(xm(1),xm(2),'kd'); set(temp,'markersize',10)
legend('Simplex','Min','Fminsearch')
pause(pse); fprintf('\n\n');


% Styblinski
fn = 'Styblinski';
fun = @(x)((x(1)^4-16*x(1)^2+5*x(1)+x(2)^4-16*x(2)^2+5*x(2))/2);
exct = [-2.903534 -2.903534];
ic = [1 2];
% straight simplex
[xs, fs, cs] = simplexmin(fun,ic,[],[],[],[],[],[],[],dbg);
% simplex function min
[x, val, cnt] = minfun(fun,ic,[],[],[],[],[],[],[],dbg);
% fmin function
[xm, fm, ~, out] = fminsearch(fun,ic);
% fm = feval(fun,xm);
fprintf('\n');
fprintf('%s	(%s)	',fn,func2str(fun));
fprintf('Exact soln:	'); fprintf('%g	',exct); fprintf('\n');
fprintf('ic:	'); fprintf('%g	',ic); fprintf('\n');
fprintf('fmin:	%g	%g evals	@	',fm,out.funcCount); fprintf('%g	',xm); 
fprintf('error:	'); fprintf('%g	',exct-xm); fprintf('\n');
fprintf('simplex:	%g	%g evals	@	',fs,cs); fprintf('%g,	',xs); 
fprintf('error:	'); fprintf('%g	',exct-xs'); fprintf('\n');
fprintf('minfun:	%g	%g evals	@	',val,cnt); fprintf('%g	',x); 
fprintf('error:	'); fprintf('%g	',exct-x');
fprintf('\n');
f = figure; set(f,'color','w'); hold on; grid on;
xlabel('X1'); ylabel('X2'); colorbar; title(char('Styblinski',char(fun)))
xv = linspace(min([ic(1) x(1) xs(1) xm(1)]),max([ic(1) x(1) xs(1) xm(1)]),np)';
yv = linspace(min([ic(2) x(2) xs(2) xm(2)]),max([ic(2) x(2) xs(2) xm(2)]),np);
[xx, yy] = meshgrid(xv,yv);
zz = xx;
for ii = 1:np^2
	zz(ii) = feval(fun,[xx(ii) yy(ii)]);
end
imagesc(xv,yv,zz); axis image
temp = plot(xs(1),xs(2),'k*'); set(temp,'markersize',10)
temp = plot(x(1),x(2),'ko'); set(temp,'markersize',10)
temp = plot(xm(1),xm(2),'kd'); set(temp,'markersize',10)
legend('Simplex','Min','Fminsearch')
pause(pse); fprintf('\n\n');



%% bounded
fprintf('Constrained test - https://en.wikipedia.org/wiki/Test_functions_for_optimization\n\n\n')



% allocation
out.funcCount = 0;
fm = Inf;


% rosenbrock 1
fn = 'rosenbrock 1';
fun = @(x)((1 - x(1))^2 + 100*(x(2) - x(1)^2)^2);
Ane = {@(x)((x(1)-1)^3 - x(2) + 1); @(x)(x(1) + x(2) - 2)};
Bne = [0; 0];
Ae = [];
Be = [];
exct = [1 1];
ic = [0.25 2.0];
ub = [1.5 2.5]; lb = [-1.5 -0.5];
% straight simplex
[xs, fs, cs] = simplexmin(fun,ic,[],Ane,Bne,Ae,Be,lb,ub,dbg);
% simplex function min
[x, val, cnt] = minfun(fun,ic,[],Ane,Bne,Ae,Be,lb,ub,dbg);
% fmin function
options = optimoptions('fmincon','Display','none','Algorithm','sqp');
xm = NaN*x';
[xm, fm , ~, out] = fmincon(fun,ic,[],[],[],[],lb',ub',@rosenbrook_const1,options);
% fm = feval(fun,xm);
fprintf('\n');
fprintf('%s	(%s)\n',fn,func2str(fun));
fprintf('ub:	'); fprintf('%g	',ub); fprintf('\n');
fprintf('lb:	'); fprintf('%g	',lb); fprintf('\n');
fprintf('Exact soln:	'); fprintf('%g	',exct); fprintf('\n');
fprintf('ic:	'); fprintf('%g	',ic); fprintf('\n');
fprintf('fmincon:	%g	%g evals	@	',fm,out.funcCount); fprintf('%g	',xm); 
fprintf('error:	'); fprintf('%g	',exct-xm); fprintf('\n');
fprintf('simplex:	%g	%g evals	@	',fs,cs); fprintf('%g,	',xs); 
fprintf('error:	'); fprintf('%g	',exct-xs'); fprintf('\n');
fprintf('minfun:	%g	%g evals	@	',val,cnt); fprintf('%g	',x); 
fprintf('error:	'); fprintf('%g	',exct-x');
fprintf('\n');

f = figure; set(f,'color','w'); hold on; grid on;
xlabel('X1'); ylabel('X2'); colorbar; title(char('rosenbrock 1',char(fun)))
xv = linspace(lb(1),ub(1),np)'; yv = linspace(lb(2),ub(2),np);
[xx, yy] = meshgrid(xv,yv);
zz = xx;
for ii = 1:np^2
	zz(ii) = feval(fun,[xx(ii) yy(ii)]);
end
imagesc(xv,yv,zz); axis image
temp = plot(xs(1),xs(2),'k*'); set(temp,'markersize',10)
temp = plot(x(1),x(2),'ko'); set(temp,'markersize',10)
temp = plot(xm(1),xm(2),'kd'); set(temp,'markersize',10)
legend('Simplex','Min','Fmincon')
pause(pse); fprintf('\n\n');



% rosenbrock 2
fn = 'rosenbrock 2';
fun = @(x)((1 - x(1))^2 + 100*(x(2) - x(1)^2)^2);
Ane = {@(x)((x(1))^2 + x(2) - 2)};
Bne = [0]; %#ok<*NBRAK>
Ae = [];
Be = [];
exct = [1 1];
ic = [-1.25 -1.25];
ub = [1.5 1.5]; lb = [-1.5 -1.5];
% straight simplex
[xs, fs, cs] = simplexmin(fun,ic,[],Ane,Bne,Ae,Be,lb,ub,dbg);
% simplex function min
[x, val, cnt] = minfun(fun,ic,[],Ane,Bne,Ae,Be,lb,ub,dbg);
% fmin function
options = optimoptions('fmincon','Display','none','Algorithm','sqp');
xm = NaN*x';
[xm, fm , ~, out] = fmincon(fun,ic,[],[],[],[],lb',ub',@rosenbrook_const2,options);
% fm = feval(fun,xm);
fprintf('\n');
fprintf('%s	(%s)\n',fn,func2str(fun));
fprintf('ub:	'); fprintf('%g	',ub); fprintf('\n');
fprintf('lb:	'); fprintf('%g	',lb); fprintf('\n');
fprintf('Exact soln:	'); fprintf('%g	',exct); fprintf('\n');
fprintf('ic:	'); fprintf('%g	',ic); fprintf('\n');
fprintf('fmincon:	%g	%g evals	@	',fm,out.funcCount); fprintf('%g	',xm); 
fprintf('error:	'); fprintf('%g	',exct-xm); fprintf('\n');
fprintf('simplex:	%g	%g evals	@	',fs,cs); fprintf('%g,	',xs); 
fprintf('error:	'); fprintf('%g	',exct-xs'); fprintf('\n');
fprintf('minfun:	%g	%g evals	@	',val,cnt); fprintf('%g	',x); 
fprintf('error:	'); fprintf('%g	',exct-x');
fprintf('\n');
f = figure; set(f,'color','w'); hold on; grid on;
xlabel('X1'); ylabel('X2'); colorbar; title(char('rosenbrock 2',char(fun)))
xv = linspace(lb(1),ub(1),np)'; yv = linspace(lb(2),ub(2),np);
[xx, yy] = meshgrid(xv,yv);
zz = xx;
for ii = 1:np^2
	zz(ii) = feval(fun,[xx(ii) yy(ii)]);
end
imagesc(xv,yv,zz); axis image
temp = plot(xs(1),xs(2),'k*'); set(temp,'markersize',10)
temp = plot(x(1),x(2),'ko'); set(temp,'markersize',10)
temp = plot(xm(1),xm(2),'kd'); set(temp,'markersize',10)
legend('Simplex','Min','Fmincon')
pause(pse); fprintf('\n\n');



% mishra
fn = 'mishra';
fun = @(x)(sin(x(2)) * exp((1 - cos(x(1)))^2) + cos(x(1)) * exp((1 - cos(x(2)))^2) + (x(1) - x(2))^2);
Ane = {@(x)((x(1) + 5)^2 + (x(2) + 5)^2 - 25)};
Bne = [0];
Ae = [];
Be = [];
exct = [-3.1302468 -1.5821422];
ic = [-9.25 -1.25];
ub = [0 0]; lb = [-10 -6.5];
% straight simplex
[xs, fs, cs] = simplexmin(fun,ic,[],Ane,Bne,Ae,Be,lb,ub,dbg);
% simplex function min
[x, val, cnt] = minfun(fun,ic,[],Ane,Bne,Ae,Be,lb,ub,dbg);
% fmin function
options = optimoptions('fmincon','Display','none','Algorithm','sqp');
xm = NaN*x';
[xm, fm , ~, out] = fmincon(fun,ic,[],[],[],[],lb',ub',@mishra_const,options);
% fm = feval(fun,xm);
fprintf('\n');
fprintf('%s	(%s)\n',fn,func2str(fun));
fprintf('ub:	'); fprintf('%g	',ub); fprintf('\n');
fprintf('lb:	'); fprintf('%g	',lb); fprintf('\n');
fprintf('Exact soln:	'); fprintf('%g	',exct); fprintf('\n');
fprintf('ic:	'); fprintf('%g	',ic); fprintf('\n');
fprintf('fmincon:	%g	%g evals	@	',fm,out.funcCount); fprintf('%g	',xm); 
fprintf('error:	'); fprintf('%g	',exct-xm); fprintf('\n');
fprintf('simplex:	%g	%g evals	@	',fs,cs); fprintf('%g,	',xs); 
fprintf('error:	'); fprintf('%g	',exct-xs'); fprintf('\n');
fprintf('minfun:	%g	%g evals	@	',val,cnt); fprintf('%g	',x); 
fprintf('error:	'); fprintf('%g	',exct-x');
fprintf('\n');
f = figure; set(f,'color','w'); hold on; grid on;
xlabel('X1'); ylabel('X2'); colorbar; title(char('mishra',char(fun)))
xv = linspace(lb(1),ub(1),np)'; yv = linspace(lb(2),ub(2),np);
[xx, yy] = meshgrid(xv,yv);
zz = xx;
for ii = 1:np^2
	zz(ii) = feval(fun,[xx(ii) yy(ii)]);
end
imagesc(xv,yv,zz); axis image
temp = plot(xs(1),xs(2),'k*'); set(temp,'markersize',10)
temp = plot(x(1),x(2),'ko'); set(temp,'markersize',10)
temp = plot(xm(1),xm(2),'kd'); set(temp,'markersize',10)
legend('Simplex','Min','Fmincon')
pause(pse); fprintf('\n\n');



% townsend
fn = 'townsend';
fun = @(x)(-(cos((x(1) - 0.1) * x(2)))^2 - x(1) * sin(3*x(1) + x(2)));
Ane = {@(x)(x(1)^2 + x(2)^2 - ((2*cos((atan2(x(1),x(2)))) - ...
	cos(2 * (atan2(x(1),x(2))))/2 - cos(3 * (atan2(x(1),x(2))))/4 - ...
	cos(4 * (atan2(x(1),x(2))))/8)^2 + (2 * sin((atan2(x(1),x(2)))))^2))};
Bne = [0];
Ae = [];
Be = [];
exct = [2.0052938 1.1944509];
ic = [1.0 0.5];
ub = [2.5 1.75]; lb = [-2.25 -2.5];
% straight simplex
[xs, fs, cs] = simplexmin(fun,ic,[],Ane,Bne,Ae,Be,lb,ub,dbg);
% simplex function min
[x, val, cnt] = minfun(fun,ic,[],Ane,Bne,Ae,Be,lb,ub,dbg);
% fmin function
options = optimoptions('fmincon','Display','none','Algorithm','sqp');
xm = NaN*x';
[xm, fm , ~, out] = fmincon(fun,ic,[],[],[],[],lb',ub',@townsend_const,options);
% fm = feval(fun,xm);
fprintf('\n');
fprintf('%s	(%s)\n',fn,func2str(fun));
fprintf('ub:	'); fprintf('%g	',ub); fprintf('\n');
fprintf('lb:	'); fprintf('%g	',lb); fprintf('\n');
fprintf('Exact soln:	'); fprintf('%g	',exct); fprintf('\n');
fprintf('ic:	'); fprintf('%g	',ic); fprintf('\n');
fprintf('fmincon:	%g	%g evals	@	',fm,out.funcCount); fprintf('%g	',xm); 
fprintf('error:	'); fprintf('%g	',exct-xm); fprintf('\n');
fprintf('simplex:	%g	%g evals	@	',fs,cs); fprintf('%g,	',xs); 
fprintf('error:	'); fprintf('%g	',exct-xs'); fprintf('\n');
fprintf('minfun:	%g	%g evals	@	',val,cnt); fprintf('%g	',x); 
fprintf('error:	'); fprintf('%g	',exct-x');
fprintf('\n');
f = figure; set(f,'color','w'); hold on; grid on;
xlabel('X1'); ylabel('X2'); colorbar; title(char('townsend',char(fun)))
xv = linspace(lb(1),ub(1),np)'; yv = linspace(lb(2),ub(2),np);
[xx, yy] = meshgrid(xv,yv);
zz = xx;
for ii = 1:np^2
	zz(ii) = feval(fun,[xx(ii) yy(ii)]);
end
imagesc(xv,yv,zz); axis image
temp = plot(xs(1),xs(2),'k*'); set(temp,'markersize',10)
temp = plot(x(1),x(2),'ko'); set(temp,'markersize',10)
temp = plot(xm(1),xm(2),'kd'); set(temp,'markersize',10)
legend('Simplex','Min','Fmincon')
pause(pse); fprintf('\n\n');




fn = 'townsend';
fun = @(x)(-(cos((x(1) - 0.1) * x(2)))^2 - x(1) * sin(3*x(1) + x(2)));
Ane = {@(x)(x(1)^2 + x(2)^2 - ((2*cos((atan2(x(1),x(2)))) - ...
	cos(2 * (atan2(x(1),x(2))))/2 - cos(3 * (atan2(x(1),x(2))))/4 - ...
	cos(4 * (atan2(x(1),x(2))))/8)^2 + (2 * sin((atan2(x(1),x(2)))))^2))};
Bne = [0];
Ae = [];
Be = [];
exct = [2.0052938 1.1944509];
ic = [-1 0.5];
ub = [2.5 1.75]; lb = [-2.25 -2.5];
% straight simplex
[xs, fs, cs] = simplexmin(fun,ic,[],Ane,Bne,Ae,Be,lb,ub,dbg);
% simplex function min
[x, val, cnt] = minfun(fun,ic,[],Ane,Bne,Ae,Be,lb,ub,dbg);
% fmin function
options = optimoptions('fmincon','Display','none','Algorithm','sqp');
xm = NaN*x';
[xm, fm , ~, out] = fmincon(fun,ic,[],[],[],[],lb',ub',@townsend_const,options);
% fm = feval(fun,xm);
fprintf('\n');
fprintf('%s	(%s)\n',fn,func2str(fun));
fprintf('ub:	'); fprintf('%g	',ub); fprintf('\n');
fprintf('lb:	'); fprintf('%g	',lb); fprintf('\n');
fprintf('Exact soln:	'); fprintf('%g	',exct); fprintf('\n');
fprintf('ic:	'); fprintf('%g	',ic); fprintf('\n');
fprintf('fmincon:	%g	%g evals	@	',fm,out.funcCount); fprintf('%g	',xm); 
fprintf('error:	'); fprintf('%g	',exct-xm); fprintf('\n');
fprintf('simplex:	%g	%g evals	@	',fs,cs); fprintf('%g,	',xs); 
fprintf('error:	'); fprintf('%g	',exct-xs'); fprintf('\n');
fprintf('minfun:	%g	%g evals	@	',val,cnt); fprintf('%g	',x); 
fprintf('error:	'); fprintf('%g	',exct-x');
fprintf('\n');
f = figure; set(f,'color','w'); hold on; grid on;
xlabel('X1'); ylabel('X2'); colorbar; title(char('townsend',char(fun)))
xv = linspace(lb(1),ub(1),np)'; yv = linspace(lb(2),ub(2),np);
[xx, yy] = meshgrid(xv,yv);
zz = xx;
for ii = 1:np^2
	zz(ii) = feval(fun,[xx(ii) yy(ii)]);
end
imagesc(xv,yv,zz); axis image
temp = plot(xs(1),xs(2),'k*'); set(temp,'markersize',10)
temp = plot(x(1),x(2),'ko'); set(temp,'markersize',10)
temp = plot(xm(1),xm(2),'kd'); set(temp,'markersize',10)
legend('Simplex','Min','Fmincon')
pause(pse); fprintf('\n\n');




fn = 'townsend';
fun = @(x)(-(cos((x(1) - 0.1) * x(2)))^2 - x(1) * sin(3*x(1) + x(2)));
Ane = {@(x)(x(1)^2 + x(2)^2 - ((2*cos((atan2(x(1),x(2)))) - ...
	cos(2 * (atan2(x(1),x(2))))/2 - cos(3 * (atan2(x(1),x(2))))/4 - ...
	cos(4 * (atan2(x(1),x(2))))/8)^2 + (2 * sin((atan2(x(1),x(2)))))^2))};
Bne = [0];
Ae = [];
Be = [];
exct = [2.0052938 1.1944509];
ic = [-2.0 -1.5];
ub = [2.5 1.75]; lb = [-2.25 -2.5];
% straight simplex
[xs, fs, cs] = simplexmin(fun,ic,[],Ane,Bne,Ae,Be,lb,ub,dbg);
% simplex function min
[x, val, cnt] = minfun(fun,ic,[],Ane,Bne,Ae,Be,lb,ub,dbg);
% fmin function
options = optimoptions('fmincon','Display','none','Algorithm','sqp');
xm = NaN*x';
[xm, fm , ~, out] = fmincon(fun,ic,[],[],[],[],lb',ub',@townsend_const,options);
% fm = feval(fun,xm);
fprintf('\n');
fprintf('%s	(%s)\n',fn,func2str(fun));
fprintf('ub:	'); fprintf('%g	',ub); fprintf('\n');
fprintf('lb:	'); fprintf('%g	',lb); fprintf('\n');
fprintf('Exact soln:	'); fprintf('%g	',exct); fprintf('\n');
fprintf('ic:	'); fprintf('%g	',ic); fprintf('\n');
fprintf('fmincon:	%g	%g evals	@	',fm,out.funcCount); fprintf('%g	',xm); 
fprintf('error:	'); fprintf('%g	',exct-xm); fprintf('\n');
fprintf('simplex:	%g	%g evals	@	',fs,cs); fprintf('%g,	',xs); 
fprintf('error:	'); fprintf('%g	',exct-xs'); fprintf('\n');
fprintf('minfun:	%g	%g evals	@	',val,cnt); fprintf('%g	',x); 
fprintf('error:	'); fprintf('%g	',exct-x');
fprintf('\n');
f = figure; set(f,'color','w'); hold on; grid on;
xlabel('X1'); ylabel('X2'); colorbar; title(char('townsend',char(fun)))
xv = linspace(lb(1),ub(1),np)'; yv = linspace(lb(2),ub(2),np);
[xx, yy] = meshgrid(xv,yv);
zz = xx;
for ii = 1:np^2
	zz(ii) = feval(fun,[xx(ii) yy(ii)]);
end
imagesc(xv,yv,zz); axis image
temp = plot(xs(1),xs(2),'k*'); set(temp,'markersize',10)
temp = plot(x(1),x(2),'ko'); set(temp,'markersize',10)
temp = plot(xm(1),xm(2),'kd'); set(temp,'markersize',10)
legend('Simplex','Min','Fmincon')
pause(pse); fprintf('\n\n');



% simionescu
fn = 'simionescu';
fun = @(x)(0.1 * x(1) * x(2));
Ane = {@(x)(x(1)^2 + x(2)^2 - (1 + 0.2 * cos(8 * atan(x(1)/x(2))))^2)};
Bne = [0];
Ae = [];
Be = [];
exct = [0.84852813 -0.84852813];
ic = [0.5 -1.0];
ub = [1.25 1.25]; lb = [-1.25 -1.25];
% straight simplex
[xs, fs, cs] = simplexmin(fun,ic,[],Ane,Bne,Ae,Be,lb,ub,dbg);
% simplex function min
[x, val, cnt] = minfun(fun,ic,[],Ane,Bne,Ae,Be,lb,ub,dbg);
% fmin function
options = optimoptions('fmincon','Display','none','Algorithm','sqp');
xm = NaN*x';
[xm, fm , ~, out] = fmincon(fun,ic,[],[],[],[],lb',ub',@simi_const,options);
% fm = feval(fun,xm);
fprintf('\n');
fprintf('%s	(%s)\n',fn,func2str(fun));
fprintf('ub:	'); fprintf('%g	',ub); fprintf('\n');
fprintf('lb:	'); fprintf('%g	',lb); fprintf('\n');
fprintf('Exact soln:	'); fprintf('%g	',exct); fprintf('\n');
fprintf('ic:	'); fprintf('%g	',ic); fprintf('\n');
fprintf('fmincon:	%g	%g evals	@	',fm,out.funcCount); fprintf('%g	',xm); 
fprintf('error:	'); fprintf('%g	',exct-xm); fprintf('\n');
fprintf('simplex:	%g	%g evals	@	',fs,cs); fprintf('%g,	',xs); 
fprintf('error:	'); fprintf('%g	',exct-xs'); fprintf('\n');
fprintf('minfun:	%g	%g evals	@	',val,cnt); fprintf('%g	',x); 
fprintf('error:	'); fprintf('%g	',exct-x');
fprintf('\n');
f = figure; set(f,'color','w'); hold on; grid on;
xlabel('X1'); ylabel('X2'); colorbar; title(char('simionescu',char(fun)))
xv = linspace(lb(1),ub(1),np)'; yv = linspace(lb(2),ub(2),np);
[xx, yy] = meshgrid(xv,yv);
zz = xx;
for ii = 1:np^2
	zz(ii) = feval(fun,[xx(ii) yy(ii)]);
end
imagesc(xv,yv,zz); axis image
temp = plot(xs(1),xs(2),'k*'); set(temp,'markersize',10)
temp = plot(x(1),x(2),'ko'); set(temp,'markersize',10)
temp = plot(xm(1),xm(2),'kd'); set(temp,'markersize',10)
legend('Simplex','Min','Fmincon')
pause(pse); fprintf('\n\n');



return



%% non linear constraints
function [c,ceq] = rosenbrook_const1(x) %#ok<*DEFNU>
c(1) = ((x(1)-1)^3 - x(2) + 1);
c(2) = (x(1) + x(2)-2);
ceq = [ ];

function [c,ceq] = rosenbrook_const2(x)
c(1) = (x(1)^2 + x(2)^2 - 2);
ceq = [ ];

function [c,ceq] = mishra_const(x)
c(1) = (x(1) + 5)^2 + (x(2) + 5)^2 - 25;
ceq = [ ];

function [c,ceq] = townsend_const(x)
c(1) = (x(1)^2 + x(2)^2 - ((2*cos((atan2(x(1),x(2)))) - ...
	cos(2 * (atan2(x(1),x(2))))/2 - cos(3 * (atan2(x(1),x(2))))/4 - ...
	cos(4 * (atan2(x(1),x(2))))/8)^2 + (2 * sin((atan2(x(1),x(2)))))^2));
ceq = [ ];

function [c,ceq] = simi_const(x)
c(1) = x(1)^2 + x(2)^2 - (1 + 0.2 * cos(8 * atan(x(1)/x(2))))^2;
ceq = [ ];
