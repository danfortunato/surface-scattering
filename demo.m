function demo

close all
figure(1), rundemo('laplace',   flat=true)
figure(2), rundemo('yukawa',    flat=true)
figure(3), rundemo('helmholtz', flat=true)
figure(4), rundemo('helmholtz', flat=false)

end

function rundemo(pde, opts)

arguments
    pde
    opts.flat = false
end

if ( opts.flat )
    d = 2;
    z = @(u,v) 0*u;
else
    d = 3;
    z = @(u,v) 0.1*(1-erf(25*(sqrt(u.^2+v.^2)-0.3)));
end

clf
switch lower(pde)
    case 'helmholtz'
        zk = 40;
        q = @(x,y,z) 2*exp(-160*((x+0.2).^2+(y-0.1).^2));
        f = @(x,y,z) 1.5*zk^2*exp(-200*((x-0.3).^2+(y+0.1).^2));
        inc = [3 1];
        tiledlayout(2, 4);
        solve_and_plot(z=z, zk=zk);
        solve_and_plot(z=z, zk=zk, q=q);
        solve_and_plot(z=z, zk=zk, f=f);
        solve_and_plot(z=z, zk=zk, q=q, f=f);
        solve_and_plot(z=z, zk=zk, inc=inc);
        solve_and_plot(z=z, zk=zk, inc=inc, q=q);
        solve_and_plot(z=z, zk=zk, inc=inc, f=f);
        solve_and_plot(z=z, zk=zk, inc=inc, q=q, f=f);
    case 'laplace'
        zk = 0;
        f = @(x,y,z) 2*exp(-200*((x-0.3).^2+(y+0.1).^2));
        tiledlayout(1, 2);
        solve_and_plot(z=z, zk=zk);
        solve_and_plot(z=z, zk=zk, f=f);
    case 'yukawa'
        zk = 8i;
        q = @(x,y,z) -600*exp(-230*((x-0.32).^2+(y+0.12).^2));
        f = @(x,y,z) 1.5*zk^2*exp(-200*((x-0.3).^2+(y+0.1).^2));
        tiledlayout(1, 4);
        solve_and_plot(z=z, zk=zk);
        solve_and_plot(z=z, zk=zk, q=q);
        solve_and_plot(z=z, zk=zk, f=f);
        solve_and_plot(z=z, zk=zk, q=q, f=f);
end

ax = findall(gcf, 'type', 'axes');
view(ax, d);
shg

end

function solve_and_plot(opts)

arguments
    opts.q = 0;
    opts.f = 0;
    opts.inc = 0;
    opts.zk = 0;
    opts.z = @(u,v) 0*v;
    opts.test = 1
end

n = 16;
nref = 3;
rect = [-0.5 0.5 -0.5 0.5];
L = surfacescatterer(n, z=opts.z, q=opts.q, zk=opts.zk, nref=nref, rect=rect);
u = L.solve(opts.f, inc=opts.inc);

% Evaluate in the exterior
m = 100;
padx = 0.2*diff(rect(1:2));
pady = 0.2*diff(rect(3:4));
[xx, yy] = meshgrid(linspace(rect(1)-padx, rect(2)+padx, m), ...
                    linspace(rect(3)-pady, rect(4)+pady, m));
ii = xx < rect(1) | xx > rect(2) | yy < rect(3) | yy > rect(4);
vv = nan(m);
vv(ii) = u.ext(xx(ii),yy(ii));

% Plot
nexttile
plot(real(u.int)), hold on
pcolor(xx, yy, real(vv)), hold off
shading interp
colormap turbo
camlight
colorbar

% Add a title
t = [];
if ( isnumeric(opts.q) && isscalar(opts.q) && opts.q == 0 )
    t = [t 'q=0, '];
else
    t = [t 'q≠0, '];
end
if ( isnumeric(opts.f) && isscalar(opts.f) && opts.f == 0 )
    t = [t 'f=0, '];
else
    t = [t 'f≠0, '];
end
if ( isnumeric(opts.inc) && isscalar(opts.inc) && opts.inc == 0 )
    t = [t 'inc=0'];
else
    t = [t 'inc≠0'];
end
title(t)

end
