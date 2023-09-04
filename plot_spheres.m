exact = importdata('exact_grid');
theta = exact(1,2:end);
phi   = exact(2:end,1);

data = exact(2:end,2:end);

[tgrid,pgrid] = meshgrid(theta,phi);
r = ones(size(tgrid));

[x,y,z] = sph2cart(tgrid,pgrid,r);

surf(x,y,z,data);shg


for i=1:3
 fname = ['loop_' num2str(i)];
exact = importdata(fname);
theta = exact(1,2:end);
phi   = exact(2:end,1);

data = exact(2:end,2:end);

[tgrid,pgrid] = meshgrid(theta,phi);
r = ones(size(tgrid));

[x,y,z] = sph2cart(tgrid,pgrid,r);
figure;
surf(x,y,z,data);shg
end

