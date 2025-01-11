load('system_adiff.mat'); 
[number_xy,~] = size(xy);
u_exact = zeros(number_xy,1);
for i = 1:number_xy
    a1 = xy(i,1);
    b1 = xy(i,2);
    u_exact(i,1) = cos(4*pi*a1)*cos(4*pi*b1); 
end

u_div = u_exact-x_gal;
u_error = u_div'*M*u_div;
u_error_l2_norm = sqrt(u_error);