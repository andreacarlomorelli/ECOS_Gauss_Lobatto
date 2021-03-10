function [x_cart] = transform_state_cyl_to_cart(x_cyl, th)

th = th';

r = x_cyl(:,1); w = x_cyl(:,2); vr = x_cyl(:,3); 
vth = x_cyl(:,4); vw = x_cyl(:,5);

[x, y, w] = pol2cart(th,r,w);

vx = vr.*cos(th) - vth.*sin(th);
vy = vr.*sin(th) + vth.*cos(th);

x_cart = [x y w vx vy vw];

end

