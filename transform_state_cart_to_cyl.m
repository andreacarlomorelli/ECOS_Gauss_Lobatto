function [x_cyl] = transform_state_cart_to_cyl(x_cart)

x = x_cart(:,1); y = x_cart(:,2); w = x_cart(:,3); 
vx = x_cart(:,4); vy = x_cart(:,5); vw = x_cart(:,6);

[th, r, w] = cart2pol(x,y,w);

vr = (x.*vx + y.*vy)./sqrt(x.^2 + y.^2);

thdot = (x.*vy - y.*vx)./(x.^2 + y.^2);

vth = r.*thdot;

x_cyl = [r w vr vth vw];

end

