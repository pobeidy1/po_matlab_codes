
% gauss_2D for 2-D gaussians
% call as Y = gauss_2D(parameters, domain) with
% parameters = [Amp center_x center_y sigma_x sigma_y background]
 
function gauss = gauss_2DSymmetric(parameters, domain)
 
Amp = parameters(1);
center_x = parameters(2);
center_y = parameters(3);
sigma_x = parameters(4);
Background = parameters(5);


if length(domain)>2 %square domain
[D_x, D_y] = meshgrid(domain);
elseif length(domain)==2
[D_x, D_y] = meshgrid(domain{1},domain{2});
end

 
 
gauss = Amp*exp(-((((D_x - center_x).^2)./(2*sigma_x.^2)) + (((D_y - center_y).^2)./(2*sigma_x.^2))))+Background ;
