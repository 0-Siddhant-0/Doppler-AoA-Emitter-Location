function ellipse(C, k, style, offset)
%
%  USAGE: ellipse(C, k, style, offset)
%
%  Draws error ellipse for a specific covariance matrix C
%
%  Outputs: Plots the ellipse
%
%  Inputs: C = 2x2 covariance of the elements whose ellipse is to be plotted
%          k = value to set size of ellipse
%               (k = -2 ln(1-Pe) sets ellipse so that Pe is
%               probability for a location estimate to be inside the
%               ellipse - see Torrieri's paper)
%          style = text describing line style for plot
%          offset = [xe ye] location to center ellipse there (default: [0 0])

% Handle optional offset parameter
if nargin < 4
    offset = [0 0];
end

J_hat = inv(C);

% Compute eigenvectors and eigenvalues of J_hat
[V, D] = eig(J_hat);

lambda1 = min(D(1,1), D(2,2));
lambda2 = max(D(1,1), D(2,2));

l1 = sqrt(k/lambda1); 
l2 = sqrt(k/lambda2);

u1 = sqrt(1/D(1,1))*V(:,1); 
u2 = sqrt(1/D(2,2))*V(:,2);

% Compute ellipse that has semi-axes along coordinate axes
x = linspace(-l1, l1, 200);
y1 = l2*sqrt(1 - (x/l1).^2);
y2 = -y1;

% Transform ellipse into correct orientation of semi-axes
if norm(u1) > norm(u2)
  theta = atan2(u1(2), u1(1));
else
  theta = atan2(u2(2), u2(1));
end

R = [cos(theta) -sin(theta); sin(theta) cos(theta)];

old1 = [x; y1]; 
old2 = [x; y2];

new1 = R*old1;
new2 = R*old2;

% Apply offset and plot ellipse
h = plot(new1(1,:)+offset(1), new1(2,:)+offset(2), style, new2(1,:)+offset(1), new2(2,:)+offset(2), style);
set(h, 'linewidth', 2);

end