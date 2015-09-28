function [mu,Sigma] = updateStatistics(mu,Sigma,theta,i,d,r)

% Updating of mu and Sigma
%mu = mu + (theta-mu) * (1+d*i)/(i+1+d*i);
%Sigma = i/(i+1+d*i)*(Sigma + (theta-mu)*(theta-mu)' * (1+d*i)/(i+1+d*i));
mu = mu + (theta-mu) * (1+d*i)/(i+1+d*i);
Sigma = i/(i+1+d*i)*Sigma + (1+d*i)/(i+1+d*i)*(theta-mu)*(theta-mu)';

% Regularisation
[~,p] = cholcov(Sigma,0);
if p ~= 0
    Sigma = Sigma + r*eye(size(dtheta));
end

end

