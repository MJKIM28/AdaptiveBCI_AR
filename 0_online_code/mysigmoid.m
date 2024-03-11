function G = mysigmoid(U,V)
% Sigmoid kernel function with slope gamma and intercept c
global gamma
c = 0;
G = tanh(gamma*U*V' + c);
end