function A = RSMGenerator(n)
%% Generates real symmetric nxn matrices with entries bounded by +-1

A = 2*rand(n)-1;
A = 0.5*(A+A');

end