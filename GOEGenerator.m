function G = GOEGenerator(n)
%% Generates a random nxn matrix from the GOE.

G = randn(n)/sqrt(n);

G = (G + G')/sqrt(2);

end