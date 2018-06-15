function [hiveFail,failedRhombuses,totalDefect] = rhombusCheck(Hijk,n)
%% Function for checking rhombus inequalities for AWHives
%outputs number of failed rhombuses and the failed rhombus indices as an f x 8 matrix.
%n is the vector space dimension, and Hijk is a sparse hive matrix.

totalDefect = 0;
    
%% Setup storage

hiveFail = 0;
failedRhombuses = zeros(n*(n-1)*3/2,8);

%% Loop over all possible indices and check the 3 rhombuses at each point

for k = 1:(n+1)

    for i = 1:(n+2-k)

        %% Right Inequality

            if i <= (n-k)

                temp = Hijk(k,i) + Hijk(k+1,i+1) - Hijk(k,i+1) - Hijk(k+1,i);

                if temp > 0

                    hiveFail = hiveFail + 1;
                    failedRhombuses(hiveFail,:) = [k i k+1 i+1 k i+1 k+1 i];
                    totalDefect = totalDefect + temp;

                end

            end


        %% Top Inequality

            if  i > 2

                temp = Hijk(k,i) + Hijk(k+1,i-2) - Hijk(k,i-1) - Hijk(k+1,i-1);

                if temp > 0
                    hiveFail = hiveFail + 1;
                    failedRhombuses(hiveFail,:) = [k i k+1 i-2 k i-1 k+1 i-1];
                    totalDefect = totalDefect + temp;

                end

            end

        %% Left Inequality

            if  k > 2

                temp = Hijk(k,i) + Hijk(k-2,i+1) - Hijk(k-1,i) - Hijk(k-1,i+1);

                if temp > 0
                    hiveFail = hiveFail + 1;
                    failedRhombuses(hiveFail,:) = [k i k-2 i+1 k-1 i k-1 i+1];
                    totalDefect = totalDefect + temp;

                end

            end

    end

end

%% Remove any unused rows in the failedRhombuses matrix

failedRhombuses(~any(failedRhombuses,2),:) = [];

end
