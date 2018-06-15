function hiveEnsembleAveragesFctn()

% Test Script for exploring Hives, uses parallel sampling
% dim must be >2
% Computes the mean hive surfaces for various matrix ensembles, as well as
% average curvatures. See helper function below for setting the ensembles

clear
close all
clc

%% Vector space dimension
m = 15;

%% Total number of samples
sampleMax = 3*10;

%% Tollerance
tol = 10^-7;

%% Establish storage
hivesStorage = cell(1,sampleMax);
GCStorage = cell(1,sampleMax);
MCStorage = cell(1,sampleMax);

%% Flags must be unset%
flagData.waitbar = 0; %
flagData.verb = 0;    %
flagData.parallel = 0;%
%%%%%%%%%%%%%%%%%%%%%%%
    
%% Gather Samples

% Parallel loop over samples

 %% Setup Job
myCluster = parcluster; % Use default profile
j = createJob(myCluster,'Name','CollectHives','AttachedFiles','/home/jlombs/Research/Mixing/MixingCode/hives/taskStartup.m');

%% Setup Tasks
inputParams = cell(1,sampleMax);

for i = 1:sampleMax

    inputParams{i} = {m,flagData,tol};

end

t = createTask(j,@parallelWork,3,inputParams);

%% Submit job, wait for completion, retrive and organize data

submit(j);
wait(j)

data = fetchOutputs(j);

for i = 1:sampleMax
   
    hivesStorage{i} = data{i,1};
    GCStorage{i} = data{i,2};
    MCStorage{i} = data{i,3};
    
end

%% Delete data off of workers
delete(j)

%% Build mean hive and compute average GC and MC

meanHive = hivesStorage{1};
GCAvg = GCStorage{1};
MCAvg = MCStorage{1};
for i = 2:sampleMax
    
    meanHive = meanHive + hivesStorage{i};
    GCAvg = GCAvg + GCStorage{i};
    MCAvg = MCAvg + MCStorage{i};
    
end

Hijk = meanHive./sampleMax;
GCAvg = GCAvg./sampleMax;
MCAvg = MCAvg./sampleMax;

%% Plots--code largely reproduced from plotHive.m with tweaks to the labelings

% Plot numerical Hive and Rhombus Failures if Applicable

[numFailures,rhombuses,totalDefect] = rhombusCheck(Hijk,m,tol);

%Shift can be optimized for pretty pictures depending on scale of image and
%resolution of viewing screen.

textCentering = -.05;

%Get relevent coordinates and values
[k,i,v] = find(Hijk);

%Tweak values to lie on a nice triangle with the right orientation
for q = 1:numel(k)

    i(q) = i(q)-1;
    k(q) = k(q)-1 + i(q)*.5;

end

%Plot values in the triangle
plot(k+textCentering,i,'w*');
hold on
ln = findobj('type','line');
set(ln,'marker','.','markers',14,'markerfa','w')

%Set formatting properly
for q = 1:nnz(Hijk)

    if round(v(q))-v(q) < tol
        
        v(q) = round(v(q));
        textFormat = '%d';
        
    else
        
        textFormat = '%3.2f';
        
    end
    text(k(q)+textCentering,i(q),sprintf(textFormat,v(q)))

end

%Plot the 0 normalized value
text(textCentering,0,'0')

%Label
xlabel('K')
ylabel('I')

if numFailures == 0

    title(sprintf('%d-Average Hive of Order %d is a Good Hive',sampleMax,m))

else

    title(sprintf('%d Average Hive Hive of Order %d has \n %d Failures and %3.2e Total Defect',sampleMax,m,numFailures,totalDefect))

end

axis([-1 (m+1) -1 (m+1)])

%Highlight any failed rhombus inequalities by looping over them
for q = 1:numFailures

    %Create the vertex list from rhombuses
    vert = zeros(4,2);
    count = 1;
    for p = 1:2:7

        vert(count,:) = [rhombuses(q,p) rhombuses(q,p+1)];
        count = count + 1;
    end
    
    %Again shift for good coordinates
    for p = 1:4

        vert(p,2) = vert(p,2)-1;
        vert(p,1) = vert(p,1)-1 + vert(p,2)*.5;

    end

    %Set adjacency relationship between all rhombuses and their 4 vertices
    faces = [1 3 2 4 1];
    
    %Plot the semiopaque rhombus
    pb = patch('Faces',faces,'Vertices',vert,'FaceColor','r','EdgeColor','k');
    hold on
    alpha(pb,.1);

end

%% Plot Hive Surface

figure()

% Add origin point
k = [k;0];
i = [i;0];
v = [v;0];

% Compute triangulation
DT = delaunay([k+textCentering,i]);

% Remove exterior triangulations
DTRem = false(1,size(DT,1));
for q = 1:size(DT,1)
   
    if (k(DT(q,1)) - k(DT(q,2)))^2 +(i(DT(q,1)) - i(DT(q,2)))^2 > 2 || ...
        (k(DT(q,1)) - k(DT(q,3)))^2 +(i(DT(q,1)) - i(DT(q,3)))^2 > 2 || ...
        (k(DT(q,3)) - k(DT(q,2)))^2 +(i(DT(q,3)) - i(DT(q,2)))^2 > 2
        
        DTRem(q) = true;
        
    end
    
end

DT(DTRem,:) = [];

tr = triangulation(DT,k,i,v);
trisurf(tr)


%Highlight any failed rhombus inequalities by looping over them
for q = 1:numFailures

    %Create the vertex list from rhombuses
    vert = zeros(4,3);
    count = 1;
    for p = 1:2:7

        vert(count,1:2) = [rhombuses(q,p) rhombuses(q,p+1)];
        vert(count,3) = Hijk(vert(count,1),vert(count,2));
        count = count + 1;
    end
    
    %Again shift for good coordinates
    for p = 1:4

        vert(p,2) = vert(p,2)-1;
        vert(p,1) = vert(p,1)-1 + vert(p,2)*.5;

    end

    %Set adjacency relationship between all rhombuses and their 4 vertices
    faces = [1 3 2 4 1];
    
    %Plot the semiopaque rhombus
    pb = patch('Faces',faces,'Vertices',vert,'FaceColor','r','EdgeColor','k');
    hold on
    alpha(pb,.5);

end


if numFailures == 0

    title(sprintf('Average of %d %d-Hive Surfaces',sampleMax,m))

else

    title(sprintf('Failed %d Average Hive Surface of Order %d with \n %d Failures and %3.2e Total Defect',sampleMax, m,numFailures,totalDefect))

end

%% Plot the Mean and Gaussian Curvatures for the Average Hive Surface

%{
[GC, MC]= curvatures(k,i,v,DT);

figure()
patch('Faces',DT,'Vertices',[k,i,v],'FaceVertexCData',GC,'FaceColor','interp','EdgeColor','none')
colormap jet
colorbar
title(sprintf('Estimated Gaussian Curvature for the %d Average Hive',sampleMax));


figure()
patch('Faces',DT,'Vertices',[k,i,v],'FaceVertexCData',MC,'FaceColor','interp','EdgeColor','none')
colormap jet
colorbar
title(sprintf('Estimated Mean Curvature for the %d Average Hive',sampleMax));

%}

%% Plot the Average Mean and Gaussian Curvatures

figure()
patch('Faces',DT,'Vertices',[k,i,v],'FaceVertexCData',GCAvg,'FaceColor','interp','EdgeColor','none')
colormap jet
colorbar
title(sprintf('Estimated Average Gaussian Curvature for %d %d-Hives',sampleMax,m));


figure()
patch('Faces',DT,'Vertices',[k,i,v],'FaceVertexCData',MCAvg,'FaceColor','interp','EdgeColor','none')
colormap jet
colorbar
title(sprintf('Estimated Average Mean Curvature for %d %d-Hives',sampleMax,m));

%% Graph the leading boundary edge of the mean hive

bndVals = v(i==0);
bndVals = [bndVals(end);bndVals];
bndVals(end) = [];
figure()
plot(1:numel(bndVals),bndVals)
title(sprintf('Projected Boundary Curve for the %d-Average Hive Surface of Order %d',sampleMax,m));

end

function [hive,GC,MC] = parallelWork(m,flagData,tol)

    while 1
%% Setup matrix ensembles
        
        %{
        % Diagonally Dominant SPD matrices. sames
        L = RSMGenerator(m);
        L = L + m*eye(m);
        M = L;
        %}

        %{
        % Normally distributed SPD, sames
        L = randn(m);
        L = L'*L;
        M = L;
        %}
        
        %{
        % Diagonal matrices with weakly decreasing positive integer values 
        M = diag(sort(randi(50,[1,m]),'descend'));
        L = diag(sort(randi(50,[1,m]),'descend'));
        %}

        %{
        % Diagonal matrices with weakly decreasing real values less than abs(x)
        x = 500;
        M = diag(sort(x*(2*rand([m,1])-1),'descend'));
        L = diag(sort(x*(2*rand([m,1])-1),'descend'));
        %}
        
        %%{
        %GOE matrices, sames
        L = GOEGenerator(m);
        M = L;
        %}


        %% Perform optimization with up to 4 re-runs for the same boundary data
        hiveLogical = 0;
        trials = 0;
        
        while hiveLogical == 0 && trials < 4

            Hijk = AWHiveParallel(L,M,flagData);
            
            if ~isempty(Hijk)
                
                hiveLogical = rhombusCheckF(Hijk,m,tol);
                
            elseif trials == 4

                hiveLogical = false;
                
            end
            
            trials = trials + 1;

        end
        
        %% Compute triangulation and curvatures and store the results if it's a good hive
        if hiveLogical
            
            hive = Hijk;
            
            textCentering = -.05;
            [k,j,v] = find(Hijk);
            for q = 1:numel(k)

                j(q) = j(q)-1;
                k(q) = k(q)-1 + j(q)*.5;

            end
            k = [k;0];
            j = [j;0];
            v = [v;0];
            DT = delaunay([k+textCentering,j]);
            
            DTRem = false(1,size(DT,1));
            for q = 1:size(DT,1)

                if (k(DT(q,1)) - k(DT(q,2)))^2 +(j(DT(q,1)) - j(DT(q,2)))^2 > 2 || ...
                    (k(DT(q,1)) - k(DT(q,3)))^2 +(j(DT(q,1)) - j(DT(q,3)))^2 > 2 || ...
                    (k(DT(q,3)) - k(DT(q,2)))^2 +(j(DT(q,3)) - j(DT(q,2)))^2 > 2

                    DTRem(q) = true;

                end

            end

            DT(DTRem,:) = [];

            [GC, MC]= curvatures(k,j,v,DT);
            
            break
            
        end
        
    end

end