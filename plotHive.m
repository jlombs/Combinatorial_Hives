function plotHive(Hijk,m,tol)

%% Plot numerical Hive and Rhombus Failures if Applicable

[numFailures,rhombuses,totalDefect] = rhombusCheck(Hijk,m);

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

    if abs(round(v(q))-v(q)) < tol
        
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

    title(sprintf('Good Hive of Order %d',m))

else

    title(sprintf('Failed Hive of Order %d with \n %d Failures and %3.2e Total Defect',m,numFailures,totalDefect))

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
return
%% Plot Hive Surface

figure()

% Add in the origin
k = [k;0];
i = [i;0];
v = [v;0];

% Compute the triangulation
DT = delaunay([k+textCentering,i]);

% Remove exterior triangulations by an edge length cutoff
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

    title(sprintf('Good Hive Surface of Order %d',m))

else

    title(sprintf('Failed Hive Surface of Order %d with \n %d Failures and %3.2e Total Defect',m,numFailures,totalDefect))

end

%% Plot the Mean and Gaussian Curvatures for the Hive Surface

[GC, MC]= curvatures(k,i,v,DT);

figure()
patch('Faces',DT,'Vertices',[k,i,v],'FaceVertexCData',GC,'FaceColor','interp','EdgeColor','none')
colormap jet
colorbar
title('Estimated Gaussian Curvature');


figure()
patch('Faces',DT,'Vertices',[k,i,v],'FaceVertexCData',MC,'FaceColor','interp','EdgeColor','none')
colormap jet
colorbar
title('Estimated Mean Curvature');

end