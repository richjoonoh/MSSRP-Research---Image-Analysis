function oindex = oindexCalc(filepathI,mask,filepathR)
% PURPOSE:
%  image analysis program that
%  1) calculates orientation index of fibers in the image from filepathI
%  2) generates corresponding scattering curve
% 
% FUNCTION OUTPUTS:
%  oindex:    orientation index
%  figure:    scattering curve plot, which when integrated and normalized
%             yields oindex
% 
% FUNCTION INPUTS:
%  filepathI: filepath of .png post-processed image for analysis; image is
%             post-processed by imageJ plugin TWOMBLI; image should only
%             display a region of interest (manually selected);
%             filepathI file type = string
%  mask:      binary matrix whose size is same as .png image from
%             filepathI; 1's in the matrix represent portions of the .png
%             image that are analyzed while 0's represent excluded portions
%  filepathR: filepath of .zip file containing coordinates of polygonal
%             ROIs that outline chondrocytes and other artifacts; the ROI
%             .zip file obtained from imageJ; used for ROI exclusion for
%             image analysis; filepathR file type = string
% 
% AUTHOR COMMENTS:
%  this image analysis program simulates light scattering as if the image
%  is a real histological section through which light passes; calculation
%  of orientation index is exact given the image since the program scatters
%  simulated light pixel by pixel

% ----------------------------------------------------------------------
% find points surrounding fibers - in effect outlining fibers
% goal is to locate the points where light is intercepted and scatters
% ----------------------------------------------------------------------
I = imbinarize(imread(filepathI));
I(~mask) = 1;
scatterpoints = I;
[row,col] = find(scatterpoints == 0);
scatterpoints = zeros(size(I));
for i = 1:length(row)
    if row(i) > 2 && col(i) > 2 && row(i) < size(I,1) && col(i) < size(I,2)
        scatterpoints(row(i)-1:row(i)+1,col(i)-1:col(i)+1) = ...
            scatterpoints(row(i)-1:row(i)+1,col(i)-1:col(i)+1)+1;
    else
        r1 = 1; r2 = size(I,1); c1 = 1; c2 = size(I,2);
        if row(i)-1 > r1
            r1 = row(i)-1;
        end
        if row(i)+1 < r2
            r2 = row(i)+1;
        end
        if col(i)-1 > c1
            c1 = col(i)-1;
        end
        if col(i)+1 < c2
            c2 = col(i)+1;
        end
        scatterpoints(r1:r2,c1:c2) = scatterpoints(r1:r2,c1:c2)+1;
    end
end
scatterpoints = (scatterpoints < 9) - (scatterpoints == 0);
scatterpoints = scatterpoints.*I;
% remove points that are from chondrocytes and artifacts, and not fibers
rois = ReadImageJROI(filepathR);
minusscatterpoints = ones(size(I));
for iroi = 1:length(rois)
    bw = poly2mask([rois{1,iroi}.vnRectBounds(2) ...
                    rois{1,iroi}.vnRectBounds(2) ...
                    rois{1,iroi}.vnRectBounds(4) ...
                    rois{1,iroi}.vnRectBounds(4)],...
                    [rois{1,iroi}.vnRectBounds(1) ...
                    rois{1,iroi}.vnRectBounds(3) ...
                    rois{1,iroi}.vnRectBounds(3) ...
                    rois{1,iroi}.vnRectBounds(1)],size(I,1),size(I,2));
    I(logical(bw)) = 1;
    minusscatterpoints(logical(bw)) = 0;
end
scatterpoints = scatterpoints.*minusscatterpoints;
% remove points along the image boundaries as light scattering at these
% points will contain boundary effects that are not physical
[row,col] = find(scatterpoints == 1);
k = find(row == 1); k = [k; find(row == size(scatterpoints,1))];
col(k) = []; row(k) = [];
k = find(col == 1); k = [k; find(col == size(scatterpoints,2))];
col(k) = []; row(k) = [];

% ----------------------------------------------------------------------
% NOTE: light scattering mechanism shown below is discrete and linear,
% and yields an orientation index (calculated as 2*<cos^2(alpha)> - 1)
% that still spans from 0 to 1 densely and uniformly
% ----------------------------------------------------------------------
vectors = [1 1/sqrt(2) 0 -1/sqrt(2) -1 -1/sqrt(2) 0 1/sqrt(2);... %x-dir
           0 -1/sqrt(2) -1 -1/sqrt(2) 0 1/sqrt(2) 1 1/sqrt(2)]; %y-dir
%corresponding position relative to fiber point going clockwise from left
pos = [-1 -1 0 1 1 1 0 -1;0 1 1 1 0 -1 -1 -1];
scatteredvectors = zeros(2,length(row));
for i = 1:length(row)
    vectorsum = vectors;
    for count = 1:size(pos,2)
        if I(row(i)-pos(2,count),col(i)+pos(1,count)) == 1
            vectorsum(:,count) = 0;
        end
    end
        scatteredvectors(:,i) = sum(vectorsum,2);
end
% find angles of scattered light point source
scatteredangles = scatteredvectors(1,:);
for i = 1:size(scatteredvectors,2)
    if scatteredvectors(1,i) >= 0
        scatteredangles(i) = rad2deg(atan(scatteredvectors(2,i)/...
                                            scatteredvectors(1,i)));
    elseif scatteredvectors(2,i) >= 0
        scatteredangles(i) = rad2deg(pi+atan(scatteredvectors(2,i)/...
                                               scatteredvectors(1,i)));
    else
        scatteredangles(i) = rad2deg(-pi+atan(scatteredvectors(2,i)/...
                                                scatteredvectors(1,i)));
    end
end
% NOTE: oindex = 2*<cos^2(alpha)> - 1, where alpha is angle between a fiber
% and the mean fiber axis (denoted as meanaxisangle)
meanaxisangle = mode(scatteredangles(scatteredangles >= 0));
k_aboutmeanaxis = find(scatteredangles <= meanaxisangle+90 & ...
                        scatteredangles > meanaxisangle-90);
oindex = 2*sum(cos(deg2rad(scatteredangles(k_aboutmeanaxis)-...
                            meanaxisangle)).^2)/length(k_aboutmeanaxis)-1;
% scattering curve generation
figure
[N,edges] = histcounts(scatteredangles(~isnan(scatteredangles)));
edges = edges(2:end) - (edges(2)-edges(1))/2;
plot(edges, N)
title(['orientation index=' num2str(oindex) ', mean axis angle=' ...
        num2str(meanaxisangle)])
xlabel('angle from +x (deg)')
ylabel('count')
end