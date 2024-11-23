% -----------------------------------------------------------------------------
% Copyright (C) 2024 Marco Cattaneo
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the Free
% Software Foundation, either version 3 of the License, or (at your option)
% any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of  MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
% more details.
%
% You should have received a copy of the GNU General Public License along with
% this program.  If not, see <http://www.gnu.org/licenses/>.
% -----------------------------------------------------------------------------


close all; clear all; clc

PixelSize = 0.16e-6;

NumberFrames = 120;

for i = 1:NumberFrames

    n = num2str(i,'%03.f');

    % Load the image 
    name = ['Frame_',n,'.tiff'];
    image = imread(name);
    
    % Convert the image to grayscale
    grayImage = im2gray(image);
    
    % Filter the image 
    filteredImage = medfilt2(grayImage,[15 15]);
    
    % Binarise the image 
    binaryImage = imbinarize(filteredImage(40:190,140:280), 'adaptive','ForegroundPolarity','dark','Sensitivity',0.5); % Crop and sensitivity depends on the specific video under study
    
    % Remove small objects in the image
    binaryImage = imcomplement(binaryImage);
    binaryImage = bwareaopen(binaryImage, 30);

    % Fill the image
    filledImage = imfill(binaryImage, 'holes');
    for n = 1 : size(filledImage,1)
        BeginLine = find(filledImage(n,:) == 1, 1, 'first');
        EndLine = find(filledImage(n,:) == 1, 1, 'last');
        if ~isempty(BeginLine)
            filledImage(n,BeginLine:EndLine) = 1;
        end
    end
      
    % Find the bubble contour
    contourImage = edge(filledImage, 'Canny', [0.2, 0.3]);

    % Calculate the equivalent radius of the bubble
    V = 0;
    for n = 1 : size(contourImage,1)
        BeginLine = find(contourImage(n,:) == 1, 1, 'first');
        EndLine = find(contourImage(n,:) == 1, 1, 'last');
        if ~isempty(BeginLine)
            V = V + pi/4*(EndLine-BeginLine+1)^2;
        end
    end
    Req(i) = (3/(4*pi)*V)^(1/3);
    
    % Extract the centroid position of the bubble
    stats = regionprops(filledImage);
    if ~isempty(stats)
        centroid(i,1:2) = stats.Centroid;
    else
        centroid(i,1:2) = [0,0];
    end  
end

% Plot time-radius evolution
t = 1e-7:1e-7:1e-7*NumberFrames;
Req = Req*PixelSize;
f = figure(1);
f.Position = [0 0 1120 1120/4]; 
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex');
grid on
box on
hold on
plot(t, Req)
plot(t, Req,'o','MarkerSize',7);
xlabel('$t$ [s]','interpreter','latex')
ylabel('$R$ [m]','interpreter','latex')
ay=gca; ay.YAxis.Exponent = -6;
ax=gca; ax.XAxis.Exponent = -6;

% Plot time-vertical position evolution
centroid = centroid*PixelSize;
f = figure(2);
f.Position = [0 0 1120/3 1120/2]; 
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex');
grid on
box on
hold on
plot(t, abs(centroid(:,2)-(centroid(4,2))))
plot(t, abs(centroid(:,2)-(centroid(4,2))),'o','MarkerSize',7);
xlabel('$t$ [s]','interpreter','latex')
ylabel('$y$ [m]','interpreter','latex')
ay=gca; ay.YAxis.Exponent = -6;
ax=gca; ax.XAxis.Exponent = -6;

