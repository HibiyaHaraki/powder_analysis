clear; close all; clc;
%% Import figure
I = imread("data1004450.bmp");
sz = size(I);
bk = 50;
converter = 0.1754;
%% Create figure
imshow(I(282:468,1:1023));
hold on;
for x = 1:bk:sz(2)
    line([x x],[0 sz(1)],'color','k');
end
for y = 1:bk:sz(1)
    line([0 sz(2)],[y y],'color','k');
end
hold off
axis on
ax = gca;
xticks(ax,1:bk:sz(2));
xticklabels(ax,{((1:bk:sz(2))-1).*converter});
yticks(ax,1:bk:sz(1));
yticklabels(ax,{((1:bk:sz(1))-1).*converter});