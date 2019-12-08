clc,clear;

clear all; close all; clc;

I=imread('e.jpg');
img = rgb2gray(I);
 ed=edge(img,'canny');

for i=0:1:179
r=radon(ed);     %检测直线什么的，可以投影到0-179度上
end

r2=zeros(755,180);

for i=70:110
    r2(:,i)=r(:,i);
end
I2=iradon(r2,0:179);

figure,
subplot(121),imshow(I);
subplot(122),imagesc(r);
figure,
subplot(121),imagesc(r2);
subplot(122),imshow(I2)

% figure;
% plot(r)

