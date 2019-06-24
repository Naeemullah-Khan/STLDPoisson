clc
clear
close all
mex cg3d.cpp

im=imread('1.jpg');
mask=dlmread('1.txt');

tol=0.01;
alpha_vector=[5,10,15] ; %scale paramter of poisson equation
u = cg3d( double(im), int32(alpha_vector), int32(mask),double(tol) );

imagesc(image_tiling2(u))
colormap('gray')
figure
imagesc(image_tiling2(im))
colormap('gray')