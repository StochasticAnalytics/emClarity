function [ RMSD ] = BH_returnImage2ImageRMSD(image1,image2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


RMSD = sqrt(mean(mean(mean((image1-image2).^2))));

end

