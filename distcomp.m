function [ genedist ] = distcomp( genecov, diseasemean, diseasestd )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
genedist = zeros(length(genecov));
for i= 1:length(genecov)
    for j = 1:i
    genedist(i, j) = (genecov(i, j) - diseasemean(i, j))/diseasestd(i, j);
    %for each value of parkinsons covariance, subtract the mean from that
    %value and divide by the std for that value to get the mahalanobis
    %distance for each brain region 
    end;
end

genedist = genedist;
end

