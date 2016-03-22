function [ diseasemean, diseasestd ] = cellmeanfinder2( diseasemat, diseasecell, n )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

diseasemean = [];
diseasestd = [];
temp = zeros(length(cov(diseasemat)));

    for i = 1:length(cov(diseasemat))
        %runs length of each control rows
        for j = 1:i
            for k = length(diseasecell)
            %runs for each row, proper num of col
            temp(i, j) = temp(i, j) + diseasecell{1, k}(i, j);         
            end;
        diseasemean(i, j) = temp(i, j)/n;
        diseasestd(i, j) = sqrt((temp(i, j) - diseasemean(i, j)).^2/(n-1));
        end;
end;

diseasemean = diseasemean;
diseasestd = diseasestd;
end

