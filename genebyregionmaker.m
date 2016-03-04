function [ genebyregion ] = genebyregionmaker( diseasegenes, braintable, brainregionmat)
%GENEBYREGIONMAKER creates matrix for gene by region data
%   must have information for diseasegenes, numgenes
braingenestemp = 0
temploc = 0
for i = 1:height(diseasegenes)
    braingenestemp = strncmpi(diseasegenes{i, 1}, braintable{:, 1}, 10); 
    matches = max(braingenestemp);
    if matches == 1;
    %if there is a disease gene match, then:
   
            temploc = find(braingenestemp == 1);
            %finds location of max of braingenestemp (where match occurred)
            
            genebyregion(i, :) = brainregionmat(temploc, :);
            %makes the ith row of genebyregion into the temploc row of
            %brainregionmat, which is arranged the same as braintable
            %(where the temploc index is from)
        end;
end;