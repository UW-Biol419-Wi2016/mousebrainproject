function [ genebyregion ] = genebyregionmaker( diseasegenes, numgenes, braintable )
%GENEBYREGIONMAKER creates matrix for gene by region data
%   must have information for diseasegenes, numgenes

for i = 1:numgenes
    braingenestemp = strncmpi(diseasegenes{i, 1}, braintable{:, 1}, 10); 
    matches = max(braingenestemp);
    if matches == 1;
    %if there is a disease gene match, then:
   
            temploc = find(max(braingenestemp));
            %finds location of max of braingenestemp (where match occurred)
            
            genebyregion(i, :) = brainregionmat(temploc, :);
            %makes the ith row of genebyregion into the temploc row of
            %brainregionmat, which is arranged the same as braintable
            %(where the temploc index is from)
        end;
end;
