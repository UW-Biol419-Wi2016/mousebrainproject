function [ C ] = covvalcell( diseasemat, numcontrols, braintable, brainregionmat, bdt)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    covvalcell = {};

controlgenebyregion = [];
for i = 1:numcontrols

    numgenes = length(diseasemat');
    %number of genes in control shuffle

    randomperm = randperm(9669);
    randomsel = randomperm(1:numgenes);
    %random indices of genes from bdt


    randomgenetable = table([], [], 'VariableNames', {'GeneName', 'DiseaseInfo'});
    randomtablepre = bdt(randomsel,:); 
    randomgenetable = [randomtablepre]; %adds shuffled
    %genes to table

    controlgenebyregion = [];
    for i = 1:height(randomgenetable)
    braingenestemp = strncmpi(randomgenetable{i, 1}, braintable{:, 1}, 10); 
    matches = max(braingenestemp);
    if matches == 1;
    %if there is a disease gene match, then:
   
            temploc = find(braingenestemp == 1);
            %finds location of max of braingenestemp (where match occurred)
            
            controlgenebyregion(i, :) = brainregionmat(temploc(1), :);
            %makes the ith row of genebyregion into the temploc row of
            %brainregionmat, which is arranged the same as braintable
            %(where the temploc index is from)
        end;
end;


    controlcov = cov(controlgenebyregion');
    %covariance matrix of ith control
  
    
    covvalcell{i} = controlcov; 
    %stores each successive controlcov in a cell for later analysis
    
end;

C = covvalcell

end