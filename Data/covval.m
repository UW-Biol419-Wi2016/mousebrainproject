function [ C ] = covval( diseasemat, numcontrols, braintable, brainregionmat, bdt )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
diseasecell= {};
diseasemat = diseasemat;

[r, c] = size(bdt);

controlgenebyregion = [];
for i = 1:numcontrols

    numgenes = length(diseasemat');
    %number of genes in control shuffle

    randomperm = randperm(r);
    randomsel = randomperm(1:numgenes);
    %random indices of genes from bdt


    randomgenetable = table([], [], 'VariableNames', {'GeneName', 'DiseaseInfo'});
    randomtablepre = bdt(randomsel,:); 
    randomgenetable = [randomtablepre]; %adds shuffled
    %genes to table

    controlgenebyregion = genebyregionmaker(randomgenetable, braintable, brainregionmat);

    controlcov = cov(controlgenebyregion');
    %covariance matrix of ith control
  
    
    diseasecell{i} = controlcov; 
    %stores each successive controlcov in a cell for later analysis
    
end;

C = diseasecell;
end

