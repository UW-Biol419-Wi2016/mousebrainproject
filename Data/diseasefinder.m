function [ genestable ] = diseasefinder( bdt, str, strnum )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

genestable = table([], [], 'VariableNames', {'GeneName', 'DiseaseInfo'});
findpark = find(strncmpi(bdt{:, 2}, str, strnum)); %indices within bdt for parkinson disease genes
    %adding the parkinson genes to their table
        %calls gene name (column 1) at all findpark locations
        GeneName = bdt{findpark, 1};
        %calls disease info (column 5) at same indices
        DiseaseInfo = bdt{findpark, 2};
        %adds rows to the table, bdt, carrying the values from above

        genestable = [genestable; table(GeneName, DiseaseInfo)];

end

