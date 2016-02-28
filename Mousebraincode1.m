%Part 1 - setting up github


%% reading data

braintable = readtable('Copy of BrainRegionTotalDataset_Log2FoldChange.xlsx');
diseasetable = readtable('MGIdisease.txt', 'delimiter', 'tab', 'readvariablenames', 0);

% testing if names are the same


%% string comparisons
bdt = table([], [], 'VariableNames', {'GeneName', 'DiseaseInfo'});
%initializes and empty table with two columns, one for gene names, and one
%for disease info

for i = 1:8715
    %for loop to compare each element of braintable (gene names column)
    %with all gene names in disease table to see if there is a match
    
    braingenestemp = strncmpi(braintable{i, 1}, diseasetable{:, 1}, 10);
    matches = max(braingenestemp);
    %finds if there is a match (if there is at least one value of 1, and
    %there will not be more matches per vector)
    
    if matches == 1
        %if there is a match then:
        temploc = find(braingenestemp == 1);
        %index location where the match occurred 
        
        GeneName = diseasetable{temploc, 1};
        %calls gene name (column 1) at that location
        
        DiseaseInfo = diseasetable{temploc, 5};
        %calls disease info (column 5) at that location 
        
        bdt = [bdt; table(GeneName, DiseaseInfo)];
        %adds a row to the table, bdt, carrying the values from above
    end;
end;

%% Sorting the brain-disease table (finding disease-related values only)

bdt2 = bdt
%initializes table that only will contain genes with disease info

for i = 1:9669
    %searches the whole bdt table
    if strncmpi(bdt{i, 2}, '', 2);
        %if bdt disease info is empty, then
        bdt([i],:) = [];
        %row i is deleted, leaving only disease info containing rows
    end;
end;







