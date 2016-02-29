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
    %compares strings of ith gene of braintable and all gene names of
    %disease table to find a match
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



%% Pulling out a diseases genes: Parkinsons
parkinsonsgenes = table([], [], 'VariableNames', {'GeneName', 'DiseaseInfo'});

findpark = find(strncmpi(bdt{:, 2}, 'Parkinson', 9));
%indices within bdt for parkinson disease genes

    %adding the parkinson genes to their table

        GeneName = bdt{findpark, 1};
        %calls gene name (column 1) at all findpark locations
        
        DiseaseInfo = bdt{findpark, 2};
        %calls disease info (column 5) at same indices 
        
        parkinsonsgenes = [parkinsonsgenes; table(GeneName, DiseaseInfo)];
        %adds rows to the table, bdt, carrying the values from above
        
%% Making brainregion matrix

%brain table sorting (turning into a matrix and tagging genes that are
%correlated with the chosen disease. First converting to double, then
%choosing subset of genes

brainregiontable = braintable{:, 25:33};
%isolates brain region expression info into a separate table

brainregionmat = [];
%double that stores expression data in logical

for j = 1:8715
    %searches all rows of brain region table
    
    for i = 1:9
        %searches the brain region row
    
        if strncmpi(brainregiontable{j, i}, '+', 3) == 1
            %if the value of brtable is 1 then:
        
            brainregionmat(j, i) = 1;
            %brain region matrix value set as 1
        else
         brainregionmat(j, i) = 0;
            %if value is empty, there will be a zero
        end;
    end;
end;

%% Making the final matrix, genebyregion
diseasegenes = parkinsonsgenes;

numgenes = 11;

genebyregion = [];
%initializes matrix for gene by region data

%'Cerebellum', 'Corpus Callosum', 'Motor Cortex', 'Olfactory Bulb', 'Optic Nerve', 'Prefrontal Cortex', 'Striatum', 'Thalamus', 'Hippocampus'});
%1:9 names of new double's columns


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

%none of the parkinson's genes are considerably more expressed in one brain
%region than another

%visually looked at braintable to confirm this result
%% Sorting the brain-disease table (finding disease-related values only)
%broken code that could sort bdt 

bdt2 = bdt;
%initializes table that only will contain genes with disease info

diseasenames = unique(bdt2.DiseaseInfo);
%makes a list of disease names found in bdt2

dnames = cell2table(diseasenames);

for i = 1:9669
    %searches the whole bdt table
    %if strncmpi(bdt2{i, 2}, '', 2);
    %if isempty(bdt2) == 1
            %if bdt disease info is empty, then
            bdt2(i,:) = [];
            %row i is deleted, leaving only disease info containing rows
    end;
end;





