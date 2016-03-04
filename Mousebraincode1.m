%Part 1 - setting up github

%% reading data
braintable = readtable('BrainRegionTotalDataset_Log2FoldChange.xlsx');
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
numgenesspec = height(diseasegenes);
%number of genes in table
genebyregion = [];
%initializes matrix for gene by region data


%'Cerebellum', 'Corpus Callosum', 'Motor Cortex', 'Olfactory Bulb', 'Optic Nerve', 'Prefrontal Cortex', 'Striatum', 'Thalamus', 'Hippocampus'});
%1:9 names of new double's columns
for i = 1:numgenesspec
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

%none of the parkinson's genes are considerably more expressed in one brain
%region than another

%visually looked at braintable to confirm this result

parkinsonsmat = genebyregionmaker(parkinsonsgenes, braintable, brainregionmat);

%% Making random shuffles of genes to make sets of genebyregion matrices (for control)
oneby1 = [];
oneby2 = [];
oneby3 = [];
oneby4 = [];
oneby5 = [];
oneby6 = [];
oneby7 = [];
oneby8 = [];
oneby9 = [];
twoby2 = [];
twoby3 = [];
twoby4 = [];
twoby5 = [];
twoby6 = [];
twoby7 = [];
twoby8 = [];
twoby9 =[];
threeby3 = [];
threeby4 = [];
threeby5 = [];
threeby6 = [];
threeby7 = [];
threeby8 = [];
threeby9 = [];
fourby4 = [];
fourby5 = [];
fourby6 = [];
fourby7 = [];
fourby8 = [];
fourby9 = [];
fiveby5 = [];
fiveby6 = [];
fiveby7 = [];
fiveby8 = [];
fiveby9 = [];
sixby6 = [];
sixby7 = [];
sixby8 = [];
sixby9 = [];
sevenby7 = [];
sevenby8 = [];
sevenby9 = [];
eightby8 = [];
eightby9 = [];
nineby9 = [];
%initalizing vectors to hold covariance values for each potential
%combination of brain regions(numbers here are the labels for each region
%as specified below)
%'Cerebellum', 'Corpus Callosum', 'Motor Cortex', 'Olfactory Bulb', 'Optic Nerve', 'Prefrontal Cortex', 'Striatum', 'Thalamus', 'Hippocampus'});
%1:9 names of new double's columns

numcontrols = 10000;
for i = 1:numcontrols

numgenes = 11;

randomperm = randperm(9669);
randomsel = randomperm(1:numgenes);


randomgenetable = table([], [], 'VariableNames', {'GeneName', 'DiseaseInfo'});
randomtablepre = bdt(randomsel,:); %shuffled genes
randomgenetable = [randomtablepre]; %adds shuffled
%genes to table

controlgenebyregion = genebyregionmaker(randomgenetable, braintable, brainregionmat);

controlcov = cov(controlgenebyregion);

oneby1(i) = controlcov(1, 1);
oneby2(i) = controlcov(1, 2);
oneby3(i) = controlcov(1, 3);
oneby4(i) = controlcov(1, 4);
oneby5(i) = controlcov(1, 5);
oneby6(i) = controlcov(1, 6);
oneby7(i) = controlcov(1, 7);
oneby8(i) = controlcov(1, 8);
oneby9(i) = controlcov(1, 9);

twoby2(i) = controlcov(2,2);
twoby3(i) = controlcov(2,3);
twoby4(i) = controlcov(2,4);
twoby5(i) = controlcov(2,5);
twoby6(i) = controlcov(2,6);
twoby7(i) = controlcov(2,7);
twoby8(i) = controlcov(2,8);
twoby9(i) = controlcov(2,9);

threeby3(i) = controlcov(3, 3);
threeby4(i) = controlcov(3, 4);
threeby5(i) = controlcov(3, 5);
threeby6(i) = controlcov(3, 6);
threeby7(i) = controlcov(3, 7);
threeby8(i) = controlcov(3, 8);
threeby9(i) = controlcov(3, 9);

fourby4(i) = controlcov(4,4);
fourby5(i) = controlcov(4,5);
fourby6(i) = controlcov(4,6);
fourby7(i) = controlcov(4,7);
fourby8(i) = controlcov(4,8);
fourby9(i) = controlcov(4,9);

fiveby5(i) = controlcov(5,5);
fiveby6(i) = controlcov(5,6);
fiveby7(i) = controlcov(5,7);
fiveby8(i) = controlcov(5,8);
fiveby9(i) = controlcov(5,9);

sixby6(i) = controlcov(6, 6);
sixby7(i) = controlcov(6, 7);
sixby8(i) = controlcov(6, 8);
sixby9(i) = controlcov(6, 9);

sevenby7(i) = controlcov(7, 7);
sevenby8(i) = controlcov(7, 8);
sevenby9(i) = controlcov(7, 9);

eightby8(i) = controlcov(8,8);
eightby9(i) = controlcov(8, 9);

nineby9(i) = controlcov(9, 9);

end;
    
    
%visually looked at braintable to confirm this result

%% Histograms of controls

hist45 = histogram(fourby5);

meanhist = mean(fourby5)
%% Measuring covariance of test set

covdisease = cov(genebyregion);


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



%% experimenting with tables and cells
covariancetable = {'oneby1', 
'oneby2',
'oneby3', 
'oneby4', 
'oneby5',
'oneby6',
'oneby7',
'oneby8', 
'oneby9',
'twoby2',
'twoby3',
'twoby4',
'twoby5',
'twoby6',
'twoby7',
'twoby8',
'twoby9',
'threeby3',
'threeby4',
'threeby5',
'threeby6',
'threeby7',
'threeby8',
'threeby9',
'fourby4',
'fourby5',
'fourby6',
'fourby7',
'fourby8',
'fourby9',
'fiveby5',
'fiveby6',
'fiveby7',
'fiveby8',
'fiveby9',
'sixby6',
'sixby7',
'sixby8',
'sixby9',
'sevenby7',
'sevenby8',
'sevenby9',
'eightby8',
'eightby9',
'nineby9'};

%for i = 1:45
    %searches the table of names
    %covariancetable{i}
    
   % eval(covariancetable(1))
    
    




