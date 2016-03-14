%Part 1 - setting up github

%% reading data
braintable = readtable('/Users/esteligarcia/GitHub/Data/BrainRegionTotalDataset_Log2FoldChange.xlsx');
diseasetable = readtable('/Users/esteligarcia/GitHub/Data/MGIdisease.txt',...
    'delimiter', 'tab','readvariablenames', 0);
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

%% Pulling out a diseases genes: Parkinsons & Huntingtons
parkinsonsgenes = table([], [], 'VariableNames', {'GeneName', 'DiseaseInfo'});
findpark = find(strncmpi(bdt{:, 2}, 'Parkinson', 9)); %indices within bdt for parkinson disease genes
    %adding the parkinson genes to their table
        %calls gene name (column 1) at all findpark locations
        GeneName = bdt{findpark, 1};
        %calls disease info (column 5) at same indices
        DiseaseInfo = bdt{findpark, 2};
        %adds rows to the table, bdt, carrying the values from above
        parkinsonsgenes = [parkinsonsgenes; table(GeneName, DiseaseInfo)];

huntingtonsgenes = table([], [], 'VariableNames', {'GeneName', 'DiseaseInfo'});
    findhunt = find(strncmpi(bdt{:, 2}, 'Huntington',10));
    genenamehunt = bdt{findhunt, 1};
    diseaseinfohunt = bdt{findhunt, 2};
    huntingtonsgenes = [parkinsonsgenes; table(genenamehunt, diseaseinfohunt)];
%% Making brainregion matrix
%brain table sorting (turning into a matrix and tagging genes that are
%correlated with the chosen disease. First converting to double, then
%choosing subset of genes
brainregiontable = braintable{:, 24:33};
%isolates brain region expression info into a separate table

brainregionmat = [];
%double that stores expression data in logical

for j = 1:8715
    %searches all rows of brain region table
    
    for i = 1:10
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

%'brain stem', 'Cerebellum', 'Corpus Callosum', 'Motor Cortex', 'Olfactory Bulb', 'Optic Nerve', 'Prefrontal Cortex', 'Striatum', 'Thalamus', 'Hippocampus'});
%1:9 names of new double's columns
parkinsonsmat = genebyregionmaker(diseasegenes, braintable, brainregionmat);
%insert other genes here

%% Making random shuffles of genes to make sets of genebyregion matrices (for control)
 n = 10 ;
    covvalcell = cell(n, 1) ;
    %makes a cell matrix to store all matrices from the loop

%initalizing vectors to hold covariance values for each potential
%combination of brain regions(numbers here are the labels for each region
%as specified below)
%'Cerebellum', 'Corpus Callosum', 'Motor Cortex', 'Olfactory Bulb', 'Optic Nerve', 'Prefrontal Cortex', 'Striatum', 'Thalamus', 'Hippocampus'});
%1:9 names of new double's columns
controlgenebyregion = [];
numcontrols = 100;
for i = 1:numcontrols

    numgenes = 11;
    %number of genes in control shuffle

    randomperm = randperm(9669);
    randomsel = randomperm(1:numgenes);
    %random indices of genes from bdt


    randomgenetable = table([], [], 'VariableNames', {'GeneName', 'DiseaseInfo'});
    randomtablepre = bdt(randomsel,:); 
    randomgenetable = [randomtablepre]; %adds shuffled
    %genes to table

    controlgenebyregion = genebyregionmaker(randomgenetable, braintable, brainregionmat);

    controlcov = cov(controlgenebyregion);
    %covariance matrix of ith control
    
    genecontrolcov = cov(controlgenebyregion');
    %covariance genebygene
    
    covvalcell{i} = controlcov; 
    %stores each successive controlcov in a cell for later analysis
    
    genecovvalcell{i} = genecontrolcov;
end;
    
%visually looked at braintable to confirm this result

%% Same as below (obtaining histograms and others) but from gene by gene cov matrices
num = {};

for i=1:length(genecontrolcov);
   eval(sprintf('temp%d', i));
end;
%loop that creates variable names "temp" for each row of gene by gene table
%and stores them in a cell for later recall


genecovtemp = zeros(length(genecontrolcov));



%%
for j = 1:length(genecovvalcell)
    %runs control
    
    for i = 1:length(genecontrolcov)
        %runs for each gene of gene by gene
        eval(sprintf('temp%d', i));
        
        temp1(j, 1:i) =  genecovvalcell{j, 1}(i, 1:i);
  
      
        %takes the values of covvalmatrix and puts them in temp variables of name
        %respective to row number in gene by gene region
    end;
    
end;




%% Histograms of controls
%lol = covvalcell{1,1}(1,1); %returns 0 = OK
%indexing pattern for retrieving values within arrays within a cell

%matrices holding each successive rows combinations data for histograms
covtemp = zeros(10);

temp1 = [];
temp2 = [];
temp3 = []; 
temp4 = [];
temp5 = [];
temp6 = [];
temp7 = [];
temp8 = [];
temp9 = [];
temp10 = [];

%how to read temp matrices. second brain region is number on matrix name.
%first brain region is column of matrix.

for j = 1:length(covvalcell)
    for i= 1:10
        covtemp(i, 1:i) = covvalcell{j, 1}(i, 1:i);
        temp1(j, 1:1) = covtemp(1, 1:1);
        temp2(j, 1:2) = covtemp(2, 1:2);
        temp3(j, 1:3) = covtemp(3, 1:3);
        temp4(j, 1:4) = covtemp(4, 1:4);
        temp5(j, 1:5) = covtemp(5, 1:5);
        temp6(j, 1:6) = covtemp(6, 1:6);
        temp7(j, 1:7) = covtemp(7, 1:7);
        temp8(j, 1:8) = covtemp(8, 1:8);
        temp9(j, 1:9) = covtemp(9, 1:9);
        temp10(j, 1:10) = covtemp(10, 1:10);
    end;
    
end;


%%


num = {temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10};

for i = 1:10
[trow(i), tcol(i)] = size(num{i});
end;
%vector of column numbers for each temp matrix

%% Finding Control Mean and Std

controlstd = [];

controlmean = [];


for j = 1:10
    %running for each temp matrix
    for i = 1:tcol(j)
        %running for each column of the given matrix
        formatSpec = 'Region %d by region %d';
        A1 = i;
        A2 = j;
        str = sprintf(formatSpec,A1,A2);
        %names each histogram the correct value
        %figure;
        %histogram(num{j}(:, i), 10);
        %title(str);
        %graphs histogram of all rows (controls) and i columns (one
        %histogram for each brain region comparison
        controlmean(j, i) = mean(num{j}(:, i));
        %adds into matrix the means of each combo (lower left triangle
        %holds means in same order as controlcov (regionwise))
        controlstd(j, i) = std(num{j}(:, i));
    end;
end; 


%% Calculating std from mean for parkinsons vs controls
%parkinsonsmat
%controlstd
%controlmean
parkinsonsdist = [];
parkinsonscov = cov(parkinsonsmat);
for i= 1:10
    parkinsonsdist(i, 1:i) = (controlmean(i, 1:i) - parkinsonscov(i,1:i))/controlstd(i, 1:i);
    %for each value of parkinsons covariance, subtract the mean from that
    %value and divide by the std for that value to get the mahalanobis
    %distance for each brain region 
end;

%% Calculating std from mean for parkinsons vs controls
%parkinsonsmat
%controlstd
%controlmean
parkinsonsdist = [];
parkinsonscov = cov(parkinsonsmat);
for i= 1:10
    parkinsonsdist(i, 1:i) = (controlmean(i, 1:i) - parkinsonscov(i,1:i))/controlstd(i, 1:i);
    %for each value of parkinsons covariance, subtract the mean from that
    %value and divide by the std for that value to get the mahalanobis
    %distance for each brain region 
end;

imagesc(parkinsonsdist)
colorbar
title('Brain Region x Brain Region Interactions Relating to Parkinsons Gene expression')
xlabel('brain regions')
ylabel('brain regions')
%% Gene by region for-loop
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
%region than another, visually looked at braintable to confirm this result
%region than another, visually looked at braintable to confirm this result