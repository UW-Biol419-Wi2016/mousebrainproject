%Part 1 - setting up github

%% reading data
braintable = readtable('/Users/esteligarcia/GitHub/mousebrainproject/Data/BrainRegionTotalDataset_Log2FoldChange.xlsx');
diseasetable = readtable('/Users/esteligarcia/GitHub/mousebrainproject/Data/MGIdisease.txt',...
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
%% Huntington's

huntingtonsgenes = table([], [], 'VariableNames', {'GeneName', 'DiseaseInfo'});

findhunt = find(strncmpi(bdt{:, 2}, 'Huntington', 10));
%indices within bdt for parkinson disease genes

    %adding the parkinson genes to their table
        GeneName = bdt{findhunt, 1};
        %calls gene name (column 1) at all findpark locations
        
        DiseaseInfo = bdt{findhunt, 2};
        %calls disease info (column 5) at same indices 
   
        huntingtonsgenes = [huntingtonsgenes; table(GeneName, DiseaseInfo)];
        %adds rows to the table, bdt, carrying the values from above


%% Alzheimer's

alzgenes = table([], [], 'VariableNames', {'GeneName', 'DiseaseInfo'});

findalz = find(strncmpi(bdt{:, 2}, 'Alzheimer', 9));
%indices within bdt for parkinson disease genes

    %adding the parkinson genes to their table
        GeneName = bdt{findalz, 1};
        %calls gene name (column 1) at all findpark locations
        
        DiseaseInfo = bdt{findalz, 2};
        %calls disease info (column 5) at same indices 
   
        alzgenes = [alzgenes; table(GeneName, DiseaseInfo)];
        %adds rows to the table, bdt, carrying the values from above

%% Major Depressive Disorder

depressiongenes = table([], [], 'VariableNames', {'GeneName', 'DiseaseInfo'});

finddepression = find(strncmpi(bdt{:, 2}, 'Major Depressive', 16));
%indices within bdt for parkinson disease genes

    %adding the parkinson genes to their table
        GeneName = bdt{finddepression, 1};
        %calls gene name (column 1) at all findpark locations
        
        DiseaseInfo = bdt{finddepression, 2};
        %calls disease info (column 5) at same indices 
   
        depressiongenes = [depressiongenes; table(GeneName, DiseaseInfo)];
        %adds rows to the table, bdt, carrying the values from above
        
        %ONLY 2 genes, shouldn't do
        
%% Autism

autismgenes = table([], [], 'VariableNames', {'GeneName', 'DiseaseInfo'});

findaut = find(strncmpi(bdt{:, 2}, 'autism', 6));
%indices within bdt for parkinson disease genes

    %adding the parkinson genes to their table
        GeneName = bdt{findaut, 1};
        %calls gene name (column 1) at all findpark locations
        
        DiseaseInfo = bdt{findaut, 2};
        %calls disease info (column 5) at same indices 
   
        autismgenes = [autismgenes; table(GeneName, DiseaseInfo)];
        %adds rows to the table, bdt, carrying the values from above
        
        %Holy moly, there are 35 autism-linked genes!!! Great!
        
%% Schizophrenia

schizogenes = table([], [], 'VariableNames', {'GeneName', 'DiseaseInfo'});

findschizo = find(strncmpi(bdt{:, 2}, 'schizophrenia', 13));
%indices within bdt for parkinson disease genes

    %adding the parkinson genes to their table
        GeneName = bdt{findschizo, 1};
        %calls gene name (column 1) at all findpark locations
        
        DiseaseInfo = bdt{findschizo, 2};
        %calls disease info (column 5) at same indices 
   
        schizogenes = [schizogenes; table(GeneName, DiseaseInfo)];
        %adds rows to the table, bdt, carrying the values from above
        
        
        %33 genes for this one!! Wow!

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

%initializes matrix for gene by region data

%'brain stem', 'Cerebellum', 'Corpus Callosum', 'Motor Cortex', 'Olfactory Bulb', 'Optic Nerve', 'Prefrontal Cortex', 'Striatum', 'Thalamus', 'Hippocampus'});
%1:10 names of new double's columns
parkinsonsmat = genebyregionmaker(parkinsonsgenes, braintable, brainregionmat);
%insert other genes here

schizomat = genebyregionmaker(schizogenes, braintable, brainregionmat);

autmat = genebyregionmaker(autismgenes, braintable, brainregionmat);

huntmat = genebyregionmaker(huntingtonsgenes, braintable, brainregionmat);

alzmat = genebyregionmaker(alzgenes, braintable, brainregionmat);

  
%% Making random shuffles of genes to make sets of genebyregion matrices (for control)
 n = 10 ;
    covvalcell = cell(n, 1) ;
    %makes a cell matrix to store all matrices from the loop

%initalizing vectors to hold covariance values for each potential
%combination of brain regions(numbers here are the labels for each region
%as specified below)
%'Cerebellum', 'Corpus Callosum', 'Motor Cortex', 'Olfactory Bulb', 'Optic Nerve', 'Prefrontal Cortex', 'Striatum', 'Thalamus', 'Hippocampus'});
%1:10 names of new double's columns
controlgenebyregion = [];
numcontrols = 1000;
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
    

    %covariance genebygene
    
    covvalcell{i} = controlcov; 
    %stores each successive controlcov in a cell for later analysis
    
end;
    


%visually looked at braintable to confirm this result

%% All cells for diseases
autcell = covval(autmat, 1000, braintable, brainregionmat, bdt);
%makes controls for genebygene autism mat

alzcell = covval(alzmat, 1000, braintable, brainregionmat, bdt);

huntcell = covval(huntmat, 1000, braintable, brainregionmat, bdt);

schizocell = covval(schizomat, 1000, braintable, brainregionmat, bdt);

parkcell = covval(parkinsonsmat, 1000, braintable, brainregionmat, bdt);

%% Same as below (obtaining histograms and others) but from gene by gene cov matrices
[autmean, autstd] = cellmeanfinder(autmat, autcell, 1000);
%finds mean and std for each combo for autism

[alzmean, alzstd] = cellmeanfinder(alzmat, alzcell, 1000);

[parkmean, parkstd] = cellmeanfinder(parkinsonsmat, parkcell, 1000);

[schizomean, schizostd] = cellmeanfinder(schizomat, schizocell, 1000);

[huntmean, huntstd] = cellmeanfinder(huntmat, huntcell, 1000);

%% comparing cov matrix to control cov matrices
geneparkinsonscov = cov(parkinsonsmat');
geneparkinsonsdist = distcomp(geneparkinsonscov, parkmean, parkstd);

genealzcov = cov(alzmat');
genealzdist = distcomp(genealzcov, alzmean, alzstd);


genehuntcov = cov(huntmat');
genehuntdist = distcomp(genehuntcov, huntmean, huntstd);

geneautcov = cov(autmat');
geneautdist = distcomp(geneautcov, autmean, autstd);

geneschizocov = cov(schizomat');
geneschizodist = distcomp(geneschizocov, schizomean, schizostd);

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


%% Calculating std from mean for parkinsons/diseases vs controls
%parkinsonsmat
%controlstd
%controlmean
parkinsonsdist = [];
parkinsonscov = cov(parkinsonsmat);
alzdist = [];
alzcov = cov(alzmat);
huntdist = [];
huntcov = cov(huntmat);
autdist = [];
autcov = cov(autmat);
schizodist = [];
schizocov = cov(schizomat);


for i= 1:10
    for j = 1:i
    parkinsonsdist(i, j) = (parkinsonscov(i, j)- controlmean(i, j))/controlstd(i, j);
    %for each value of parkinsons covariance, subtract the mean from that
    %value and divide by the std for that value to get the mahalanobis
    %distance for each brain region 
    alzdist(i, j) = (controlmean(i, j) - alzcov(i,j))/controlstd(i, j);
    autdist(i, j) = (controlmean(i, j) - autcov(i,j))/controlstd(i, j);
    schizodist(i, j) = (controlmean(i, j) - schizocov(i, j))/controlstd(i, j);
    huntdist(i, j) = (controlmean(i, j) - huntcov(i, j))/controlstd(i, j);
    end;
end;


%% Visualizing the data results (dist matrices)
figure;
hold on
imagesc(parkinsonsdist);
colorbar
title('1A Parkinsons Gene Expression by Region');
xlabel('Brain Region');
ylabel('Brain Region');
textStrings = num2str(parkinsonsdist(:),'%0.2f'); %make strings from matrix
textStrings = strtrim(cellstr(textStrings)); % trims
[x,y] = meshgrid(1:length(parkinsonsdist));   %# Create x and y coordinates for the strings
hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
              'HorizontalAlignment','center');

figure;
imagesc(geneparkinsonsdist);
colorbar
title('1B Parkinsons Gene Expression by Gene');
xlabel('Genes');
ylabel('Genes');
textStrings = num2str(geneparkinsonsdist(:),'%0.2f'); %make strings from matrix
textStrings = strtrim(cellstr(textStrings)); % trims
[x,y] = meshgrid(1:length(geneparkinsonsdist));   %# Create x and y coordinates for the strings
hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
              'HorizontalAlignment','center');

figure;
imagesc(alzdist);
colorbar
title('2A Alzheimers Gene Expression by Region');
xlabel('Brain Region');
ylabel('Brain Region');
textStrings = num2str(alzdist(:),'%0.2f'); %make strings from matrix
textStrings = strtrim(cellstr(textStrings)); % trims
[x,y] = meshgrid(1:length(alzdist));   %# Create x and y coordinates for the strings
hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
              'HorizontalAlignment','center');

figure;
imagesc(genealzdist);
colorbar
title('2B Alzheimers Gene Expression by Gene');
xlabel('Genes');
ylabel('Genes');
textStrings = num2str(genealzdist(:),'%0.2f'); %make strings from matrix
textStrings = strtrim(cellstr(textStrings)); % trims
[x,y] = meshgrid(1:length(genealzdist));   %# Create x and y coordinates for the strings
hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
              'HorizontalAlignment','center');


figure;
imagesc(huntdist);
colorbar
title('3A Huntingtons Gene Expression by Region');
xlabel('Brain Region');
ylabel('Brain Region');
textStrings = num2str(huntdist(:),'%0.2f'); %make strings from matrix
textStrings = strtrim(cellstr(textStrings)); % trims
[x,y] = meshgrid(1:length(huntdist));   %# Create x and y coordinates for the strings
hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
              'HorizontalAlignment','center');

figure;
imagesc(genehuntdist);
colorbar
title('3B Huntingtons Gene Expression by Genes');
xlabel('Genes');
ylabel('Genes');
textStrings = num2str(genehuntdist(:),'%0.2f'); %make strings from matrix
textStrings = strtrim(cellstr(textStrings)); % trims
[x,y] = meshgrid(1:length(genehuntdist));   %# Create x and y coordinates for the strings
hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
              'HorizontalAlignment','center');

figure;
imagesc(autdist);
colorbar
title('4A Autism Gene Expression by Region');
xlabel('Region');
ylabel('Region');
textStrings = num2str(autdist(:),'%0.2f'); %make strings from matrix
textStrings = strtrim(cellstr(textStrings)); % trims
[x,y] = meshgrid(1:length(autdist));   %# Create x and y coordinates for the strings
hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
              'HorizontalAlignment','center');

figure;
imagesc(geneautdist);
colorbar
title('4B Autism Gene Expression by Genes');
xlabel('Genes');
ylabel('Genes');
textStrings = num2str(geneautdist(:),'%0.2f'); %make strings from matrix
textStrings = strtrim(cellstr(textStrings)); % trims
[x,y] = meshgrid(1:length(geneautdist));   %# Create x and y coordinates for the strings
hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
              'HorizontalAlignment','center');

figure;
imagesc(schizodist);
colorbar
title('5A Schizophrenia Gene Expression by Region');
xlabel('Region');
ylabel('Region');
textStrings = num2str(schizodist(:),'%0.2f'); %make strings from matrix
textStrings = strtrim(cellstr(textStrings)); % trims
[x,y] = meshgrid(1:length(schizodist));   %# Create x and y coordinates for the strings
hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
              'HorizontalAlignment','center');

figure;
imagesc(geneschizodist);
colorbar
title('5B Schizophrenia Gene Expression by Gene');
xlabel('Genes');
ylabel('Genes');
textStrings = num2str(geneschizodist(:),'%0.2f'); %make strings from matrix
textStrings = strtrim(cellstr(textStrings)); % trims
[x,y] = meshgrid(1:length(geneschizodist));   %# Create x and y coordinates for the strings
hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
              'HorizontalAlignment','center');