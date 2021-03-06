function [ ] = fullProtein( header, outFile, lambda)
% this script runs full protein optimization of self and pair parameters

tstart = cputime;

self = 20;
pair = 400;

if isempty(lambda)
    lambda = 1000;
end

reslist = 'residue.txt';
cmap = 'contact.txt';

EAAfreq = importdata('aafreq.txt');
EAAfreq = -log(EAAfreq);
EAAfreq = EAAfreq - mean(EAAfreq(:));

% load a residue list
residueList = textscan(fopen(sprintf('%s', reslist)), '%s');
nResidues = length(residueList{1});

% load a contact map, and get all parameters to optimize
contactList = textscan(fopen(sprintf('%s', cmap)), repmat('%s ', [1, 2]));
nContacts = length(contactList{1});


%% first somehow determine an order of residues
% intuitively, the residues with more contacts should be looked first
residueNcon = zeros(nResidues,1);
for r = 1:length(residueList{1})
    Ncon = length(find(strcmp(contactList{1}, residueList{1}(r)) | strcmp(contactList{2}, residueList{1}(r))));
    residueNcon(r) = Ncon;
end
residueList = residueList{1};

[~, orderedResidues] = sort(residueNcon, 1, 'descend');

%% iterate until converge
Conv_Threshold = 0.01*nResidues;

currentParamsValues = cell(nResidues + nContacts, 1);
for i = 1: nResidues
    currentParamsValues{i} = zeros(1, self);
end
for i = 1: size(contactList{1}, 1)
    currentParamsValues{nResidues + i} = zeros(1, pair);
end

pobjTotal = Inf;

while (1)
    %% iterate over residues
    objTotal = 0;
    for i = 1: nResidues
        currentResidue = residueList{orderedResidues(i)};
        conInds = find( strcmp(contactList{1}, currentResidue) | strcmp(contactList{2}, currentResidue));
        usedCons = [contactList{1}(conInds), contactList{2}(conInds)];
        conNames = usedCons(~strcmp(usedCons, currentResidue));
        
        % load data for this position
        leftMat = [];
        try
            [leftMat, rightVec, paramsUsage] = loadEnsemble(currentResidue, conNames, header); 
        catch e
            disp('');
        end
        if isempty(leftMat)
            continue;
        end
        % prior background values;
        EAAfreqVec = repmat(EAAfreq, size(rightVec, 1)/20, 1);
        
        % optimize parameters for this position, should return the value of
        % objective function and current parameters
        % fitModel maybe kept the same
%         defaultValue = zeros(self + length(conNames) * pair, 1);
        
        % true if the current residue is the first one in the pair
        firstInPair = logical(strcmp(usedCons(:,1), currentResidue));
        
        % get current values, can be results from previous round
        currentValue = zeros(1, self + pair*length(conInds));
        currentValue(1:self) = currentParamsValues{orderedResidues(i)};
        for j = 1:length(conNames)
            currentPairValue = currentParamsValues{nResidues + conInds(j)};
            if ~firstInPair(j) % if currentResidue is not the first one in pair
                currentPairValue = reshape(reshape(currentPairValue, self, self)', 1, pair);
            end
            currentValue( self+1+(j-1)*pair : self+j*pair ) = currentPairValue;
        end
            
        [optParams, obj] = ffitModel(leftMat, rightVec, currentValue', currentValue, EAAfreqVec, paramsUsage, lambda);
        
        currentParamsValues{orderedResidues(i)} = optParams(1:self)';
        for j = 1:length(conNames)
            optParamsforJ = optParams(self + (j-1)*pair + 1: self + j*pair)';
            if ~firstInPair(j)
                optParamsforJ = reshape(reshape(optParamsforJ, self, self)', 1, pair);
            end  
            currentParamsValues{nResidues + conInds(j)} = optParamsforJ';
        end
        objTotal = objTotal + obj;
    end
    save(strcat(outFile, '.mat'), 'currentParamsValues');
    
    disp(sprintf('objective from previous grand cycle is %f, and now is %f', pobjTotal, objTotal));

    if pobjTotal - objTotal < Conv_Threshold
        break;
    end
    pobjTotal = objTotal;
end

telapse = cputime - tstart;
disp(telapse);

end

