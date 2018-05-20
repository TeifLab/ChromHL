% This function is called to calculate the 
% protein binding and chromatin state probability at
% TestSiteNumber

function [c,teta] = PointOfMapOfBindingCalc(TestSiteNumber)

global fNumberOfLigands rank eNumberOfChromatinStates Lpolymer

% Set the boundary conditions for the vectors
% LeftVect and RightVect and the derivatives
% with respect to K (for protein concentration)
% and with respect to E (for the chromatin state)
LeftVect=ones([1 rank]);
RightVect=ones([rank 1]);
LeftVectdK=zeros([fNumberOfLigands 1 rank]);
RightVectdK=zeros([fNumberOfLigands rank 1]);
LeftVectdE=zeros([eNumberOfChromatinStates 1 rank]);
RightVectdE=zeros([eNumberOfChromatinStates rank 1]);

% These will hold the point of binding map
c = zeros([fNumberOfLigands+1 1]);
teta = zeros([eNumberOfChromatinStates 1]);

% Initialize the variables:
StatSum = 0;
dStatdK = zeros([fNumberOfLigands 1]);
dStatdE = zeros([eNumberOfChromatinStates 1]);

% We have to calculate the partition function Z
% and the derivatives dZ/dK and dZ/dE for ther
% binding map 
for ThisSiteNumber = 1:Lpolymer
    
    % if the function gets above this threshold
    % then normalise by this and count how
    % often this happens
    % this helps with numerical stability of the algorithm
    max=1.e25;
    NormCount=0;
    
    % Calculate the L1 norm of the LeftVector
    Trace1 = sum(LeftVect(1,:));
    
    % If we are above the threshold (either in the vector, or in the 
    % partition function itself), then normalise
    if (StatSum>max || Trace1> max)
        LeftVectdK = LeftVectdK/max;
        LeftVectdE = LeftVectdE/max;
        LeftVect= LeftVect/max;
		StatSum = StatSum/max;
        NormCount = NormCount+1;
    end
    
    % Get the transfer matrix and the derivatives of the transfer matrix
    [Q, dQdK, dQdE] = MatrixInitMicrodomain(ThisSiteNumber,TestSiteNumber);
    
    % For each ligand
    for g=1:fNumberOfLigands
        
	% generate some temporary vectors to ensure the indices are correct
	% in the below multiplication
	LVDKg1z = shiftdim(LeftVectdK(g,1,:))';
	dQdKgzz = shiftdim(dQdK(g,:,:));
		
        % update the derivative of the vector (d/dK) for each ligand
        LeftVectdK(g,1,1:rank) = LVDKg1z * Q ...
            +LeftVect(1,:) * dQdKgzz;
    end % g
    
    % for each chromain state
    for e=1:eNumberOfChromatinStates
        
	% generate the temp vectors again to ensure correct indices in the
	% below multiplication
	
	LVDEe1z = shiftdim(LeftVectdE(e,1,:))';
        dQdEezz = shiftdim(dQdE(e,:,:));
	
	% update the derivative of the vector (d/dE) for each chromatin state
        LeftVectdE(e,1,1:rank) = LVDEe1z * Q ...
            +LeftVect(1,:) * dQdEezz;
    end
    
    LeftVect=LeftVect * Q;
    
end %ThisSiteNumber

% Find the partition function and derivatives
% by multiplying leftVect by RightVect
% and updating function
% also update dZ/dK 
% and dZ/dE
for i=1:rank
    for g=1:fNumberOfLigands
        dStatdK(g)=dStatdK(g) + LeftVectdK(g,1,i)*RightVect(i,1) ...
            +LeftVect(1,i)*RightVectdK(g,i,1);
    end
    
    for e=1:eNumberOfChromatinStates
        dStatdE(e)=dStatdE(e) + LeftVectdE(e,1,i)*RightVect(i,1) ...
            +LeftVect(1,i)*RightVectdE(e,i,1);
    end
    
    
    StatSum=StatSum + LeftVect(1,i)*RightVect(i,1);
end %i

% Use partition function to calculate probability of protein binding
for g=1:fNumberOfLigands
    c(g+1)=dStatdK(g)/StatSum; 
end %g

% Use partition function to calculate probability of chromatin state
for e=1:eNumberOfChromatinStates
    teta(e)=dStatdE(e)/StatSum; 
end %e

% Set unbound protein binding probability
c(1)=1.-1./StatSum*(max^NormCount);

end

