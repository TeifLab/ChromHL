%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Subroutine MapOfBindingCalc calculates the map ofbinding
% for a given c0(g), g=1,2,3
% ...........dStatda to be added later..............
%
% Stat - statsum
% max - after Trace of the current statsum value reaches max, StatSum
%       is being normilized by max
% NormCount - number of normalizations
% noLeftOverhang=.true. prohibits ligand overhang from the left DNA end
% noRightOverhang=.true. prohibits ligand overhang from the right DNA end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [c,teta] = PointOfMapOfBindingCalc(TestSiteNumber)

global fNumberOfLigands rank eNumberOfChromatinStates Lpolymer

% Boundary conditions:

LeftVect=ones([1 rank]);
RightVect=ones([rank 1]);
LeftVectdK=zeros([fNumberOfLigands 1 rank]);
RightVectdK=zeros([fNumberOfLigands rank 1]);
LeftVectdE=zeros([eNumberOfChromatinStates 1 rank]);
RightVectdE=zeros([eNumberOfChromatinStates rank 1]);

c = zeros([fNumberOfLigands+1 1]);
teta = zeros([eNumberOfChromatinStates 1]);

% Initialization:

StatSum = 0;
dStatdK = zeros([fNumberOfLigands 1]);
dStatdE = zeros([eNumberOfChromatinStates 1]);

% We have to calculate Partition Function for each ThisSiteNumber

for ThisSiteNumber = 1:Lpolymer
    
    %normalization:
    max=1.e25;
    NormCount=0;
    
    Trace1 = sum(LeftVect(1,:));
    
    if (StatSum>max || Trace1> max)
        LeftVectdK = LeftVectdK/max;
        LeftVectdE = LeftVectdE/max;
        LeftVect= LeftVect/max;
		StatSum = StatSum/max;
        NormCount = NormCount+1;
    end
    
    %StatSum calculation:
    
    [Q, dQdK, dQdE] = MatrixInitUnwrap(ThisSiteNumber,TestSiteNumber); %%CHECK THIS
    
    for g=1:fNumberOfLigands
        %LVDKg1z = squeeze(LeftVectdK(g,1,:))';
        %dQdKgzz = squeeze(dQdK(g,:,:));
		
		LVDKg1z = shiftdim(LeftVectdK(g,1,:))';
		dQdKgzz = shiftdim(dQdK(g,:,:));
		
        
        %disp(size(LVDKg1z))
        %disp(size(Q))
        
        LeftVectdK(g,1,1:rank) = LVDKg1z * Q ...
            +LeftVect(1,:) * dQdKgzz;
    end % g
    
    for e=1:eNumberOfChromatinStates
        %LVDEe1z = squeeze(LeftVectdE(e,1,:))';
        %dQdEgzz = squeeze(dQdE(e,:,:));
		
		LVDEe1z = shiftdim(LeftVectdE(e,1,:))';
        dQdEezz = shiftdim(dQdE(e,:,:));
		
        LeftVectdE(e,1,1:rank) = LVDEe1z * Q ...
            +LeftVect(1,:) * dQdEezz;
    end
    
    LeftVect=LeftVect * Q;
    
    
    %for debug----------------------------
    %     if (ThisSiteNumber==TestSiteNumber) then
    %     	fprintf(fid7,'n= %s\n',ThisSiteNumber);
    % 	fprintf(fid7,'%s\n','Nonzero Q(i,j) elements:');
    % 	for i=1:rank
    % 	    for j=1:rank
    %     		if Q(i,j)~=0
    % 		    fprintf(fid7,'Q(%d,%d) %f\n',i,j,Q(i,j));
    % 		end
    % 	    end %j
    % 	end %i
    %
    % 	for i=1:rank
    % 	    for j=1:rank
    % 		if dQdK(i,j)~=0
    % 		   fprintf(fid7,'dQdK(%d,%d) %f\n',i,j,dQdK(i,j));
    % 		end
    % 	    end %j
    % 	end %i
    %
    %
    %     end %(ThisSiteNumber.eq.TestSiteNumber)
    
    % for debug----------------------------
    
    
end %ThisSiteNumber

%find statsum multiplying LeftVect by RightVect
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


%find the probability of protein binding
for g=1:fNumberOfLigands
    c(g+1)=dStatdK(g)/StatSum; % "derivatives" were not divided by K(n,g), so I don't multiply here %+1 as
end %g

%find the probability of chromatin state e
for e=1:eNumberOfChromatinStates
    teta(e)=dStatdE(e)/StatSum; % "derivatives" were not divided by s(e), so I don't multiply here
end %e

c(1)=1.-1./StatSum*(max^NormCount);

end

