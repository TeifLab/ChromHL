%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subroutine ConstantsInitUnwrap
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = ConstantsInitUnwrap

global fNumberOfLigands rank  eNumberOfChromatinStates R
global m mStart % ligand lengthes: integer, dimension (1:10) 
global nMaxGap % max interaction lengthes : integer, dimension (1:10) 
global MAXnMaxGap MINnMaxGap % max(nMaxGap), min(nMaxGap)
global iGap % states corresponding to g-j-g gaps integer, dimension (1:10000,0:10) 
global iLeftFreeEnd iRightFreeEnd % polymer end matrix state numbers
global noGaps 


%=================================================================
%-----------------------------------------------------------------
%-----------------------------------------------------------------
%-----------------------------------------------------------------
%-----------------------------------------------------------------
%defining the states for the first ligand-DNA contact and the matrix rank:


MAXnMaxGap = max(nMaxGap(:));
MINnMaxGap = min(nMaxGap(:));

mStart = zeros([fNumberOfLigands 1]);
iGap = zeros([MAXnMaxGap fNumberOfLigands+1]);

for g=1:fNumberOfLigands
    if (g==1)
        mStart(g)=1;
    end
    if (g>1)
        mStart(g)=mStart(g-1)+m(g-1);
    end
end % g
iLeftFreeEnd=mStart(fNumberOfLigands)+m(fNumberOfLigands);
iRightFreeEnd=iLeftFreeEnd+1;

if(noGaps)
    rank=iRightFreeEnd+1;
else
    rank=iRightFreeEnd;
    %g1-j-g2 gaps:
    for g=1:fNumberOfLigands
        for j=1:nMaxGap(g)
            if (nMaxGap(g)>0)
                rank=rank+1;
                iGap(j,g+1)=rank;
            end
        end %j
    end %g
    
    %first unit in a gap longer then the interaction length:
    rank=rank+1;
    iGap(1,1)=rank;
    
    %the other units in a gap longer then the interaction length:
    for j=2:MAXnMaxGap+1
        rank=rank+1;
        iGap(j,1)=rank;
    end
    
    
end %(noGaps)

R = rank;
rank = eNumberOfChromatinStates*rank;

end %subroutine ConstantsInitUnwrap



