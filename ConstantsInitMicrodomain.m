% This function sets the numbers for the states in the individual blocks
% in the matrix (of size RxR), for the block-structured transfer matrix
% which has overall size (R * #states) x R * #states)

function [] = ConstantsInitMicrodomain

global fNumberOfLigands rank  eNumberOfChromatinStates R
global m mStart % ligand lengthes: integer, dimension (1:10) 
global nMaxGap % max interaction lengthes : integer, dimension (1:10) 
global MAXnMaxGap MINnMaxGap % max(nMaxGap), min(nMaxGap)
global iGap % states corresponding to g-j-g gaps integer, dimension (1:10000,0:10) 
global iLeftFreeEnd iRightFreeEnd % polymer end matrix state numbers
global noGaps 

% Set the maximum and minimum allowedgaps between proteins
MAXnMaxGap = max(nMaxGap(:));
MINnMaxGap = min(nMaxGap(:));

% Count states, mStart holds where each ligand binding states start
mStart = zeros([fNumberOfLigands 1]);
iGap = zeros([MAXnMaxGap fNumberOfLigands+1]);

% First enumerate the ligands
for g=1:fNumberOfLigands
    if (g==1)
        mStart(g)=1;
    end
    if (g>1)
        mStart(g)=mStart(g-1)+m(g-1);
    end
end % g

% Then set states for left free end and right free end
iLeftFreeEnd=mStart(fNumberOfLigands)+m(fNumberOfLigands);
iRightFreeEnd=iLeftFreeEnd+1;

% Then go through each potential gap between proteins and enumerate
% those states
if(noGaps) % If there are no gaps, only one more state which is a "not-bound" state
    rank=iRightFreeEnd+1;
else % enumerate gaps between proteins
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

% Set R, number of states per chromatin state and
% rank, the total number of states (which is R * # chromatin states)
R = rank;
rank = eNumberOfChromatinStates*rank;

end 
