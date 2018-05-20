% This calculates the binding map by calling PointOfMapOfBindingCalc 
% for each unit in the lattice, returning it in cMap and tetaMap
%
% Note that this and PointOfMapOfBindingCalc can be sped up by
% pre-calculating the transfer matrices and derivatives

function[cMap,tetaMap] = MapOfBindingCalc()

global Lpolymer fNumberOfLigands eNumberOfChromatinStates

% These variables will hold the protein binding map and
% chromatin state probability map
cMap=zeros([fNumberOfLigands+1  Lpolymer]);
tetaMap=zeros([eNumberOfChromatinStates  Lpolymer]);


% Work along the polymer from left-end
for TestSiteNumber=1:Lpolymer
    
    % Output to stdout which site we are currently working with
    fprintf('(site,g,e)={%d,%s,%s}\n',TestSiteNumber,'all','all');
            
    % Calculate the binding at this site
    [c,teta]=PointOfMapOfBindingCalc(TestSiteNumber);
     
    % fill in the binding map point at this site 
    cMap(2:end,TestSiteNumber)=c(2:end);
    tetaMap(:,TestSiteNumber)=teta(:);
    
end 

end
