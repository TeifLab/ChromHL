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

function[cMap,tetaMap] = MapOfBindingCalc()

global Lpolymer fNumberOfLigands eNumberOfChromatinStates
global fid7

cMap=zeros([fNumberOfLigands+1  Lpolymer]);
tetaMap=zeros([eNumberOfChromatinStates  Lpolymer]);

TIME=tic;

%write(7,*) ' Calculating the map of binding: g=', g
for TestSiteNumber=1:Lpolymer
    %for g=1:fNumberOfLigands
    %    for e=1:eNumberOfChromatinStates
            
            fprintf('(site,g,e)={%d,%s,%s}\n',TestSiteNumber,'all','all');
            
            [c,teta]=PointOfMapOfBindingCalc(TestSiteNumber);
            
            cMap(2:end,TestSiteNumber)=c(2:end);
            tetaMap(:,TestSiteNumber)=teta(:);
            
            
            % 			%stiff boundary conditions
            % 			if (TestSiteNumber+m(g)-1.le.Lpolymer) then
            %
            % 				for i=0,m(g)-1
            % 					cMap(g, TestSiteNumber+i)=cMap(g, TestSiteNumber+i)+c(g)
            % 				end %i
            % 			end
            %
            
            
            %time estimation
            if (TestSiteNumber==1 )%&& g==1 && e==1)
                time2 = toc(TIME);
%                fprintf(fid7,'%s\n',['One point:' num2str(time2) 'sec']);
%                fprintf(fid7,'%s\n',['Estimated calculation time:' ...
%                    num2str(time2*Lpolymer/60) ' min']);
            end %time estimation
            
      %  end % g
    %end %e
end % TestSiteNumber



end

