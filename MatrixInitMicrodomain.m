% This function calculates the matrix for the i'th polymer unit
%
% Q(rank,rank)                            : is the matrix of statistical weights
% dQdK(fNumberOfLigands,rank,rank)        : is the matrix of derivatives of Q by c0(g), g=1..fNumberOfLigands
% dQdE(eNumberOfChromatinStates,rank,rank): is the matrix of derivatives of Q by s, f=1..eNumberOfChromatinStates
%
% c0(g) is the bulk molar concentration of g-type ligand in solution
% ThisSiteNumber is the number of the site we calculate the matrix
% TestSiteNumber if the number of the site for calculating the derivatives
% KKK is the combined binding constant
% g1,g2 ligand type
% h1,h2 number of unwrapped protein units
% e1,e2 chromatin state

function [Q,dQdK,dQdE] = MatrixInitMicrodomain(ThisSiteNumber, TestSiteNumber)

global Lpolymer fNumberOfLigands rank eNumberOfChromatinStates R
global m mStart % ligand lengthes: integer, dimension (1:10)
global nMaxGap % max interaction lengthes : integer, dimension (1:10)
global MAXnMaxGap MINnMaxGap % max(nMaxGap), min(nMaxGap)
global w %cooperativity parameters : double precision w(0:10000,0:10,0:10, 0:5)
global KKK % double precision(10000,10,147,5)
global Unwrap % double precision Unwrap(0:200,1:10)
global  c0 % c0(10)
global s sigma % s(1:10000,1:5), sigma(1:5,1:5)
global iGap % states corresponding to g-j-g gaps integer, dimension (1:10000,0:10)
global iLeftFreeEnd iRightFreeEnd % polymer end matrix state numbers
global noGaps
global lig1mod lig2mod lig3mod lig4mod lig5mod lig6mod

% Check if ligands are modified
if(lig2mod)
    c0(2)=c0(1);
end
if(lig3mod)
    c0(2)=c0(1);
    c0(3)=c0(1);
end
if(lig4mod)
    c0(2)=c0(1);
    c0(3)=c0(1);
    c0(4)=c0(1);
end
if(lig5mod)
    c0(2)=c0(1);
    c0(3)=c0(1);
    c0(4)=c0(1);
    c0(5)=c0(1);
end
if(lig6mod)
    c0(2)=c0(1);
    c0(3)=c0(1);
    c0(4)=c0(1);
    c0(5)=c0(1);
    c0(6)=c0(1);
end

% generate some aliases and initialise transfer matrix Q
% and derivatives dQdK, dQdE
n=ThisSiteNumber;
f=fNumberOfLigands;
eMax=eNumberOfChromatinStates;
Q = zeros(rank,rank);
dQdK = zeros(fNumberOfLigands,rank,rank);
dQdE = zeros(eNumberOfChromatinStates,rank,rank);

% Go through each case and compute the non-zero elements of the transfer matrix
% Qij, and its derivatives dQdK and dQdE

% 1) 1st unit of g-type protein followed by the 2nd unit:
for e1=1:eMax
    for e2=1:eMax
        for g1=1:f
            if (m(g1)>1)
                if (n<=(Lpolymer - (m(g1) - 1)))
                    i=R*(e1-1)+mStart(g1);
                    j=R*(e2-1)+i+1;
                    if (n>1)
                        Q(i,j) = s(n,e1)*sigma(e1,e2)*KKK(n,g1,1,e1);
                    end
                    if (n==1)
                        Q(i,j) = s(n,e1)*sigma(e1,e2)*KKK(n,g1,1,e1)*c0(g1);
                    end
                    if (n==TestSiteNumber)
                        dQdK(g1,i,j) = Q(i,j);
                        dQdE(e1,i,j) = Q(i,j);
                    end
                end %ThisSiteNumber
            end %m(g1)
        end %g1
    end
end


% 2) Unit h+1 of g-type protein followed by unit  h+2:
for e1=1:eMax
    for e2=1:eMax
        for g1=1:f
            if (m(g1)>2)
                for h1=1:(m(g1)-2)
                    if (n>h1 && n<=(Lpolymer - (m(g1) - h1-1)))
                        i=R*(e1-1)+mStart(g1)+h1;
                        j=R*(e2-1)+i+1;
                        Q(i,j) = s(n,e1)*KKK(n,g1,h1+1,e1);
                        if (n==TestSiteNumber)
                            dQdK(g1,i,j) = Q(i,j);
                            dQdE(e1,i,j) = Q(i,j);
                        end
                    end % ThisSiteNumber
                end % h
            end % (m(g1)>2)
        end % g1
    end
end


% 3) Unit (mg - h1) of g-type protein followed by a right free DNA end:
for e1=1:eMax
    for e2=1:eMax
        for g1=1:f
            for h1=0:(m(g1)-1)
                if (n>=m(g1)-h1 && n<=Lpolymer-h1)
                    i=R*(e1-1)+mStart(g1)+m(g1)-1-h1;
                    j=R*(e2-1)+iRightFreeEnd;
                    if (m(g1)>1)
                        if(n>1)
                            Q(i,j) = s(n,e1)*sigma(e1,e2)*KKK(n,g1,m(g1)-h1,e1)*Unwrap(h1+1,g1);
                        end
                        if(n==1)                           
                            Q(i,j) = s(n,e1)*sigma(e1,e2)*KKK(n,g1,m(g1)-h1,e1)*Unwrap(h1+1,g1)*c0(g1);
                        end
                        if (n==TestSiteNumber)
                            dQdK(g1,i,j) = Q(i,j);
                            dQdE(e1,i,j) = Q(i,j);
                        end
                    end
                    if (m(g1)==1)
                        if (n>1)
                            Q(i,j) = s(n,e1)*sigma(e1,e2)*KKK(n,g1,1,e1)*w(1,g1+1,1,e1+1);
                        end
                        if (n==1)
                            Q(i,j) = s(n,e1)*sigma(e1,e2)*KKK(n,g1,1,e1)*w(1,g1+1,1,e1+1)*c0(g1);
                        end
                        if (n==TestSiteNumber)
                            dQdK(g1,i,j) = Q(i,j);
                            dQdE(e1,i,j) = Q(i,j);
                        end
                    end
                end %n
            end %h
        end %g
    end
end


% 4) Unit (mg1 - h1) of g1-protein followed by unit 1 of g2-protein (no gap):

for e1=1:eMax
    for e2=1:eMax
        for g1=1:f
            for g2=1:f
                for h1=0:m(g1)-1
                    i=R*(e1-1)+mStart(g1)+m(g1)-1-h1;
                    j=R*(e2-1)+mStart(g2);
                    if (n>=(m(g1)-h1) && n<=Lpolymer-m(g2))
                        if (m(g1)>1)
                            if(n>1)
                                Q(i,j) = s(n,e1)*sigma(e1,e2)*w(1,g1+1,g2+1,e1+1)*KKK(n,g1,m(g1)-h1,e1)*Unwrap(h1+1,g1)*c0(g2);
                            end
                            if(n==1)
                                Q(i,j) = s(n,e1)*sigma(e1,e2)*w(1,g1+1,g2+1,e1+1)*KKK(n,g1,m(g1)-h1,e1)*Unwrap(h1+1,g1)*c0(g2)*c0(g1);
                            end
                            if (n==TestSiteNumber)
								dQdK(g1,i,j) = Q(i,j);
								dQdE(e1,i,j) = Q(i,j);
                            end
                        end
                        if (m(g1)==1)
                            if(n>1)
                                Q(i,j) = s(n,e1)*sigma(e1,e2)*KKK(n,g1,1,e1)*w(1,g1+1,g2+1,e1+1)*c0(g2);
                            end
                            if(n==1)
                                Q(i,j) = s(n,e1)*sigma(e1,e2)*KKK(n,g1,1,e1)*w(1,g1+1,g2+1,e1+1)*c0(g2)*c0(g1);
                            end
                            if (n==TestSiteNumber)
								dQdK(g1,i,j) = Q(i,j);
								dQdE(e1,i,j) = Q(i,j);
                            end
                        end % m(g1)
                    end %n
                end %h1
            end %g2
        end %g1
    end
end


% 5) Last unit of g1-protein followed by unit h2+1 of g2-protein (no gap):
for e1=1:eMax
    for e2=1:eMax
        for g1=1:f
            for g2=1:f
                for h2=1:m(g2)-1
                    i=R*(e1-1)+mStart(g1)+m(g1)-1;
                    j=R*(e2-1)+mStart(g2)+h2;
                    if (n>=m(g1) && n<=Lpolymer-(m(g2)-h2))
                        if (m(g1)>1)
                            Q(i,j) = s(n,e1)*sigma(e1,e2)*w(1,g1+1,g2+1,e1+1)*KKK(n,g1,m(g1),e1)*Unwrap(h2+1,g2)*c0(g2);
                            if (n==TestSiteNumber)
                                dQdK(g1,i,j) = Q(i,j);
                                dQdE(e1,i,j) = Q(i,j);
                            end
                        end
						
                        if (m(g1)==1)
                            if(n>1)
                                Q(i,j) = s(n,e1)*sigma(e1,e2)*KKK(n,g1,1,e1)*w(1,g1+1,g2+1,e1+1)*c0(g2);
                            end
                            if(n==1)
                                Q(i,j) = s(n,e1)*sigma(e1,e2)*KKK(n,g1,1,e1)*w(1,g1+1,g2+1,e1+1)*c0(g2)*c0(g1);
                            end
                            if (n==TestSiteNumber)
								dQdK(g1,i,j) = Q(i,j);
								dQdE(e1,i,j) = Q(i,j);
                            end
                        end % m(g1)
                    end %n
                end %h2
            end %g2
        end %g1
    end
end


% 6) Left free end continues:
for e1=1:eMax
    for e2=1:eMax
        i=R*(e1-1)+iLeftFreeEnd;
        j=R*(e2-1)+iLeftFreeEnd;
        Q(i,j) = s(n,e1)*sigma(e1,e2)*w(1,1,1,e1+1);
        if (n==TestSiteNumber)
            dQdE(e1,i,j) = Q(i,j);
        end;
    end %e1
end %e2

% 7) Right free end continues:
for e1=1:eMax
    for e2=1:eMax
        i=R*(e1-1)+iRightFreeEnd;
        j=R*(e2-1)+iRightFreeEnd;
        if(n>1)
            Q(i,j) = s(n,e1)*sigma(e1,e2)*w(1,1,1,e1+1);
        end
        if (n==TestSiteNumber)
            dQdE(e1,i,j) = Q(i,j);
        end
    end %e1
end %e2

% 8) Left free end followed by unit h+1 of g-type protein:
for e1=1:eMax
    for e2=1:eMax
        for g2=1:f
            for h2=0:m(g2)-1
                i=R*(e1-1)+iLeftFreeEnd;
                j=R*(e2-1)+mStart(g2) + h2;
                if (n<=Lpolymer-(m(g2)-h2))
                    Q(i,j) = s(n,e1)*sigma(e1,e2)*w(1,1,g2+1,e2+1)*Unwrap(h2+1,g2)*c0(g2);
                    if (n==TestSiteNumber)
                        dQdE(e1,i,j) = Q(i,j);
                    end
                end %ThisSiteNumber
            end %h2
        end %g
    end
end


if(~noGaps)
    
    %non-interacting gaps:
    
    % 9) Unit (mg - h) of g-type protein followed by a non-interacting gap longer Vg:
    for e1=1:eMax
        for e2=1:eMax
            for g1=1:f
                for h1=0:m(g1)-1
                    if (n>=(m(g1)-h1) && n<Lpolymer)  % hanging out prohibited
                        i=R*(e1-1)+mStart(g1)+m(g1)-1-h1;
                        j=R*(e2-1)+iGap(1,1);
                        if (m(g1)>1)
                            if(n>1)
                                Q(i,j) = s(n,e1)*sigma(e1,e2)*w(nMaxGap(g1)+2,g1+1,1,e1+1)*KKK(n,g1,m(g1)-h1,e1)*Unwrap(h1+1,g1);
                            end
                            if(n==1)
                                Q(i,j) = s(n,e1)*sigma(e1,e2)*w(nMaxGap(g1)+2,g1+1,1,e1+1)*KKK(n,g1,m(g1)-h1,e1)*Unwrap(h1+1,g1)*c0(g1);
                            end
                            if (n==TestSiteNumber)
                                dQdK(g1,i,j) = Q(i,j);
                                dQdE(e1,i,j) = Q(i,j);
                            end
                        end
                        
                        if (m(g1)==1)
                            if(n>1)
                                Q(i,j) = s(n,e1)*sigma(e1,e2)*KKK(n,g1,1,e1)*w(nMaxGap(g1)+2,g1+1,1,e1+1);
                            end
                            if(n==1)
                                Q(i,j) = s(n,e1)*sigma(e1,e2)*KKK(n,g1,1,e1)*w(nMaxGap(g1)+2,g1+1,1,e1+1)*c0(g1);
                            end
                            if (n==TestSiteNumber)
                                dQdK(g1,i,j) = Q(i,j);
                                dQdE(e1,i,j) = Q(i,j);
                            end %(n==TestSiteNumber)
                        end %(m(g1)==1)
                    end %n
                end %h1
            end %g1
        end
    end
    
    % 10) Large non-interacting gap s, units before max(Vg):
    for e1=1:eMax
        for e2=1:eMax
            for j=1:MAXnMaxGap
                i = R*(e1-1)+iGap(j,1);
                jj= R*(e2-1)+iGap(j+1,1);            
                if (n>1 && n<Lpolymer)
                    Q(i,jj) = s(n,e1)*sigma(e1,e2);
                    if (n~=TestSiteNumber)
                        dQdE(e1,i,jj) = Q(i,jj);
                    end
                    
                end
            end % j
        end
    end
    
    
    % 11) Large non-interacting gap s, units after max(Vg):
    for e1=1:eMax
        for e2=1:eMax
            i=R*(e1-1)+iGap(MAXnMaxGap+1,1);
            j=R*(e2-1)+iGap(MAXnMaxGap+1,1);
            if (n>1 && n<Lpolymer)
                Q(i,j) = s(n,e1)*sigma(e1,e2);
                if (n==TestSiteNumber)
                    dQdE(e1,i,j) = Q(i,j);
                end
            end
        end
    end
    
    % 12) Large non-interacting gap followed by unit h+1 of g2-type protein:
    for e1=1:eMax
        for e2=1:eMax
            for g2=1:f
                for h2=0:m(g2)-1
                    for jj=MINnMaxGap+1: MAXnMaxGap+1
                        if(jj>nMaxGap(g2))
                            i=R*(e1-1)+iGap(jj,1);
                            j=R*(e2-1)+mStart(g2)+h2;
                            if (n>1 && n<=Lpolymer-(m(g2)-h2))
                                Q(i,j) = s(n,e1)*sigma(e1,e2)*w(nMaxGap(g2)+2,1,g2+1,e1+1)*Unwrap(h2+1,g2)*c0(g2);
                            end %n
                            if (n==TestSiteNumber)
                                dQdE(e1,i,j) = Q(i,j);
                            end
                        end %jj
                    end % jj
                end %h2
            end %g
        end
    end
    
    
    
    %the gaps covered by protein-protein interactions:
    
    if(nMaxGap(g2)<nMaxGap(g1))
        nMaxGap(g2)=nMaxGap(g1);
    end
    
    % 13) Unit (mg1 - h) of g1-type protein followed by g1-L-g2 gap:
    for e1=1:eMax
        for e2=1:eMax
            for g1=1:f
                for g2=1:f
                    for h1=0:m(g1)-1
                        for jj=1:nMaxGap(g2)
                            if (n>=(m(g1)-h1) && n<Lpolymer && nMaxGap(g2)>0)  % hanging out prohibited
                                i=R*(e1-1)+mStart(g1)+m(g1)-1-h1;
                                j=R*(e2-1)+iGap(jj,g2+1);
                                if (m(g1)>1)
                                    if(n>1)
                                        Q(i,j) = s(n,e1)*sigma(e1,e2)*w(jj+1,g1+1,g2+1,e1+1)*KKK(n,g1,m(g1)-h1,e1)*Unwrap(h1+1,g1);
                                    end
                                    if(n==1)
                                        Q(i,j) = s(n,e1)*sigma(e1,e2)*w(jj+1,g1+1,g2+1,e1+1)*KKK(n,g1,m(g1)-h1,e1)*Unwrap(h1+1,g1)*c0(g1);
                                    end
                                    if (n==TestSiteNumber)
                                        dQdK(g1,i,j)= Q(i,j);
                                        dQdE(e1,i,j)= Q(i,j);
                                    end
                                end %(m(g1)>1)
                                
                                if (m(g1)==1)
                                    if(n>1)
                                        Q(i,j) = s(n,e1)*sigma(e1,e2)*KKK(n,g1,1,e1)*w(jj+1,g1+1,g2+1,e1+1);
                                    end
                                    if(n==1)
                                        Q(i,j) = s(n,e1)*sigma(e1,e2)*KKK(n,g1,1,e1)*w(jj+1,g1+1,g2+1,e1+1)*c0(g1);
                                    end
                                    if (n==TestSiteNumber)
                                        dQdK(g1,i,j)= Q(i,j);
                                        dQdE(e1,i,j)= Q(i,j);
                                    end
                                end
                            end %n
                        end %h1
                    end %jj
                end %g2
            end %g1
        end
    end
    
    
    % 14) g1-L-g2 gap s (L free units before g2-type protein):
    for e1=1:eMax
        for e2=1:eMax
            for g2=1:f
                for j=2:nMaxGap(g2)
                    i = R*(e1-1)+iGap(j,g2+1);
                    jj =R*(e2-1)+iGap(j-1,g2+1);
                    if (n>1 && n<Lpolymer && nMaxGap(g2)>0)
                        Q(i,jj) = s(n,e1)*sigma(e1,e2);
                        if (n==TestSiteNumber)
                            dQdE(e1,i,jj) = Q(i,jj);
                        end
                    end
                end % j
            end % g2
        end
    end
    
    % 15) g1-L-g2 gap followed by unit h2+1 of g2 protein:
    for e1=1:eMax
        for e2=1:eMax
            for g2=1:f
                for h2=0:m(g2)-1
                    if (nMaxGap(g2)>0)
                        i=R*(e1-1)+iGap(1,g2+1);
                        j=R*(e2-1)+mStart(g2)+h2;
                        if (n>1 && n<=Lpolymer-(m(g2)-h2))
                            Q(i,j) = s(n,e1)*sigma(e1,e2)*Unwrap(h2+1,g2)*c0(g2);
                            if (n==TestSiteNumber)
                                dQdE(e1,i,j)= Q(i,j);
                            end
                        end %n
                    end %(nMaxGap(g2)>0)
                end %g2
            end %g
        end
    end
    
end

end

