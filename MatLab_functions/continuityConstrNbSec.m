function [isNbFeasible]=continuityConstrNbSec(nbcc,nb3)

for i=1:3
    if nbcc(i)==0
        isNbFeasible=true;
    elseif all([nbcc(i)>0,nbcc(i)==2])
        isNbFeasible=true;
	elseif all([nbcc(i)>2,mod(nbcc(i),2)==0])
		if nb3(i)==nbcc(i)
			isNbFeasible=true;
		elseif nb3(i)==nbcc(i)+nbcc(i)-1
			isNbFeasible=true;
		else
			isNbFeasible=false;
		end
    elseif all([nbcc(i)>0,mod(nbcc(i),2)==1]) % number of rebar in the first layer is odd
        if mod(nb3(i),2)==0
            isNbFeasible=false;
        else
            isNbFeasible=true;
        end
    end
end

