function [nb1,nb2,nb3,isfeasible]=nb3DSec1(nbmax,Aos,ab,nbcc)
    
    Abmax1 = nbmax(1) * ab(1) ;
    Abmax2 = nbmax(2) * ab(2) ;
    Abmax3 = nbmax(3) * ab(3) ;
    
    if sum(nbcc.*ab)>=Aos
        isfeasible=true;
        nb1=nbcc(1);
        nb2=nbcc(2);
        nb3=nbcc(3);
    else
        if Abmax1 >= Aos % only one rebar layer is necessary
 
            isfeasible=true;
            
            if any([nbcc(1)==2,nbcc(1)==0]) 
                nb1 = min([ceil(Aos / ab(1)),nbmax(1)]) ;
                if nb1<2
                    nb1=2;
                end
            elseif nbcc(1)>2
                nb1 = nbcc(1) ;
            end
            
            Ab1=nb1*ab(1);
            if Aos-Ab1>0
                nb2=ceil((Aos-Ab1) / ab(2));
                if nb2>nb1
                    nb2=nb1;
                end
                Ab2=nb2*ab(2);
                if Aos-Ab1-Ab2 > 0
                    nb3=ceil((Aos-Ab1-Ab2) / ab(3));
                else
                    nb3=0;
                end
                if nb3>nb2
                    nb3=nb2;
                end
            else
                nb2=max([0,nbcc(2)]);
                nb3=max([0,nbcc(3)]);
            end

        elseif Abmax1 + Abmax2 >= Aos % only two rebar layers are necessary

            isfeasible=true;
            
            if any([nbcc(1)==2,nbcc(1)==0]) 
                nb1= min(ceil(Aos / ab(1)),nbmax(1)) ;
                if nb1<2
                    nb1=2;
                end
            elseif nbcc(1)>2
                nb1=nbcc(1);
            end
            
            Ab1=nb1*ab(1);
            if any([nbcc(2)==0,nbcc(2)==2])
                nb2=min([ceil((Aos-Ab1)/ab(2)),nbmax(2)]);
                if nb2>nb1
                    nb2=nb1;
                end
                Ab2=nb2*ab(2);
                if Aos - Ab1 - Ab2 > 0
                    nb3=min([ceil((Aos-Ab1-Ab2)/ab(3)),nbmax(3)]);
                else
                    nb3=nbcc(3);
                end
            elseif nbcc(2)==1
                nb2=min([ceil((Aos-Ab1)/ab(2)),nbmax(2)]);
                if mod(nb2,2)==0
                    if any([nb2==nbmax(2),nb2==nb1]) 
                        nb2=nb2-1;
                        if nbcc(3)==0
                            nb3=1;
                        else
                            nb3=nbcc(3);
                        end
                    else
                        nb2=nb2+1;
                        nb3=nbcc(3);
                    end
                else
                    nb3=nbcc(3);
                end
            elseif nbcc(2)>2
                nb2=nbcc(2);
                Ab2=nb2*ab(2);
                if Aos - Ab1 - Ab2 > 0
                    nb3=min([ceil((Aos-Ab1-Ab2)/ab(3)),nbmax(3)]);
                    if nb3>nb2
                        nb3=nb2;
                    end
                else
                    nb3=nbcc(3);
                end
            end

            if mod(nb1,2)==0 && mod(nb2,2)~=0
                if any([nb2==nbmax(2),nb2==nb1])
                    nb2=nb2-1;
                    nb3=nb3+1;
                else
                    nb2=nb2+1;
                end
            end
            if mod(nb2,2)==0 && mod(nb3,2)~=0
                nb3=nb3+1;
            end

        elseif Abmax1 + Abmax2 + Abmax3 >= Aos

            isfeasible=true;
            if any([nbcc(1)==2,nbcc(1)==0]) 
                nb1= min([ceil(Aos / ab(1)),nbmax(1)]) ;
                if nb1<2
                    nb1=2;
                end
            elseif nbcc(1)>2
                nb1=nbcc(1);
            end
            
            Ab1=nb1*ab(1);
            if any([nbcc(2)==0,nbcc(2)==2])
                nb2=min(ceil((Aos-Ab1)/ab(2)),nbmax(2));
                if nb2>nb1
                    nb2=nb1;
                end
            elseif nbcc(2)==1
                nb2=min(ceil((Aos-Ab1)/ab(2)),nbmax(2));
                if nb2>nb1
                    nb2=nb1;
                end
                if mod(nb2,2)==0
                    if any([nb2==nbmax(2),nb2==nb1])
                        nb2=nb2-1;
                        if nbcc(3)==0
                            nb3=1;
                        else
                            nb3=nbcc(3);
                        end
                    else
                        nb2=nb2+1;
                        nb3=nbcc(3);
                    end
                else
                    nb3=nbcc(3);
                end
            elseif nbcc(2)>2
                nb2=nbcc(2);
                nb3=nbcc(3);
            end
            
            Ab2=nb2*ab(2);
            if any([nbcc(3)==0,nbcc(3)==2])
                nb3 = min([ceil((Aos-Ab1-Ab2)/ab(3)),nbmax(3)]);
                if nb3>nb2
                    nb3=nb2;
                end
            elseif nbcc(3)==1
                nb3=ceil(min((Aos-Ab1-Ab2)/ab(3),nbmax(3)));
                if mod(nb3,2)==0
                    if nb3==nbmax(2)
                        nb3=nb3-1;
                    else
                        nb3=nb3+1;
                    end
                end
            elseif nbcc(3)>2
                nb3=nbcc(3);
            end
    
            if mod(nb1,2)==0 && mod(nb2,2)~=0
                nb2=nb2+1;
            end
    
            if mod(nb2,2)==0 && mod(nb3,2)~=0
                nb3=nb3+1;
            end
        else
            nb1 = 0;
            nb2 = 0;
            nb3 = 0;
            isfeasible=false;
        end
    end
end
