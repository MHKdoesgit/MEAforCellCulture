

function  p = generateElectricPulses(dur, amp, gap,pulsetype)


addval = 1;
switch lower(pulsetype)
    
    
    
    case {'cathodic-anodic','cath-an','ca','pos-neg','1',1}
        
        if length(amp)==1,    amp = [amp,amp];          end
        if length(gap)==1,    gap = [0, gap, 0];        end
        if length(dur)==1,    dur = [dur, dur];         end
        if amp(1)/abs(amp(1))==-1,   amp(1) = -amp(1);  end
        if amp(2)/abs(amp(2))== 1,    amp(2) = -amp(2); end
        p = [[0,amp(1),amp(1),0,0,amp(2),amp(2),0,0]; gap(1)+[0,addval,dur(1),dur(1)+addval,...
            dur(1)+gap(2)+addval+1,dur(1)+gap(2)+addval+2,dur(1)+gap(2)+dur(2),dur(1)+gap(2)+dur(2)+addval,...
            dur(1)+gap(2)+dur(2)+addval+gap(3)+1]];
        
    case {'anodic-cathodic','an-cath','ac','neg-pos','2',2}
        
        if length(amp)==1,    amp = [amp,amp];          end
        if length(gap)==1,    gap = [0, gap, 0];        end
        if length(dur)==1,    dur = [dur, dur];         end
        if amp(1)/abs(amp(1))== 1,   amp(1) = -amp(1);  end
        if amp(2)/abs(amp(2))==-1,    amp(2) = -amp(2); end
        
        p = [[0,amp(1),amp(1),0,0,amp(2),amp(2),0,0]; gap(1)+[0,addval,dur(1),dur(1)+addval,...
            dur(1)+gap(2)+addval+1,dur(1)+gap(2)+addval+2,dur(1)+gap(2)+dur(2),dur(1)+gap(2)+dur(2)+addval,...
            dur(1)+gap(2)+dur(2)+addval+gap(3)+1]];
        
        
    case {'anodic','an','a','neg','3',3}
        
        if length(gap)==1,          gap = [0, gap];         end
        if amp(1)/abs(amp(1))== 1,  amp(1) = -amp(1);       end
        p = [[0,amp(1),amp(1),0,0]; gap(1)+[0,addval,dur(1),dur(1)+addval,...
            dur(1)+gap(2)+addval+1]];
        
        
    case {'cathodic','cath','c','pos','4',4}
        
        if length(gap)==1,          gap = [0, gap];         end
        if amp(1)/abs(amp(1))== -1, amp(1) = -amp(1);       end
        p = [[0,amp(1),amp(1),0,0]; gap(1)+[0,addval,dur(1),dur(1)+addval,...
            dur(1)+gap(2)+addval+1]];
        
    case {'triphasic','pos-neg-pos','t','pnp','5',5}
        
        if length(amp)==1,    amp = [amp,amp,amp/2];        end
        if length(gap)==1,    gap = [0, gap,gap, 0];        end
        if length(dur)==1,    dur = [dur, dur,dur];         end        
        if amp(1)/abs(amp(1))==-1,      amp(1) = -amp(1);   end
        if amp(2)/abs(amp(2))== 1,      amp(2) = -amp(2);   end
        if amp(1)/abs(amp(3))==-1,      amp(3) = -amp(3);   end
        
        p = [[0,amp(1),amp(1),0,0,amp(2),amp(2),0,0,amp(3),amp(3),0,0]; ...
            gap(1)+[0,addval,dur(1),dur(1)+addval,dur(1)+gap(2)+addval+1,...
            dur(1)+gap(2)+addval+2,dur(1)+gap(2)+dur(2)+1,dur(1)+gap(2)+dur(2)+addval+1,...
            dur(1)+gap(2)+dur(2)+addval+gap(3)+2,dur(1)+gap(2)+dur(2)+addval+3+gap(3),...
            dur(1)+gap(2)+dur(2)+gap(3)+dur(3),dur(1)+gap(2)+dur(2)+gap(3)+dur(3)+addval,...
            dur(1)+gap(2)+dur(2)+gap(3)+dur(3)+gap(4)+2]];             
end
p = flipud(p);

end