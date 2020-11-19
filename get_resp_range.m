function range = get_resp_range(P)
    
    sum = 0.;
    avg = 0.;
    cnt = 0;
    r_start = 0;
    r_end = 0;
    
    for ii = 1:length(P)
%         P(ii)
        if P(ii) == 0
            continue;
        else
           if ii < 50
               cnt = cnt + 1;
               sum = sum + P(ii);
               avg = sum/cnt;
           else
               if P(ii) > avg*3
                   if r_start == 0
                       r_start = ii;
                   else
                       r_end = ii;
                   end
               else 
                   cnt = cnt + 1;
                   sum = sum + P(ii);
                   avg = sum/cnt;
               end
               if r_start ~= 0 && ii-r_start > 50 
                   break;
               end
           end
        end
        
    end
    
    range = r_start:r_end;

end