% input orders to calculate and specify positions pos = [x_pos y_pos z_pos]
% where x_pos, y_pos, and z_pos are each column vectors of spatial positions
%%%use vector for orders to calculate
function [harm_all] = calc_spherical_harmonics_arb_points(orders_to_calculate,pos)



N = numel(pos(:,1));


ii=0;

for oo = 1:numel(orders_to_calculate); disp(num2str(oo))
    
    n = orders_to_calculate(oo);
    m = -orders_to_calculate(oo):1:orders_to_calculate(oo);
    num_harm = numel(m);


    for mm = 1:num_harm
        ii = ii+1;
        for pp = 1:N
          

                    temp(pp) = leg_rec_harmonic(n,m(mm),pos(pp,1),pos(pp,2),pos(pp,3));
                    

          
        end
        
        harm_all(:,ii) = temp;

    end




end


