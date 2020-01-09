function [jump_Re, jump_Im, I_Re, I_Im] = check_furth(filename)
    data = load(filename);
    r = data(:,1);
    k_z = data(:,2);
    k_theta = data(:,3);
    Bpmn_r_Re = data(:,4);
    Bpmn_r_Im = data(:,5);
    furth_2_Re = data(:,6); I_Re = furth_2_Re(furth_2_Re ~= 0.0);
    furth_2_Im = data(:,7); I_Im = furth_2_Im(furth_2_Im ~= 0.0);
    
    k2 = k_z .* k_z + k_theta .* k_theta;
    furth_1_Re = r ./ k2 .* gradient(Bpmn_r_Re, r);
    furth_1_Im = r ./ k2 .* gradient(Bpmn_r_Im, r);
    
    %{
    imax_Re = find(islocalmax(furth_1_Re));
    imax_Re = imax_Re(1);
    fmax_Re = furth_1_Re(imax_Re);
    imin_Re = find(islocalmin(furth_1_Re));
    imin_Re = imin_Re(1);
    fmin_Re = furth_1_Re(imin_Re);
    imax_Im = find(islocalmax(furth_1_Im));
    imax_Im = imax_Im(1);
    fmax_Im = furth_1_Im(imax_Im);
    imin_Im = find(islocalmin(furth_1_Im));
    imin_Im = imin_Im(1);
    fmin_Im = furth_1_Im(imin_Im);
    %}
    [fmax_Re, imax_Re] = max(furth_1_Re);
    [fmin_Re, imin_Re] = min(furth_1_Re);
    [fmax_Im, imax_Im] = max(furth_1_Im);
    [fmin_Im, imin_Im] = min(furth_1_Im);
    jump_Re = (fmax_Re - fmin_Re) * sign(imax_Re - imin_Re);
    jump_Im = (fmax_Im - fmin_Im) * sign(imax_Im - imin_Im);
    
    figure;
    plot(r, furth_1_Re, '-k', r, furth_1_Im, '--r', ...
        r(imax_Re), fmax_Re, '*k', r(imin_Re), fmin_Re, '*k', ...
        r(imax_Im), fmax_Im, '*r', r(imin_Im), fmin_Im, '*r');
    xlabel('r / cm');
    ylabel('integrand 1 / G cm^{-2}');
    legend({'real', 'imag'}, 'Location', 'northwest');
end
