function sm_gen_swi(reg, row, res)
    USE_CORRECTION = 0;
    NUM_YEARS = 6; % 2009-2014
    NUM_DAYS = 365;
    
    
    if strcmp('sir',res)
        [lmask, ~] = loadsir(['landmasks/' reg '.sir.lmask']);
        [~, num_columns] = size(lmask);  
    elseif strcmp('grd',res)
        [lmask, ~] = loadsir(['landmasks/' reg '.grd.lmask']);
        [~, num_columns] = size(lmask); 
    else
        error('Resolution option not supported!');
    end
    
    [arid_mask, ~, ~, ~] = loadsir(['climate/' reg '.' res]);
    arid_mask = flipud(arid_mask == 4 | arid_mask == 5 | arid_mask == 6 | arid_mask == 7);
    
    num_rows = 1;
    % Declare storage arrays
    SWI_stor = nan(NUM_DAYS,NUM_YEARS,num_columns,1);
    ms_stor = nan(NUM_DAYS,NUM_YEARS,num_columns,1);
    sig0_dry_stor = nan(NUM_DAYS,NUM_YEARS,num_columns,1);
    
    % Open NETCDF Files
    % first get the c0 values
    basedir = '/auto/temp/lindell/soilmoisture/';
    
    c0_file = [basedir 'c0/c0_' reg '.nc'];    
    disp('Loading data from NetCDF Files');

    c0_id = netcdf.open(c0_file,'NC_NOWRITE');
    c0_dry_varid = netcdf.inqVarID(c0_id,'dry');
    c0_dry_data = netcdf.getVar(c0_id,c0_dry_varid)';

    c0_wet_varid = netcdf.inqVarID(c0_id,'wet');
    c0_wet_data = netcdf.getVar(c0_id,c0_wet_varid)';

    if USE_CORRECTION,
        load(['/home/lindell/research/soil_moisture/c0_wet/' reg '.mat']); % Load the corrected values for c0 wet (use sm_view_c0 to gen.)
        c0_wet_data = c0_wet_corrected;
    end
    
    netcdf.close(c0_id);
    
    % Now get the sigma0 time series
    ts_a = [basedir 'ts/ts_' reg '_a.nc'];
    ts_b = [basedir 'ts/ts_' reg '_b.nc'];
    
    ncdisp(ts_a);
    ncdisp(ts_b);
    disp('Loading data from NetCDF Files');
    
    ts_a_id = netcdf.open(ts_a,'NC_NOWRITE');
    ts_b_id = netcdf.open(ts_b,'NC_NOWRITE');
    
    ts_a_varid = netcdf.inqVarID(ts_a_id,'data');
    ts_b_varid = netcdf.inqVarID(ts_b_id,'data');
    
    ts_a_data = netcdf.getVar(ts_a_id,ts_a_varid,[0,0,0,row-1],[365,6,num_columns,1]);
    ts_b_data = netcdf.getVar(ts_b_id,ts_b_varid,[0,0,0,row-1],[365,6,num_columns,1]);
    
    netcdf.close(ts_a_id);
    netcdf.close(ts_b_id);

    disp('Done loading data, beginning processing.');
    
    % adjust datasets
    disp('Setting no-data values to NaN');
    c0_dry_data(c0_dry_data == 0) = nan;
    c0_wet_data(c0_wet_data == 0) = nan;
    ts_a_data(ts_a_data == 0) = nan;
    ts_a_data(abs(abs(ts_a_data)-33) < .01) = nan;
    ts_b_data(ts_b_data == 0) = nan;
    ts_b_data(abs(abs(ts_b_data)-3) < .01) = nan;
    % Get rid of outliers in the slope values
    iq = iqr(ts_b_data(:));
    av = mean(ts_b_data(:),'omitnan');
    ts_b_data(ts_b_data > av + 2*iq) = nan;
    ts_b_data(ts_b_data < av - 2*iq) = nan;


    % Scipal Disseratation p. 49
    theta_wet = 40;
    theta_dry = 25;
    
    disp('Beginning iterations');
    for col = 1:num_columns,
        fprintf('Row: %04d Col: %04d\n',row,col);

        % organize data into single time series
        ts_a_iter = squeeze(ts_a_data(:,:,col,1));
        ts_a_iter = double(reshape(ts_a_iter,1,NUM_YEARS*NUM_DAYS));

        ts_b_iter = squeeze(ts_b_data(:,:,col,1));
        fit_data = double(mean(reshape(ts_b_iter',NUM_YEARS,NUM_DAYS)));
%         fit_data = double(reshape(ts_b_iter,1,NUM_YEARS*NUM_DAYS));
        
        nan_mask = ~isnan(fit_data);
        fit_days = 1:length(fit_data);
        fit_days = fit_days(nan_mask)';

        % Check if we have data from this time series
        if length(fit_days) < 15,
            continue;
        end
        
        sinfit = fit(fit_days,fit_data(nan_mask)','fourier3');
        sinfit = repmat(sinfit(1:365),6,1);
        sinfit = fit((1:length(sinfit)).',sinfit,'fourier2');
        
        % Get c0_dry, c0_wet
        c0_dry = c0_dry_data(row,col);
        c0_wet = c0_wet_data(row,col);

        if (isnan(c0_wet) || isnan(c0_dry))
            continue;
        end

        sigma0_dry = c0_dry - (sinfit(1:length(ts_a_iter)))*(theta_dry - 40);
        sigma0_wet = c0_wet - (sinfit(1:length(ts_a_iter)))*(theta_wet - 40);
        
        
        if (arid_mask(row,col) == 1)
            sigma0_sens = sigma0_wet - sigma0_dry;
            meas_mask = sigma0_sens < 5;
            sigma0_wet(meas_mask) = sigma0_dry(meas_mask) + 5;
        end
        
%         sigma0_dry = repmat(c0_dry,[1,2190])';
        
        % for the time series file, get the moisture for each point
        ms = (ts_a_iter' - sigma0_dry) ./ (sigma0_wet - sigma0_dry);

        % calculate soil water index
        T = 20; % scipal p. 75
        SWI = nan(length(ts_a_iter),1);

        data_days = find(~isnan(ms));
        for ii = 1:length(data_days),
            SWI(data_days(ii)) = sum(ms(data_days(1:ii)).* ... 
                exp(-(data_days(ii)-data_days(1:ii))/T)) ./ ... 
                sum(exp(-(data_days(ii)-data_days(1:ii))/T));
        end
        
        % Save in buffer
        SWI_stor(:,:,col,1) = reshape(SWI,NUM_DAYS,NUM_YEARS); 
        ms_stor(:,:,col,1) = reshape(ms,NUM_DAYS,NUM_YEARS); 
        sig0_dry_stor(:,:,col,1) = reshape(sigma0_dry,NUM_DAYS,NUM_YEARS);
    end
    
    
    disp('Saving NetCDF Files...');
    mode = netcdf.getConstant('NETCDF4');
    mode = bitor(mode,netcdf.getConstant('CLOBBER'));
    swi_ncid = netcdf.create([basedir 'swi/swi_' reg '_' sprintf('%04d',row) '.nc'],mode);

    dimids(4) = netcdf.defDim(swi_ncid,'row',num_rows);
    dimids(3) = netcdf.defDim(swi_ncid,'column',num_columns);
    dimids(2) = netcdf.defDim(swi_ncid,'year',NUM_YEARS);
    dimids(1) = netcdf.defDim(swi_ncid,'day',NUM_DAYS);
    
    swi_varid = netcdf.defVar(swi_ncid,'swi','NC_FLOAT',dimids);
    ms_varid = netcdf.defVar(swi_ncid,'ms','NC_FLOAT',dimids);
    sig0_dry_varid = netcdf.defVar(swi_ncid,'dry','NC_FLOAT',dimids);
    
    netcdf.endDef(swi_ncid);
    
    netcdf.putVar(swi_ncid,swi_varid,SWI_stor);
    netcdf.putVar(swi_ncid,ms_varid,ms_stor);
    netcdf.putVar(swi_ncid,sig0_dry_varid,sig0_dry_stor);
    netcdf.close(swi_ncid);
   
    disp('Done.');

    
end
