function [] = recScript(wanted_op)

modBin=10;    % this could be changed, but works well in most cases and is ill advised.
modThick=300; % this could be changed, particularly if your sample is thicker
              % than 3000 pixels. But then it is probably a bad candidate for
              % high resolution tomography anyhow ( assuming your pixel size > 1Ang) 

if strcmpi(wanted_op, 'build')
    % make a bin10 directory and create bin 10 tomos for evaluation.
    mkdir('bin10');

    if ~isfolder('aliStacks')
        error('No aliStacks folder found. Please run the alignment script first.');
    end

    stack_names = dir('aliStacks/*.fixed');
    stack_names = {stack_names.name};
    nParProcesses = 4;
    try
        EMC_parpool(nParProcesses+1)
    catch
        delete(gcp('nocreate'))
        EMC_parpool(nParProcesses+1)
    end

    parfor iStack = 1:length(stack_names)
        fprintf('Working on stack %s\n', stack_names{iStack});
        bin_cmd = sprintf('newstack -bin %d aliStacks/%s bin10/%s_bin10.fixed', modBin, stack_names{iStack}, stack_names{iStack}(1:end-11));
                           rec_cmd = sprintf('tilt -input bin10/%s_bin10.fixed -output bin10/%s_bin10.rec -TILTFILE fixedStacks/%s.tlt -RADIAL 0.15,0.05 -UseGPU 0 -THICKNESS %d -RotateBy90',  ...
         stack_names{iStack}(1:end-11), stack_names{iStack}(1:end-11), stack_names{iStack}(1:end-11), modThick);
        errmsg = system(sprintf('%s > /dev/null', bin_cmd));
        if errmsg
            fprintf('Error binning the stack');
            system(sprintf('%s', bin_cmd));
            error('Error binning the stack');
        end
        errmsg = system(sprintf('%s > /dev/null', rec_cmd));
        if errmsg
            fprintf('Error reconstructing the stack');
            system(sprintf('%s', rec_cmd));
            error('Error reconstructing the stack');
        end
    end
elseif strcmpi(wanted_op, 'recon')

    mkdir('recon');

    if ~isfolder('bin10')
        error('No bin10 folder found. Please run recScript(build) first.');
    end

    mod_names = dir('bin10/*.mod');
    mod_names = {mod_names.name};
    if isempty(mod_names)
        error('No mod files found in the bin10 folder');
    end
    
    for iMod = 1:length(mod_names)
        mod_basename = mod_names{iMod}(1:end-10);
        if ~isfile(sprintf('aliStacks/%s_ali1.fixed', mod_basename))
            error('Did not find the aligned, fixed stack at %s_ali1.fixed', mod_basename);
        end
        header = getHeader(MRCImage(sprintf('aliStacks/%s_ali1.fixed', mod_basename), 0));
        
        if ~isfolder('fixedStacks')
            error('No fixedStacks folder found. Please run the alignment script first.');
        end

        if isfile(sprintf('fixedStacks/%s.local', mod_basename))
            local_file = sprintf('-LOCALFILE fixedStacks/%s.local', mod_basename);
        else
            local_file = '';
        end

        cmd = sprintf('model2point -contour bin10/%s bin10/mod.tmp', mod_names{iMod});
        msg = system(sprintf('%s > /dev/null', cmd));
        if msg
            fprintf('Error running model2point');
            system(sprintf('%s', cmd));
            error('Error running model2point');
        end
        mod_f = importdata('bin10/mod.tmp');
        mod_f = sortrows(mod_f, 1);
        if mod(size(mod_f,1), 6) ~= 0
            error('The number of contours is not a multiple of 6 %s', mod_names{iMod});
        end
        rmpath('bin10/mod.tmp');
        numParticles = size(mod_f, 1)/6;
        fid = fopen(sprintf('recon/%s_recon.coords', mod_basename), 'w');
        fprintf(fid, '%s\n', mod_basename);
        fprintf(fid, '%d\n', numParticles);

        for iRegion = 1:numParticles
            xLt = mod_f((iRegion-1)*6+1, 2)*modBin;
            xHt = mod_f((iRegion-1)*6+2, 2)*modBin;
            yLt = mod_f((iRegion-1)*6+3, 3)*modBin;
            yHt = mod_f((iRegion-1)*6+4, 3)*modBin;
            zLt = mod_f((iRegion-1)*6+5, 4)*modBin;
            zHt = mod_f((iRegion-1)*6+6, 4)*modBin;

            xL = min(xLt, xHt);
            xH = max(xLt, xHt);
            yL = min(yLt, yHt);
            yH = max(yLt, yHt);
            zL = min(zLt, zHt);
            zH = max(zLt, zHt);

            width = int32(xH-xL);
            height = int32(yH-yL);
            yLow = int32(yL);
            yTop = int32(yLow + height - 1);
            oX = int32((xL + (width/2.0)) - header.nX/2.0);
            thickness = int32(sqrt((zH - zL)^2));
            oZ = int32((zL+thickness/2.0) - (modThick*modBin)/2.0);
            oY = int32(yL + height/2.0 - header.nY/2.0);

            fprintf(fid, '%d\n%d\n%d\n%d\n%d\n%d\n', width, height, thickness, oX, oY, oZ);

        end
    end

else
    error('Invalid operation. Please use "build" or "recon"');
end
    

    


end