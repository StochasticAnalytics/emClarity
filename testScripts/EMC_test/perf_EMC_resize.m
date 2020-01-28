function perf_EMC_resize

inSizes = {[100,100], [1000,1000], [4000,4000], [8000,8000], [50,50,50], [200,200,200]};
precision = {'single', 'double'};
limits = {[0,0,0,0], [30,30,30,30], [-30,-30,-30,-30], [30,-30,-30,30]};
options = {{'taper', false}; ...
           {'taper', true}; ...
           {'taper', false; 'origin', -1}; ...
           {'taper', true;  'origin', -1}};

for iPre = 1:length(precision)
  	for iSize = 1:length(inSizes)
        for iLim = 1:length(limits)
            for iOp = 1:length(options)
                
                b = limits{iLim};
                c = options{iOp};
                
                % is 3d
                if numel(inSizes{iSize}) == 3
                    extension = '3d';
                    a = ones(inSizes{iSize}(1,1), inSizes{iSize}(1,2), inSizes{iSize}(1,3));
                else
                    extension = '2d';
                    a =	ones(inSizes{iSize}(1,1), inSizes{iSize}(1,2));
                end
                
                % is fft
 
                b = limits{iLim};
                c = options{iOp};

                if strcmp(precision{iPre}, 'single')
                    a = single(a);
                end

                % perf cpu
                out = EMC_resize(a,b,c);
                ctime = timeit(@()EMC_resize(a,b,c));
                
                % perf gpu
                ga = gpuArray(a);
                gout = EMC_resize(ga,b,c);
                gtime = gputimeit(@()EMC_resize(ga,b,c));

                gerr = max(max(abs(gather(gout)-out)));

                fprintf('%dx%d\t%s \n', inSizes{iSize}(1,1), inSizes{iSize}(1,2), precision{iPre});
                disp(['Execution time on CPU = ',num2str(ctime)]);
                disp(['Execution time on GPU = ',num2str(gtime)]);
                fprintf('Maximum absolute error = %s \n', num2str(gerr));

                clear a out ga gout
            end
        end
    end
end

end