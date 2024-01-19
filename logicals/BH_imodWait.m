function [  ] = BH_imodWait( fileName )
%Wait on a file that is slowly written to disk
%   system should not return until the process is finished, but in some
%   instances, IMOD functions return before this is true.
%   1) in observed cases the header is not written so use this as a check
%   for completion
%   2) Make sure the file is actually growing or throw an error.
fSize = 0;
failed = 1;
while (failed)
  [failed, ~] = system(sprintf('header  %s',fileName));
  fInfo = dir(fileName);
  if fInfo.bytes > fSize;
    fSize = fInfo.bytes
    pause(5);
  else
    % give it another longer pause
    fprintf('System seems slow, extended pause.\n')
    pause(30)
    [failedMeTwice, ~] = system(sprintf('header  %s',fileName));
    if (failedMeTwice)
      if fInfo.bytes > fSize
        fSize = fInfo.bytes
        pause(5)
      else
        error('##### system impatience.')
      end
    end
  end
end

end

