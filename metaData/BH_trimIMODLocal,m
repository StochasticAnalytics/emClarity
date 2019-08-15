function [ ] = BH_trimIMODLocal(fileName,nTilts,skip)

system(sprintf('cp %s %s_orig',fileName,fileName));

fileID = fopen(sprintf('%s_orig',fileName),'r');

p = textscan(fileID, '%s', 'CommentStyle',{'%'},'Delimiter','\n','TreatAsEmpty',{' '});

fOUT = fopen(fileName,'w');

nLines = length(p{1});

readingTilts = true;
nRead = 0;

fprintf(fOUT,'    %s\n',p{1}{1});
nAnglesPerLine = length(strsplit(p{1}{2}));
nPrinted = 0;

for iLine = 2:nLines
  
  thisLine = strsplit(p{1}{iLine});
  nFields = length(thisLine);


   if (readingTilts)
     % fprintf('\n');remove one record at a time
     for iField = 1:nFields
       nRead = nRead + 1;
       if ~ismember(nRead,skip)
        fprintf(fOUT,' % 5.2f ',str2double(thisLine{iField}));
        newLinePrinted = false;
        nPrinted = nPrinted + 1;
        if nPrinted == nAnglesPerLine
          newLinePrinted = true;
          fprintf(fOUT,'\n');
          nPrinted = 0;
        end
       else

       end   
       
      if nRead == nTilts
        needNewLine = false;
        readingTilts = false;
        nRead = 0;
        nPrinted = 0;
         break
      end
     end
     
     if (readingTilts && newLinePrinted)
       fprintf(fOUT,'\n');
     end
   else
     if ~(newLinePrinted)
       fprintf(fOUT,'\n');
       newLinePrinted = true;
     end
     nRead = nRead + 1;
     if ~ismember(nRead,skip)
 
      fprintf(fOUT,'  %s\n',p{1}{iLine});
   
     end
     if nRead == nTilts
       readingTilts = true;
       newLinePrinted = false;
       nRead = 0;
     end    
 
   
  end

  
end
fclose(fileID);
fclose(fOUT);

end
