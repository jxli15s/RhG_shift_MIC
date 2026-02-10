 function write_hr(g,filename)
    %open a file to write wannier hoppings
    fileID=fopen(filename,'w');
    
    [numBands,~,numRpts]=size(g.ham);

    % Get the current date and time
    currentTime = datetime('now');
    formattedTime = datestr(currentTime, 'mm/dd/yyyy at HH:MM:SS');

    % Write the header with the current data and time to the file
    fprintf(fileID,' write on %s\n', formattedTime);
    % Write the orbital numbers to the file
    fprintf(fileID,'\t %d\n', numBands);
    % Write the sites numbers to the file
    fprintf(fileID,'\t %d\n', numRpts);
    
    % numbers per line
    numsPerLine = 15;
    
    % get the site numbers of hoppings
    len=size(g.hopr,1);

    degeneracy = ones(1,len);

    % Write the data into the file with 15 numbers per line
    for i = 1:numsPerLine:len
        % Determine the index range of data for the current line
        endIdx = min(i+numsPerLine-1,len);
        fprintf(fileID, '%5d', degeneracy(i:endIdx));
        fprintf(fileID, '\n');
    end

    for i = 1:numRpts
        for j = 1:numBands
            for k = 1:numBands
                fprintf(fileID, '%5d %5d %5d %5d %5d %12.6f %12.6f\n',g.hopr(i,:),j,k,real(g.ham(j,k,i)),imag(g.ham(j,k,i)));
            end
        end
    end

    
 end
