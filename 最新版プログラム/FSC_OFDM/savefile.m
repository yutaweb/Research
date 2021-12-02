function savefile(filename, m, headerstr, dlm)
%   SAVEFILE Write ASCII delimited file for a matrix with a header string.
%   SAVEFILE(FILENAME,M,HEADERSTR,DLM) writes matrix M into FILENAME using the
%   character DLM as the delimiter.  Specify '\t' to produce 
%   tab-delimited files.
%   Default delimiter is a comma.
%   The headerstr is added to the first line of the file for indicating the columns
%   SAVEFILE(FILENAME,M) will save M with no header using comma separation

%   Modified from DLMWRITE script by Chang-Jun Ahn 19/8/2000
%   Brian M. Bourgault 10/22/93
%   Copyright (c) 1984-98 by The MathWorks, Inc.
%   $Revision: 5.14 $  $Date: 1997/11/21 23:35:06 $
%
%   test for proper filename

if ~isstr(filename),
    error('FILENAME must be a string.');
end;

if nargin < 3, 
   error('Requires at least 2 input arguments.'); 
end

NEWLINE = sprintf('\n');
if strncmp(computer,'PCWIN',5)
   NEWLINE = sprintf('\r\n');
end

% delimiter defaults to Comma for CSV
if nargin < 4, dlm = ','; end
dlm = sprintf(dlm); % Handles special characters.


if nargin < 3, headerstr = ''; end
   
% open the file
if strncmp(computer,'MAC',3)
  fid = fopen(filename ,'wt');
else
  fid = fopen(filename ,'wb');
end

if fid == (-1), 
   error(['Could not open file ' filename]); 
end

fwrite(fid,headerstr,'char');
fwrite(fid, NEWLINE, 'char');
% dimensions size of matrix
[br,bc] = size(m);


% start dumping the array, for now number format float
for i = 1:br
    for j = 1:bc
            str = num2str(m(i,j));
            fwrite(fid, str, 'uchar');    
        if(j < bc)
            fwrite(fid, dlm, 'uchar');    
        end
    end
    fwrite(fid, NEWLINE, 'char'); % this may \r\n for DOS 
end

% close files
fclose(fid);
