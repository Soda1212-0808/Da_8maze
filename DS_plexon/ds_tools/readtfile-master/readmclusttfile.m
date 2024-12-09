function [ timestamp, numSpikes, hdr ] = readmclusttfile( sFilePath )
%READMCLUSTTFILE   Reads a cluster t-file produced by MClust.
%   Inputs:
%     sFilePath -   This is the full path to the t-file.
%
%   Outputs:
%     timestamp:        A list of the timestamps of the spikes in increasing
%                       order in units of 10^-4 seconds.
%     numSpikes:        Number of spikes in the cluster.
%     hdr:                Header of the .t file.

%  23 April 2002, C. Higginson, created.

% Check integer type. The .t files might be saved as either 32 bit integers
% (with .t extension) or 64 bit integers (with .t64 extension).
[ p, file, ext ] = fileparts( sFilePath );

flag64 = regexp( ext, '64' );

if ~isempty( flag64 )
    intType = 'uint64';
    
else
    intType = 'uint32';
    
end

% Quick check to see if being called correctly.
if nargin ~= 1
    error( 'There should be exactly 1 input argument.' );
    
end;

% Open for binary read access
iClusterFileID = fopen( sFilePath, 'r', 'b' );

% Throw an error if there is a problem reading the file.
if ( iClusterFileID == -1 )
    error( [ 'Error opening cluster file:' sFilePath ] );
    
end;

% Find the end of the header.
fseek( iClusterFileID, 0, 'bof' );
beginheader = '%%BEGINHEADER';
endheader = '%%ENDHEADER';
iH = 1;
hdr = { };
curfpos = ftell( iClusterFileID );
headerLine = fgetl( iClusterFileID );
if strcmp( headerLine, beginheader )
    hdr{ 1 } = headerLine;
    
    while ~feof( iClusterFileID ) && ~strcmp( headerLine, endheader )
        headerLine = fgetl( iClusterFileID );
        iH = iH + 1;
        hdr{ iH } = headerLine;
        
        if strcmp( headerLine, endheader )
            break;
            
        end;
        
    end;
    
end;

% Read all of the timestamps.
[ timestamp, numSpikes ] = fread( iClusterFileID, inf, intType );
timestamp = timestamp';

% It is important that the timestamps are sequential for later analysis.
delta = diff( timestamp );
if min( delta ) < 0
    error( [...
        'Spike timestamps in video file must be non-decreasing. '...
        sFilePath ] );
    
end;

% Give a warning if the cell has no spikes.
if (numSpikes == 0)
    warning(['There are no spikes in the cluster file:' sFilePath]);
end;

% Free memory.
delta = [];

fclose(iClusterFileID);

return;
