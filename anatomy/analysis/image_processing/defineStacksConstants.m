%==========================================================================
% define some useful constants for dealing with confocal stacks
%==========================================================================

% channel indices of confocal stacks
BLUE = 1 ; 
BFIELD = 1 ; 

RED = 2 ; 
SYNAPTO = 2 ; 

GREEN = 3 ; 
GRP = 3 ; 

% image type conversions
uint16_to_double = (1 /  65535) ; 
uint8_to_double = (1 / 255) ; 

uint16_to_single = single(1 /  65535) ; 
uint8_to_single = single(1 / 255) ; 