%--------------------------------------------------------------------------
%% should make accessing kinefly data a little easier

NAMES = {'CAM3 - WBA L', 'CAM3 - WBA R', 'CAM3 - Head', 'CAM3 - Aux.',...
        'CAM2 - L Dev', 'CAM2 - L Rad', 'CAM2 - R Dev', 'CAM2 - Aux.',...
        'CAM1 - R Dev', 'CAM1 - R Rad', 'CAM1 - L Dev', 'CAM1 - Aux.',...
        'WBF', 'X position', 'Y position', 'Stimulus',...
        'CW Torque', 'Thrust', 'dX', 'Stripe Position'};

UNITS = {'rad', 'rad', 'rad', 'int.',...
    'rad', 'px', 'rad', 'int.',...
    'rad', 'px', 'rad', 'int.',...
    'Hz', 'pos', 'pos', 'V',...
    '', '', 'pix/s', 'pix'};

% array indices
L_AMP   = 1 ;
R_AMP   = 2 ; 
HEAD    = 3 ; 
AUX_3   = 4 ; 
L_DEV_2 = 5 ; 
L_RAD_2 = 6 ; 
R_DEV_2   = 7 ; 
AUX_2   = 8 ; 
R_DEV_1 = 9 ; 
R_RAD_1 = 10 ; 
L_DEV_1 = 11 ; 
AUX_1   = 12 ;
WBF     = 13 ; 
X_POS   = 14 ; 
Y_POS   = 15 ; 
STIM    = 16 ; 
CW_TRQ  = 17 ; 
THRUST  = 18 ; 
DX      = 19 ; 
STIPE   = 20 ; 

%other useful params
RAD2DEG     = 180/pi ; 
DEG2RAD     = pi/180 ; 
PULSE_WIDTH = 0.1 ; %seconds. note: this is assumed by code. should check

%--------------------------------------------------------------------------
%% filter and time subsampling preferences
N_BOOT_SAMPLES = 500 ; % # of boostrap samples to use, primaril for CI
N_SUB_SAMPLE = 200 ;    % subsampling time. 
MEDFILT_WINDOW = 3001 ; % window size for median filter

SAMPLE_RATE=20000 ; % frequncy at which data is acquired
T_START = 0 ; % in seconds
T_END = 6 ;   % in seconds

% since we're using spectrogram frequency calculation, there are now some
% extrapolation errors--need to deal with those
FREQ_RANGE = [120, 280] ;
FREQ_DIFF_THRESH = 0.01 ; % Hz per time step--even this is too large
%--------------------------------------------------------------------------
%% prefences for LED strip in plots
LED_color = 'r' ; 
LED_alpha = 0.4 ; 
LED_duration = 0.1 ; % seconds (100 ms)

LED_start_time = 1.0 ; % seconds

%--------------------------------------------------------------------------
%% general plot appearance pref
CI_alpha = 0.2 ; 

%--------------------------------------------------------------------------
%% color maps for flies
cmap_struct = struct() ; 

% motor lines
cmap_struct.SS37252 = 'Purples' ; 
cmap_struct.SS37253 = 'Purples' ; 
cmap_struct.SS01062 = 'Greys' ; 
cmap_struct.SS37246 = 'Blues' ; 
cmap_struct.SS37231 = 'Blues' ; 
cmap_struct.SS40864 = 'YlOrRd' ; 
cmap_struct.SS45772 = 'GnBu' ;
cmap_struct.SS45778 = 'PuBuGn' ;
cmap_struct.SS45782 = 'PuBu' ;
cmap_struct.SS47120 = 'YlGn' ;
cmap_struct.SS47152 = 'Greens' ;
cmap_struct.SS01055 = 'Greys' ; 
cmap_struct.SS48311 = 'BuPu' ; 
cmap_struct.SS47161 = 'BuPu' ; 
cmap_struct.SS41068 = 'Oranges' ; 
cmap_struct.SS44056 = 'Reds' ; 
cmap_struct.SS41039 = 'PuBu' ; 
%{
% OLD VERSION
cmap_struct.SS37252 = 'YlOrRd' ; 
cmap_struct.SS37253 = 'Reds' ; 
cmap_struct.SS01062 = 'Greys' ; 
cmap_struct.SS37246 = 'Blues' ; 
cmap_struct.SS37231 = 'Purples' ; 
cmap_struct.SS40864 = 'Oranges' ; 
cmap_struct.SS45772 = 'GnBu' ;
cmap_struct.SS45778 = 'PuBuGn' ;
cmap_struct.SS45782 = 'PuBu' ;
cmap_struct.SS47120 = 'RdPu' ;
cmap_struct.SS47152 = 'Greens' ;
cmap_struct.SS01055 = 'Greys' ; 
cmap_struct.SS48311 = 'Oranges' ; 
cmap_struct.SS47161 = 'Oranges' ; 
cmap_struct.SS41068 = 'Purples' ; 
cmap_struct.SS44056 = 'RdPu' ; 
cmap_struct.SS41039 = 'PuBu' ; 
%}
% sensory lines
cmap_struct.SS36076 = 'Reds' ;
cmap_struct.SS40775 = 'Blues' ;
cmap_struct.SS42406 = 'Greens' ;
cmap_struct.SS42456 = 'Oranges' ;
cmap_struct.SS42461 = 'PuBuGn' ;
cmap_struct.SS43993 = 'YlOrRd' ;
cmap_struct.SS45813 = 'RdPu' ;
cmap_struct.SS47195 = 'Purples' ;

%--------------------------------------------------------------------------
%% more descriptive names for driver lines
driver_name_struct = struct() ; 
% motor lines
driver_name_struct.SS37252 = 'SS37252 (hg2)' ; 
driver_name_struct.SS37253 = 'SS37253 (hg2)' ; 
driver_name_struct.SS01062 = 'SS01062 (ctrl)' ; 
driver_name_struct.SS37246 = 'SS37246 (i2)' ; 
driver_name_struct.SS37231 = 'SS37231 (i2 weak)' ; 
driver_name_struct.SS40864 = 'SS40864 (hg3)' ; 
driver_name_struct.SS45772 = 'SS45772 (i1 + hg1)' ;
driver_name_struct.SS45778 = 'SS45778 (weak i2)' ;
driver_name_struct.SS45782 = 'SS45782 (i1 + i2)' ;
driver_name_struct.SS47120 = 'SS47120 (tp2)' ;
driver_name_struct.SS47152 = 'SS47152 (ps1)' ;
driver_name_struct.SS01055 = 'SS01055 (ctrl)' ;
driver_name_struct.SS47125 = 'SS47125 (TTMn + hg1)' ;
driver_name_struct.SS34781 = 'SS34781 (i2)' ;
driver_name_struct.SS40980 = 'SS40980 (b1 + hg4)' ;
driver_name_struct.SS41068 = 'SS41068 (DVM)' ;
driver_name_struct.SS47161 = 'SS47161 (hg1)' ;
driver_name_struct.SS48311 = 'SS48311 (hg1)' ;
driver_name_struct.b1      = 'MB258C (b1)' ;
driver_name_struct.SS41039 = 'SS41039 (i1)' ;
driver_name_struct.SS49039 = 'SS49039 (hg3)' ;

% sensory lines
driver_name_struct.SS36076 = 'SS36076 (haltere MN)' ;
driver_name_struct.SS40775 = 'SS40775 (wing mech)' ;
driver_name_struct.SS42406 = 'SS42406 (haltere mech)' ;
driver_name_struct.SS42456 = 'SS42456 (wing mech)' ;
driver_name_struct.SS42461 = 'SS42461 (wing mech)' ;
driver_name_struct.SS43993 = 'SS43993 (wing + haltere mech)' ; 
driver_name_struct.SS45813 = 'SS45813 (wing mech)' ;
driver_name_struct.SS47195 = 'SS47195 (haltere MN)' ;


% ----------------------------------------------------------------
%% brief names for driver 
mn_name_struct.SS01055 = 'ctrl' ;
mn_name_struct.SS01062 = 'ctrl' ; 

mn_name_struct.SS34781 = 'i2' ;
mn_name_struct.SS37252 = 'hg2' ; 
mn_name_struct.SS37253 = 'hg2' ; 
mn_name_struct.SS37246 = 'i2' ; 
mn_name_struct.SS37231 = 'weak i2' ; 

% mn_name_struct.SS40864 = 'SS40864 (hg3)' ; 
% mn_name_struct.SS40980 = 'b1 + hg4' ;
mn_name_struct.SS41039 = 'i1' ;
mn_name_struct.SS41068 = 'DVM' ;
mn_name_struct.SS44056 = 'DLM' ; 
mn_name_struct.SS45772 = 'i1 + hg1' ;
mn_name_struct.SS45778 = 'weak i2' ;
mn_name_struct.SS45782 = 'i1 + i2' ;
mn_name_struct.SS47120 = 'tp2' ;
mn_name_struct.SS47152 = 'ps1' ;
mn_name_struct.SS47125 = 'tt + hg1' ;
mn_name_struct.SS47161 = 'hg1' ;
mn_name_struct.SS48311 = 'hg1' ;
mn_name_struct.SS49039 = 'hg3' ;

%--------------------------------------------------------------------------
%% axis limits
t_lim = [0.5 2] ;
amp_ylim = [-10, 8] ;
fwd_dev_ylim = [-4 4] ; 
back_dev_ylim = [-4 4] ;
wbf_ylim = [-5 10] ; 

ylim_list = [amp_ylim ; back_dev_ylim ; fwd_dev_ylim ; wbf_ylim] ; 

%--------------------------------------------------------------------------
%% figure sizing
figUnits = 'inches' ; 
%fullFigPosition = [7, 2, 7.20500, 9.72400] ; 
%halfFigPosition = [7, 2, 7.20500, 9.72400] ; 
%activityUnits = 'inches' ; 
%activityPosition = [7 2 2*4.6 3.5] ;

axisLineWidth = 0.5 ; 
%-------------------------------------------------------------------------
%% figure fonts
%font sizes
% labelFontSize = 8 ;
% axisFontSize = 8 ; 
% titleFontSize = 10 ;
% clabelFontSize = 8 ; 
labelFontSize = 6 ;
axisFontSize = 6 ; 
titleFontSize = 6 ;
clabelFontSize = 6 ; 
legendFontSize = 6 ; 

%font types
fontName = 'Arial' ; 
fontWeight = 'normal' ;
%-------------------------------------------------------------------------

