#include <tcl.h>

##############################################################################

# script for quality control (QC) of DTI FA reg to T1 and AtlasTrack fseg

# Created:  01/22/2016 by Don Hagler
# Last Mod: 01/09/2020 by Don Hagler

# Scripting reference available at:
# https://surfer.nmr.mgh.harvard.edu/fswiki/TkMeditGuide/TkMeditReference/TkMeditScripting

##############################################################################

# set intensity range
SetVolumeMinMax 0 0 1
SetVolumeMinMax 1 0 255

# volume brightness contrast 0=main 1=aux
SetVolumeBrightnessContrast 0 .35 12
SetVolumeBrightnessContrast 1 .35 12

# set zoom
SetZoomLevel 1

# turn off overlay and segmentation
SetDisplayFlag 1 0
SetDisplayFlag 15 0
SetDisplayFlag 16 0

# check corpus callosum registration
SetOrientation 2
for { set slice 125 } { $slice < 135 } { incr slice } {
    SetSlice $slice
    SetDisplayFlag 1 1
    RedrawScreen
    after 150
    SetDisplayFlag 1 0
    RedrawScreen
    after 150
}

# scan through brain
SetOrientation 2
for { set slice 75 } { $slice < 165 } { incr slice } {
  SetSlice $slice
  RedrawScreen
  after 25
}

# display fseg in each plane
SetDisplayFlag 15 1
SetDisplayFlag 16 1
SetOrientation 1
SetSlice 120; RedrawScreen;
after 500;
SetOrientation 0
SetSlice 145; RedrawScreen;
after 500
SetOrientation 2
SetSlice 135; RedrawScreen;
after 500

# turn off fseg
SetDisplayFlag 15 0
SetDisplayFlag 16 0

