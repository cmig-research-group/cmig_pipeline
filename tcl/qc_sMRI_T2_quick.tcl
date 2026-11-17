#include <tcl.h>

##############################################################################

# script for quality control (QC) of sMRI T2 reg to T1

# Created:  01/29/2020 by Don Hagler
# Last Mod: 01/29/2020 by Don Hagler

# Scripting reference available at:
# https://surfer.nmr.mgh.harvard.edu/fswiki/TkMeditGuide/TkMeditReference/TkMeditScripting

##############################################################################

# set intensity range
SetVolumeMinMax 0 0 400
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

SetOrientation 2
SetSlice 135 ; RedrawScreen;

