#!/usr/bin/env python
"""
cmig_set_abcd_skipTRs.py

Set skipTRs (number of TRs to skip at the start of each run) based on
scanner Manufacturer and SoftwareVersion.

Public function
---------------
set_abcd_skipTRs(manufacturer, software_version=None,
                 skipTRs_GE_DV25=4, skipTRs_GE=15,
                 skipTRs_Philips=8, skipTRs_Siemens=8, verbose=True)
Returns skipTRs (int).

Based on:
  reference/cmig_set_abcd_skipTRs.m
    Created:  05/11/2026 by Don Hagler
    Last Mod: 05/11/2026 by Don Hagler

Created:  07/27/2026 by Don Hagler
Last mod: 07/27/2026 by Don Hagler
"""
from __future__ import print_function
import re


def set_abcd_skipTRs(manufacturer, software_version=None,
                     skipTRs_GE_DV25=4, skipTRs_GE=15,
                     skipTRs_Philips=8, skipTRs_Siemens=8, verbose=True):
    """Return skipTRs for a scan given its Manufacturer/SoftwareVersion.

    manufacturer : str or None
        DICOM/BIDS Manufacturer field (e.g. 'GE MEDICAL SYSTEMS', 'SIEMENS',
        'Philips Medical Systems'). None/NaN/empty is treated as unknown.
    software_version : str or None
        DICOM SoftwareVersion(s) field. Only consulted for GE scanners, to
        distinguish DV25/Orchestra software from older GE software. A
        missing/empty value is treated as "not DV25/Orchestra" (falls
        through to skipTRs_GE), matching cmig_set_abcd_skipTRs.m's own
        graceful degradation for scans whose SoftwareVersion was never
        captured.
    """
    manufacturer = '' if manufacturer is None else str(manufacturer)
    if manufacturer.lower() == 'nan':
        manufacturer = ''
    software_version = '' if software_version is None else str(software_version)
    if software_version.lower() == 'nan':
        software_version = ''

    if re.search('GE', manufacturer, re.IGNORECASE):
        if re.search('DV25|Orchestra', software_version):
            skipTRs = skipTRs_GE_DV25
        else:
            skipTRs = skipTRs_GE
    elif re.search('Philips', manufacturer, re.IGNORECASE):
        skipTRs = skipTRs_Philips
    elif re.search('Siemens', manufacturer, re.IGNORECASE):
        skipTRs = skipTRs_Siemens
    else:
        if verbose:
            print('set_abcd_skipTRs: WARNING: unsupported manufacturer '
                 '{!r}; setting skipTRs = 0'.format(manufacturer))
        skipTRs = 0

    return skipTRs
