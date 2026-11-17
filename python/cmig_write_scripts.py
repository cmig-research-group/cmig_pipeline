"""
Create batch job scripts for MATLAB and cluster execution.
Handles multiple data types, including vectors and cell arrays.

Created:  12/22/25 by Don Hagler
Prev Mod: 01/07/26 by Don Hagler
Prev Mod: 06/12/26 by Don Hagler
Prev Mod: 06/15/26 by Don Hagler
Prev Mod: 07/24/26 by Don Hagler
Prev Mod: 07/28/26 by Don Hagler
Last Mod: 09/14/26 by Don Hagler
"""

import os
from typing import List, Dict, Any, Optional, Union


def cmig_write_singularity_script(
    fname: str,
    cmd: str,
    args: List[Any],
    mrtdir: str,
    appdir: str,
    indir: Union[str, List[str]],
    outdir: Union[str, List[str]],
    tags: Optional[List[str]] = None,
    parms: Optional[Dict[str, Any]] = None,
    appendflag: bool = False,
    header: str = '#!/bin/csh',
    appver: str = '',
    pubsdir: str = '/usr/pubsw'
):
    """
    Create csh script to run singularity image with runtime library compiled app
    for batch processing. Handles multiple data types, including vectors and cell arrays.

    Parameters
    ----------
    fname : str
        Output file name
    cmd : str
        MATLAB function name
    args : list
        List of arguments to come before ('key', value) pairs
    mrtdir : str
        MATLAB runtime library root directory
    appdir : str
        Compiled MATLAB application root directory
    indir : str or list of str
        Input directory(ies) to be mounted as read-only
    outdir : str or list of str
        Output directory(ies) to be mounted as read-writeable
    tags : list of str, optional
        List of names of fields in parms (the 'key's). Default is None
    parms : dict, optional
        Dictionary containing parameter values (the values). Default is None
    appendflag : bool, optional
        Whether to append to existing file. Default is False
    header : str, optional
        String containing line or lines to appear at beginning of script.
        Default is '#!/bin/csh'
    appver : str, optional
        Version string for singularity image of cmd. If not supplied, will use
        first image with name matching cmd. Default is ''
    pubsdir : str, optional
        Shared software directory root directory. Default: /usr/pubsw

    Raises
    ------
    FileNotFoundError
        If singularity image not found
    IOError
        If failed to open file for writing
    """
    import glob
    
    if tags is None:
        tags = []
    if parms is None:
        parms = {}
    
    # Ensure args, indir, outdir are lists
    if args is not None and not isinstance(args, list):
        args = [args]
    elif args is None:
        args = []
    
    if indir is not None and not isinstance(indir, list):
        indir = [indir]
    elif indir is None:
        indir = []
    
    if outdir is not None and not isinstance(outdir, list):
        outdir = [outdir]
    elif outdir is None:
        outdir = []
    
    # Set subdir for images
    subdir = os.path.join(appdir, 'images')
    if not os.path.isdir(subdir):
        subdir = appdir
    
    # Get image file
    pattern = os.path.join(subdir, f'{cmd.lower()}*{appver}.sif')
    flist = glob.glob(pattern)

    if not flist:
        raise FileNotFoundError(
            f'singularity image for cmd {cmd.lower()} not found in {subdir}'
        )
    
    imgfile = flist[0]
    
    # Open file for writing or appending
    mode = 'a' if appendflag else 'w'
    try:
        with open(fname, mode) as fid:
            # Write header
            if header:
                fid.write(f'{header}\n\n')
            
            # Write commands to create output directories
            for d in outdir:
                fid.write(f'mkdir -p {d}\n')
            
            # Write command to run singularity
            fid.write('\nsingularity run -e \\\n')
            
            # Write singularity options to mount directories
            fid.write(f"  -B {pubsdir}':'/usr/pubsw':'ro \\\n")
            fid.write(f"  -B {mrtdir}':'/opt/matlabruntime/R2023a':'ro \\\n")
            
            for d in indir:
                fid.write(f"  -B {d}':'{d}':'ro \\\n")
            
            for d in outdir:
                fid.write(f"  -B {d}':'{d} \\\n")
            
            # Write image file name
            fid.write(f'  {imgfile} \\\n')
            
            # Write args
            for a, arg in enumerate(args):
                fid.write('  ')
                cmig_write_value_csh(fid, arg)
                if a < len(args) - 1 or tags:
                    fid.write(' ')
                fid.write(' \\\n')
            
            # Write tags
            if parms:
                cmig_write_tags_csh(fid, tags, parms)
            
            fid.write('\n')
    
    except IOError as e:
        raise IOError(f'failed to open file {fname} for writing: {e}')


def cmig_write_runtime_script(
    fname: str,
    cmd: str,
    args: List[Any],
    mrtdir: str,
    appdir: str,
    tags: Optional[List[str]] = None,
    parms: Optional[Dict[str, Any]] = None,
    appendflag: bool = False,
    header: str = '#!/bin/csh -f',
    appver: str = ''
):
    """
    Create csh script to run matlab runtime library compiled app for batch processing.
    Handles multiple data types, including vectors and cell arrays.
    
    Parameters
    ----------
    fname : str
        Output file name
    cmd : str
        MATLAB function name
    args : list
        List of arguments to come before ('key', value) pairs
    mrtdir : str
        MATLAB runtime library root directory
    appdir : str
        Compiled MATLAB application root directory
    tags : list of str, optional
        List of names of fields in parms (the 'key's). Default is None
    parms : dict, optional
        Dictionary containing parameter values (the values). Default is None
    appendflag : bool, optional
        Whether to append to existing file. Default is False
    header : str, optional
        String containing line or lines to appear at beginning of script.
        Default is '#!/bin/csh -f'
    appver : str, optional
        version string for compiled subdirectory of appdir
        ignored if {appdir}/{cmd} exists
        if supplied, will look in {appdir}/compiled_{appver}
        (leading mmps_ stripped, remainder lowercased: MMPS_DJH -> compiled_djh)
        if not supplied, will look in {appdir}/compiled
        Default is ''

    Raises
    ------
    FileNotFoundError
        If app subdirectory not found
    IOError
        If failed to open file for writing
    """
    import re

    if tags is None:
        tags = []
    if parms is None:
        parms = {}

    # Ensure args is a list
    if args is not None and not isinstance(args, list):
        args = [args]
    elif args is None:
        args = []

    # Set cmdstr based on appdir, appver, and cmd
    subdir = os.path.join(appdir, cmd)
    if not os.path.isdir(subdir):
        if appver:
            stripped = re.sub(r'^mmps_', '', appver.lower())
            compdir = os.path.join(appdir, f'compiled_{stripped}')
        else:
            compdir = os.path.join(appdir, 'compiled')
        subdir = os.path.join(compdir, cmd)

    if not os.path.isdir(subdir):
        raise FileNotFoundError(f'app subdirectory {subdir} not found')
    
    cmdstr = os.path.join(subdir, f'run_{cmd}.sh')
    
    # Open file for writing or appending
    mode = 'a' if appendflag else 'w'
    try:
        with open(fname, mode) as fid:
            # Write header
            if header:
                fid.write(f'{header}\n\n')
            
            # Write command to run app
            fid.write(f'{cmdstr} \\\n')
            fid.write(f"'{mrtdir}' \\\n")
            
            # Write args
            for a, arg in enumerate(args):
                cmig_write_value_csh(fid, arg)
                if a < len(args) - 1 or tags:
                    fid.write(' ')
                fid.write(' \\\n')
            
            # Write tags
            if parms:
                cmig_write_tags_csh(fid, tags, parms)
            
            fid.write('\n')
    
    except IOError as e:
        raise IOError(f'failed to open file {fname} for writing: {e}')


def cmig_write_value_csh(fid, value: Any, nested: bool = False):
    """
    Write a value to open file specified by fid.
    Handles multiple data types, including vectors and cell arrays.

    A multi-element array is written inside double quotes, as a cell array is,
    so csh passes it as one argument instead of splitting it on the spaces.
    An array nested in a cell array is already inside the cell's quotes and
    takes none of its own.
    
    Parameters
    ----------
    fid : file object
        Open file object
    value : any
        Value to write
    
    Raises
    ------
    TypeError
        If unrecognized type
    """
    from datetime import datetime
    import math
    
    # Check if value is a float or list of floats
    if isinstance(value, float) or (isinstance(value, (list, tuple)) and 
                                     len(value) > 0 and 
                                     all(isinstance(v, float) for v in value)):
        # Handle float arrays
        values = [value] if isinstance(value, float) else list(value)
        q = '' if nested else '"'

        if len(values) != 1:
            fid.write(q + '[')
        
        for v in values:
            if math.isinf(v):
                if v < 0:
                    fid.write('-Inf')
                else:
                    fid.write('Inf')
            elif v == int(v):  # Check if float is actually an integer value
                fid.write(f'{int(v)} ')
            else:
                fid.write(f'{v:.8f} ')
        
        if len(values) != 1:
            fid.write(']' + q)
    
    # Check if value is numeric (int) or logical (bool) or list of them
    elif (isinstance(value, (int, bool)) or 
          (isinstance(value, (list, tuple)) and len(value) > 0 and 
           all(isinstance(v, (int, bool)) for v in value))):
        
        values = [value] if isinstance(value, (int, bool)) else list(value)
        q = '' if nested else '"'

        if len(values) == 1:
            fid.write(f'{int(values[0])}')
        else:
            fid.write(q + '[')
            fid.write(''.join(f'{int(v)} ' for v in values))
            fid.write(']' + q)
    
    # Check if value is a string
    elif isinstance(value, str):
        fid.write(f"'{value}'")
    
    # Check if value is empty
    elif value is None or (isinstance(value, (list, tuple)) and len(value) == 0):
        fid.write(' ')
    
    # Check if value is a cell array (list with mixed types)
    elif isinstance(value, (list, tuple)) and len(value) > 0:
        fid.write('"{')
        for k, v in enumerate(value):
            cmig_write_value_csh(fid, v, nested=True)
            if k < len(value) - 1:
                fid.write(' ')
        fid.write('}"')
    
    # Check if value is a datetime
    elif isinstance(value, datetime):
        # Format datetime similar to MATLAB's datestr
        fid.write(f"'{value.strftime('%d-%b-%Y %H:%M:%S')}'")
    
    else:
        raise TypeError('unrecognized type')


def cmig_write_tags_csh(fid, tags: List[str], parms: Dict[str, Any]):
    """
    Write ('key', value) pairs to open file specified by fid.
    Handles multiple data types, including vectors and cell arrays.
    
    Parameters
    ----------
    fid : file object
        Open file object
    tags : list of str
        List of parameter names (keys)
    parms : dict
        Dictionary containing parameter values
    """
    
    if not tags:
        tags = list(parms.keys())

    for t, tag in enumerate(tags):
        value = parms.get(tag, None)

        fid.write(f"  '{tag}' ")
        cmig_write_value_csh(fid, value)

        if t < len(tags) - 1:
            fid.write(' \\\n')


def cmig_write_value_matlab(fid, value):
    """Write a MATLAB-syntax value to open file fid.

    Analogous to mmil_write_value.m.  Handles float, int, bool, str,
    None/empty, numeric arrays, cell arrays, and datetime.

    Parameters
    ----------
    fid : file object
        Open file object
    value : any
        Value to write

    Raises
    ------
    TypeError
        If unrecognized type
    """
    import math
    from datetime import datetime

    # None or empty string or empty sequence -> []
    if (value is None or value == '' or
            (isinstance(value, (list, tuple)) and len(value) == 0)):
        fid.write('[]')

    # bool (must precede int since bool is a subclass of int)
    elif isinstance(value, bool):
        fid.write('1' if value else '0')

    # int
    elif isinstance(value, int):
        fid.write('%d' % value)

    # float
    elif isinstance(value, float):
        if math.isnan(value):
            fid.write('NaN')
        elif math.isinf(value):
            fid.write('-Inf' if value < 0 else 'Inf')
        elif value == int(value):
            fid.write('%d' % int(value))
        else:
            fid.write('%0.8f' % value)

    # string
    elif isinstance(value, str):
        fid.write("'%s'" % value)

    # datetime -> MATLAB datestr format
    elif isinstance(value, datetime):
        fid.write("'%s'" % value.strftime('%d-%b-%Y %H:%M:%S'))

    # list / tuple
    elif isinstance(value, (list, tuple)):
        # numeric array when every element is a non-bool number
        if all(isinstance(v, (int, float)) and not isinstance(v, bool)
               for v in value):
            fid.write('[')
            for k, v in enumerate(value):
                cmig_write_value_matlab(fid, v)
                if k < len(value) - 1:
                    fid.write(' ')
            fid.write(']')
        else:
            # cell array
            fid.write('{')
            for k, v in enumerate(value):
                cmig_write_value_matlab(fid, v)
                if k < len(value) - 1:
                    fid.write(' ')
            fid.write('}')

    else:
        raise TypeError('cmig_write_value_matlab: unrecognized type: %s'
                        % type(value))


def cmig_write_tags_matlab(fid, tags, parms):
    """Write ('key', value) continuation lines to open file fid.

    Analogous to mmil_write_tags.m.  Writes MATLAB-syntax 'key',value pairs
    as continuation lines inside an already-opened function call.

    Parameters
    ----------
    fid : file object
        Open file object
    tags : list of str
        List of parameter names (keys)
    parms : dict
        Dictionary containing parameter values
    """
    if not tags:
        tags = list(parms.keys())

    for t, tag in enumerate(tags):
        value = parms.get(tag, None)
        fid.write("  '%s'," % tag)
        cmig_write_value_matlab(fid, value)
        if t < len(tags) - 1:
            fid.write(',...\n')


def cmig_write_matlab_script(
    fname,
    cmd,
    args,
    tags=None,
    parms=None,
    quitflag=True,
    appendflag=False,
    header=None,
):
    """Write a MATLAB .m script that calls a function for batch processing.

    Analogous to mmil_write_script.m.  Writes argument values inline in the
    function call, followed by ('key', value) pairs for any tags/parms.

    Generated script structure::

        [header]
        cmd(...
          value1,...
          value2,...
          'key1',value1,...
          'key2',value2);
        [exit]

    Parameters
    ----------
    fname : str
        Output file path
    cmd : str
        MATLAB function name to call
    args : list
        Positional argument values
    tags : list of str, optional
        Names of key-value parameters (the 'key' strings). Default: None
    parms : dict, optional
        Key-value parameter dictionary. Default: None
    quitflag : bool, optional
        Add 'exit' statement at end of script. Default: True
    appendflag : bool, optional
        Append to existing file instead of overwriting. Default: False
    header : str, optional
        Line(s) to write at beginning of script. Default: None

    Raises
    ------
    IOError
        If the file cannot be opened for writing
    """
    if tags is None:
        tags = []
    if parms is None:
        parms = {}
    if args is None:
        args = []
    if not isinstance(args, list):
        args = [args]

    mode = 'a' if appendflag else 'w'
    try:
        with open(fname, mode) as fid:
            if header:
                fid.write('%s\n\n' % header)

            fid.write('%s(...\n' % cmd)

            for a, arg in enumerate(args):
                fid.write('  ')
                cmig_write_value_matlab(fid, arg)
                if a < len(args) - 1 or tags or parms:
                    fid.write(',')
                fid.write('...\n')

            if tags or parms:
                cmig_write_tags_matlab(fid, tags, parms)

            fid.write(');\n')
            if quitflag:
                fid.write('exit\n')

    except IOError as e:
        raise IOError('failed to open file %s for writing: %s' % (fname, e))


# Example usage
if __name__ == '__main__':
    # Example parameters
    fname = 'test_script.csh'
    cmd = 'my_function'
    args = ['input.txt', 42, [1.5, 2.7, 3.9]]
    mrtdir = '/usr/local/MATLAB/MATLAB_Runtime/v99'
    appdir = '/path/to/compiled/apps'
    tags = ['thresh', 'outdir', 'verbose']
    parms = {
        'thresh': 0.5,
        'outdir': '/output/path',
        'verbose': True
    }
    
    cmig_write_runtime_script(
        fname, cmd, args, mrtdir, appdir,
        tags=tags, parms=parms, appendflag=False
    )
