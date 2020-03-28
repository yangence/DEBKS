#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
author: Zelin Liu
email: zlliu@bjmu.edu.cn
license: GPL3
detail: Analysis differential back-spliced in of rRNA-depleted RNA-seq
"""
from docopt import docopt
import sys
from .version import __version__
#__version__=1.01

helpInfo = '''
Usage: DEBKS <command> [options]

Command:
    merge            Merge circRNA junction from other software
    count            Count linear junction with BAM file
    anno             Annotate circRNA with gene annotation
    dec              Detect differentially expressed circRNA with junction information
'''

def main():
    command_log = 'DEBKS parameters: ' + ' '.join(sys.argv)
    if len(sys.argv) == 1:
        sys.exit(helpInfo)
    elif sys.argv[1] == '--version' or sys.argv[1] == '-v':
        sys.exit(__version__)
    elif sys.argv[1] == 'dec':
        from . import dec
        dec.dec(docopt(dec.__doc__, version=__version__))
    elif sys.argv[1] == 'merge':
        from . import merge
        merge.merge(docopt(merge.__doc__, version=__version__))
    elif sys.argv[1] == 'anno':
        from . import anno
        anno.anno(docopt(anno.__doc__, version=__version__))
    elif sys.argv[1] == 'count':
        from . import count
        count.count(docopt(count.__doc__, version=__version__))
    else:
        sys.exit(helpInfo)

if __name__ == '__main__':
    main()
