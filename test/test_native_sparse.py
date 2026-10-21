"""Behavioral failure and actual binary dependency checks for native-only use."""
import re
import subprocess
import sys

if sys.argv[1]=='--linkage':
    library=sys.argv[2]
    result=subprocess.run(['ldd',library],capture_output=True,text=True,check=True)
    assert not re.search(r'lib(?:umfpack|superlu|cholmod|[sdcz]?mumps)',result.stdout,re.I),result.stdout
    symbols=subprocess.run(['nm','-D','--undefined-only',library],capture_output=True,
                           text=True,check=True).stdout
    assert not re.search(r'\b(?:umf4|umfpack_|cholmod_|[sdcz]mumps|[sdcz]gssv)',symbols,re.I),symbols
    print('PASS native-only binary linkage')
elif sys.argv[1]=='--freefem-singular':
    result=subprocess.run([sys.argv[2],'-nw','-nowait',sys.argv[3]],
                          capture_output=True,text=True,timeout=1.5)
    output=result.stdout+result.stderr
    assert result.returncode!=0 and 'FortSparse native FreeFem failure 1' in output,output
    assert 'UNEXPECTED SINGULAR SUCCESS' not in output.split('sizestack')[-1],output
    print('PASS FreeFem native singular failure without fallback')
else:
    for mode,diagnostic in [('singular','FortSparse native failure 1'),
                            ('invalid','invalid sparse CSC pointer')]:
        result=subprocess.run([sys.argv[1],mode],capture_output=True,text=True,timeout=.8)
        assert result.returncode!=0 and diagnostic in result.stdout+result.stderr, result
    print('PASS native singular/input failures without fallback')
