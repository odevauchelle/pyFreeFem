
import sys
sys.path.append('./../')

import pyFreeFem as pyff

edp = '''
real a, b;
cin >> a;
cin >> b;
cout << a*b << endl;
'''
print( pyff.run_FreeFem(edp, stdin = [2,3.5]) )
