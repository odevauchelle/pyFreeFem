import sys
sys.path.append('./../')
import pyFreeFem as pyff

script = pyff.edpScript('''
border Circle( t = 0, 2*pi ){ x = cos(t); y = sin(t); }
mesh Th = buildmesh( Circle(20) );
''' )

script.pprint()

script += 'plot(Th);'

script.run()
