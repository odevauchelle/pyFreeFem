from numpy import integer, floating, ndarray
from json import JSONEncoder, dumps

###################################
#
# Export mesh to json
#
####################################

class NpEncoder( JSONEncoder ):
    '''
    Encoder from numpy objects to json.

    Based on:
    https://stackoverflow.com/questions/50916422/python-typeerror-object-of-type-int64-is-not-json-serializable
    '''

    def default(self, obj):

        if isinstance( obj, integer ):
            return int(obj)

        if isinstance( obj, floating ):
            return float(obj)

        if isinstance( obj, ndarray ):
            return obj.tolist()

        return super(NpEncoder, self).default(obj)



def export_to_json( Th, keys = None ) :
    '''
    Exports mesh to json string.
    '''

    if keys is None :
        keys = [ 'x', 'y', 'triangles', 'node_labels', 'triangle_labels' ]

    Th_dict = {}

    for key in keys :
        Th_dict[key] = Th.__dict__[key]

    Th_dict['boundary_edges'] = []

    for key, value in Th.boundary_edges.items() : # because json objects can't handle tuple keys

        Th_dict['boundary_edges'] += [ list(key) + [value] ]

    return dumps( Th_dict, cls = NpEncoder )
