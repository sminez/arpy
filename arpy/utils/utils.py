'''
Assorted data structures and utility functions
'''


SUPER_SCRIPTS = {'B': 'ᴮ', 'A': 'ᴬ', 'T': 'ᵀ', 'E': 'ᴱ'}
SUB_SCRIPTS = {'0': '₀', '1': '₁', '2': '₂', '3': '₃', 'p': 'ₚ'}


def Tex(obj):
    '''
    Convert the string representation of an object to TeX and print.
    '''
    print(obj.__tex__())


def reorder_allowed(allowed, order='pBtThAqE'):
    '''
    Shuffle the ordering of allowed, keeping 3-Vectors together.
    NOTE: This assumes that the input is in pBtThAqE order to start.
    '''
    groups = {
        'p': ['p'], 't': ['0'], 'A': ['1', '2', '3'],
        'B': allowed[1:4], 'T': allowed[5:8], 'E': allowed[13:16],
        'q': [allowed[12]], 'h': [allowed[8]]
    }
    new = []
    for group in order:
        new += groups[group]
    return new
