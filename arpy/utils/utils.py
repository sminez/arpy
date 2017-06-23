'''
Assorted data structures and utility functions
'''


SUPER_SCRIPTS = {'B': 'ᴮ', 'A': 'ᴬ', 'T': 'ᵀ', 'E': 'ᴱ'}
SUB_SCRIPTS = {'0': '₀', '1': '₁', '2': '₂', '3': '₃',
               'p': 'ₚ', 'i': 'ᵢ', 'j': 'ⱼ', 'k': 'ₖ'}


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
    p = ['p']
    t = ['0']
    h = [a for a in allowed if len(a) == 3 and '0' not in a]
    q = [a for a in allowed if len(a) == 4]
    B = [a for a in allowed if len(a) == 2 and '0' not in a]
    T = [a for a in allowed if len(a) == 3 and '0' in a]
    A = [a for a in allowed if len(a) == 1 and a not in ['p', '0']]
    E = [a for a in allowed if len(a) == 2 and '0' in a]

    groups = {'p': p, 't': t, 'h': h, 'q': q,
              'B': B, 'T': T, 'A': A, 'E': E}
    new = []
    for group in order:
        new += groups[group]
    return new
