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


def Zet(alpha):
    '''Return the Zet of a given alpha value.'''
    # Conditions for being a member of each Zet
    zet_map = {
        # `e` elements (NOTE: `p` is a special case)
        (1, True): 'T', (3, False): 'A', (4, True): 'E',
        # `x, y, z` elements: (len, has '0')
        (2, False): 'B', (3, True): 'T', (1, False): 'A', (2, True): 'E',
    }
    ix = alpha.index

    if ix == 'p':
        return 'B'
    else:
        return zet_map[(len(ix), 0 in ix)]


def Nat(alpha):
    '''Return the Nature of a given Alpha.'''
    # Element sets for each e,x,y,z nature
    nat_map = {
        set('p'): 'e', set('123'): 'e', set('0'): 'e', set('0123'): 'e',
        set('1'): 'x', set('23'): 'x', set('023'): 'x', set('01'): 'x',
        set('2'): 'y', set('31'): 'y', set('031'): 'y', set('02'): 'y',
        set('3'): 'z', set('12'): 'z', set('012'): 'z', set('03'): 'z',
    }

    return nat_map[set(alpha.index)]


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
