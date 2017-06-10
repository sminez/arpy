'''
Assorted utility functions
'''


def TexPrint(obj):
    '''
    Convert the string representation of an object to TeX and print.
    '''
    tex_subs = {
        'α': '\\alpha', '∂': '\\partial', 'ξ': '\\xi',
        'ₚ': '_p', '₀': '_0', '₁': '_1', '₂': '_2', '₃': '_3'}

    if hasattr(obj, '__tex__'):
        print(obj.__tex__())
    else:
        _repr = obj.__repr__()

        for char in tex_subs:
            _repr = _repr.replace(char, tex_subs[char])

        # NOTE: need to sort out things like _1_2_3... etc

        print(_repr)
