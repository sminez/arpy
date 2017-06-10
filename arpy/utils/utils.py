'''
Assorted utility functions
'''


def Tex(obj):
    '''
    Convert the string representation of an object to TeX and print.
    '''
    print(obj.__tex__())
