'''
arpy (Absolute Relativity in Python)
Copyright (C) 2016-2017 Innes D. Anderson-Morrison All rights reserved.

Tools for visualising results from calculations and investigations.
'''
from ..algebra.config import config
from ..algebra.ar_types import Alpha
from ..algebra.operations import full


def cayley(op=full, padding=6, cfg=config):
    '''
    Print current Cayley table to the terminal allowing for specification
    of the operation used to compute the table.

    Any function that accepts two Alphas can be passed as op.
    '''
    comps = (' '.join(
        [str(op(Alpha(a, cfg=cfg), Alpha(b, cfg=cfg), cfg=cfg)).rjust(padding)
            for b in cfg.allowed])
        for a in cfg.allowed
    )
    for comp in comps:
        print(comp)


def sign_cayley(op=full, cfg=config):
    '''
    Print +-1 signs for the current Cayley table to the terminal allowing
    for specification of the operation used to compute the table.

    Any function that accepts two Alphas can be passed as op.
    '''
    divider = '      ' + ''.join('+---------' for _ in range(4)) + '+'
    comps = (
        ' '.join([
            '■' if op(
                Alpha(x, cfg=cfg), Alpha(y, cfg=cfg), cfg=cfg
            ).sign == -1 else '□'
            for y in cfg.allowed])
        for x in cfg.allowed
    )

    print('          ', '         '.join(['B', 'A', 'T', 'E']))
    print(divider)

    for i, comp in enumerate(comps):
        comp = '| '.join(comp[n:n+8] for n in range(0, len(comp), 8))
        print(str(Alpha(cfg.allowed[i], cfg=cfg)).ljust(5), '|', comp, '|')
        # Divide after each 4-Set
        if (i + 1) % 4 == 0:
            print(divider)


def _4block(rows, cols, op, cfg):
    '''Visualise a 4x4 block of 4 elements acting on 4 others'''
    block = []
    for r in rows:
        comps = [
            op(Alpha(r, cfg=cfg), Alpha(c, cfg=cfg), cfg=cfg).sign
            for c in cols
        ]
        block_row = ' '.join(['□' if c == 1 else '■' for c in comps])
        block.append('|' + block_row + '|')
    return block


def sign_distribution(op=full, cfg=config):
    '''
    By calculating one term for each of the 5 grouped components of the
    cayley table, look at how each metric / allowed set of indices affects
    the overall structure of the algebra.
    '''
    allowed = cfg.allowed
    bs, xs, zs = allowed[0:16:4], allowed[1:16:4], allowed[3:16:4]

    blocks = []

    row_cols = {
        '∂b': (bs, bs), '∂Ξ': (bs, xs),
        '∇': (xs, bs), '∇•': (xs, xs), '∇x': (zs, xs)
    }

    for name in ['∂b', '∂Ξ', '∇', '∇•', '∇x']:
        rows, cols = row_cols[name]
        blocks.append((name, _4block(rows, cols, op, cfg)))

    for i in range(4):
        for block in blocks:
            name = block[0].rjust(3) if i == 0 else '   '
            print(name, block[1][i], end=' ')
        print('')
