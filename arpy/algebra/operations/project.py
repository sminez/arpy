from ...utils.concepts.dispatch import dispatch_on
from ..config import config as cfg
from ..data_types import Alpha, MultiVector, Term
from .full import POINT


@dispatch_on(index=0)
def project(element, grade, cfg=cfg):
    """
    Implementation of the grade-projection operator <A>n.
    Return only the elements of A that are of grade n.
    NOTE:: αp is a grade-0 scalar element.
    """
    raise NotImplementedError


@project.add(Alpha)
def _project_alpha(element, grade, cfg=cfg):
    if element._index == POINT:
        if grade == 0:
            return element
        else:
            return None
    elif len(element._index) == grade:
        return element
    else:
        return None


@project.add(Term)
def _project_term(element, grade, cfg=cfg):
    if element.index == POINT:
        if grade == 0:
            return element
        else:
            return None
    elif len(element.index) == grade:
        return element
    else:
        return None


@project.add(MultiVector)
def _project_multivector(element, grade, cfg=cfg):
    correct_grade = []
    if grade == 0:
        for term in element:
            if term.index == POINT:
                correct_grade.append(term)
    else:
        for term in element:
            ix = term.index
            if len(ix) == grade and ix != POINT:
                correct_grade.append(term)
    res = MultiVector(correct_grade, cfg=cfg)
    return res
