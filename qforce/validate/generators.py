from colt import Colt
from abc import ABC, abstractmethod


_implemented_generators = ["annealing", "crest", "qcg"]


def _implemented_annealers():
    annealers = {"gromacs": GromacsAnnealing, "xtb": XtbAnnealing}
    return annealers.items()


class AnnealerABC(ABC, Colt):
    _user_input = """
    # 
    t_min = 200. :: float

    # 
    t_max = 500. :: float

    #
    n_steps = 10 :: int
    """

    def __init__(self, settings):
        self.settings = settings
        self._t_min = settings.annealing.t_min
        self._t_max = settings.annealing.t_max
        self._n_steps = settings.annealing.n_steps

    @abstractmethod
    def generate_input_files(self): ...

    @property
    def t_min(self):
        return self._t_min

    @property
    def t_max(self):
        return self._t_max

    @property
    def n_steps(self):
        return self._n_steps


class GromacsAnnealing(AnnealerABC):
    _user_input = """
    #
    test_input = sometest :: str
    """

    def generate_input_files(self):
        print(self._t_min)


class XtbAnnealing(AnnealerABC):
    _user_input = """
    #
    test_input = sometest :: str
    """

    def generate_input_files(self):
        print(self._t_min)
