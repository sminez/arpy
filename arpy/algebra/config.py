"""
arpy (Absolute Relativity in Python)
Copyright (C) 2016-2018 Innes D. Anderson-Morrison All rights reserved.
"""


class ARConfig:
    """The arpy paramater configuration object"""

    def __init__(self, allowed, metric, div):
        """Bind in all parameters"""
        self._allowed = allowed
        self.original_allowed = allowed
        self._metric = self._convert_metric(metric)
        self.original_metric = metric
        self.division_type = div

        # Generate the config and bind to the calling scope
        self.update_config()
        # update_env in the __init__

    def __eq__(self, other):
        if not isinstance(other, ARConfig):
            return False

        metric = self.metric == other.metric
        allowed = self.allowed == other.allowed
        div = self.division_type == other.division_type

        return metric and allowed and div

    def __repr__(self):
        metric = "".join("+" if m == 1 else "-" for m in self.metric)
        allowed = "{α" + ", α".join(self.allowed) + "}"
        return (
            "Config Details:\n"
            "===============\n"
            "Allowed Alphas: {}\n"
            "Metric:         {}\n"
            "Division type:  {}"
        ).format(allowed, metric, self.division_type)

    def reset(self):
        """Reset the metric and allowed to their default values"""
        self.allowed = self.original_allowed
        self.metric = self.original_metric
        self.update_config()
        self.update_env(lvl=3)  # See arpy __init__ for details

    def _convert_metric(self, signs):
        """Convert the supplied metric to a tuple of ints"""
        if all(sign in ["+", "-"] for sign in signs):
            if len(signs) != 4:
                raise ValueError(
                    "metric should be a 4 element string.\n" "i.e. 'ar.metric = \"+---\"'"
                )
            metric = tuple(1 if s == "+" else -1 for s in signs)
        elif all(sign in [1, -1] for sign in signs):
            if len(signs) != 4:
                raise ValueError("Metric should be a 4-tuple of 1/-1")
            metric = signs
        else:
            raise ValueError('Invalid metric: {}\nValid examples: "+---", "(1,-1,-1,-1)"')

        return metric

    @property
    def metric(self):
        return self._metric

    @metric.setter
    def metric(self, signs):
        self._metric = self._convert_metric(signs)
        self.update_config()
        self.update_env(lvl=3)  # See arpy __init__ for details

    @property
    def allowed(self):
        return self._allowed

    @allowed.setter
    def allowed(self, allowed):
        if len(allowed) != 16:
            raise ValueError("Must provide all 16 elements for allowed")
        if not all([set(c).issubset(set("p0123")) for c in allowed]):
            raise ValueError("Invalid indices for allowed: {}".format(allowed))

        self._allowed = allowed
        self.update_config()
        self.update_env(lvl=3)  # See arpy __init__ for details

    def update_config(self):
        """
        Define algebra level data and mappings.
        NOTE: The ARConfig class is extended to include an `update_env` method
        in the main __init__.py.
        """
        self._h = [a for a in self.allowed if len(a) == 3 and "0" not in a][0]
        self._q = [a for a in self.allowed if len(a) == 4][0]
        self._B = [a for a in self.allowed if len(a) == 2 and "0" not in a]
        self._T = [a for a in self.allowed if len(a) == 3 and "0" in a]
        self._A = [a for a in self.allowed if len(a) == 1 and a not in "p0"]
        self._E = [a for a in self.allowed if len(a) == 2 and "0" in a]

        # Fast lookup of zet components in {e,x,y,z} order
        _dims = "e x y z".split()
        self.zet_comps = {
            "B": dict(zip(_dims, ["p"] + self._B)),
            "T": dict(zip(_dims, ["0"] + self._T)),
            "A": dict(zip(_dims, [self._h] + self._A)),
            "E": dict(zip(_dims, [self._q] + self._E)),
        }

        e_key = "0i" if self._E[0][0] == "0" else "i0"

        # How the 3-vector components are grouped and under what names
        # TODO: Have a way to dynamically alter these names?
        self.xi_groups = {"i": self._A, e_key: self._E, "jk": self._B, "0jk": self._T}

        self.group_to_zet = {"jk": "B", "i": "A", "0jk": "T", e_key: "E"}

        # Names to group the results of calculations under: scalars & 3-vectors
        self.allowed_groups = ["p", "0", "123", "0123"] + [g for g in self.xi_groups.keys()]


# The labelling and ordering of the 16 elements of the algebra.
# NOTE:: The order will affect the visualisation of the Cayley Table
#       but not the results of finding products.
_B = ["p", "23", "31", "12"]  # ΞB :: Magnetic Field and rest mass
_T = ["0", "023", "031", "012"]  # ΞΤ :: Angular-Momentum/Charge density
_A = ["123", "1", "2", "3"]  # ΞΑ :: Current Density and hedgehog
_E = ["0123", "01", "02", "03"]  # ΞE :: Electric Field and dual rest mass
allowed = _B + _T + _A + _E

# The space-time metric that will be used
metric = (1, -1, -1, -1)

config = ARConfig(allowed, metric, "into")
