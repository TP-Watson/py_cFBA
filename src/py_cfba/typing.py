"""Type aliases for the Python cFBA Toolbox."""

__all__ = []

from os import PathLike
from pathlib import Path
from typing import Literal, NamedTuple, Optional, Type, TypedDict, Union

import optlang.interface
from numpy.typing import NDArray

# solver-specific types
from optlang import Constraint as _Constraint
from optlang import Model as _Model
from optlang import Variable as _Variable

# solver-agnostic types
Constraint: Type[optlang.interface.Constraint] = _Constraint
Model: Type[optlang.interface.Model] = _Model
Variable: Type[optlang.interface.Variable] = _Variable


FileName = Union[str, Path, PathLike]

ConstraintType = Literal["equality", "min", "max"]
Quota = tuple[ConstraintType, str, int, float]


class AlphaOutput(NamedTuple):
    alpha: float
    prob: optlang.interface.Model


class CapacityMatrices(NamedTuple):
    A: NDArray
    B: NDArray


class KineticParamsBounds(NamedTuple):
    low_b_var: NDArray
    upp_b_var: NDArray


class FluxOutput(NamedTuple):
    fluxes: NDArray
    amounts: NDArray
    t: NDArray


class ExtractedImbalancedMetabolites(NamedTuple):
    indices_balanced: list[int]
    indices_imbalanced: list[int]
    imbalanced_mets: list[str]
    balanced_mets: list[str]
    w: NDArray
    Sb: NDArray
    Si: NDArray


class InitSMatrix(NamedTuple):
    S: NDArray
    mets: list[str]
    rxns: list[str]


class LPProblemOutput(NamedTuple):
    cons: list[optlang.interface.Constraint]
    Mk: NDArray
    imbalanced_mets: list[str]
    nm: int
    nr: int
    nt: int


class SpeciesData(TypedDict):
    compartment: str
    imbalanced: bool
    w_contribution: Optional[float]


SpeciesDict = dict[str, SpeciesData]


class ReactionData(TypedDict):
    reactants: dict[str, float]
    products: dict[str, float]
    kinetic_law: dict[str, float]
    annotation: str


ReactionDict = dict[str, ReactionData]
