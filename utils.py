from typing import Tuple, Union, Iterable, Dict, Optional, List
import uuid
from io import StringIO

from IPython.display import SVG
from openff.toolkit.topology import FrozenMolecule

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.rdDepictor import Compute2DCoords

import numpy as np

import nglview
from nglview.base_adaptor import Structure, Trajectory

from openff.models.models import DefaultModel
from openff.units import unit, Quantity
from pydantic import validator
from pydantic.color import Color

from openff.toolkit import Molecule, ForceField, Topology
from openff.toolkit.typing.engines.smirnoff.parameters import BondType, ParameterType
from openff.toolkit.utils.exceptions import ParameterLookupError


AMBER_COLOR = np.array([255, 192, 0])
SAGE_COLOR = np.array([141, 208, 169])
NEW_PARAM_COLOR = np.array([191, 192, 238])
WHITE = np.array([255, 255, 255])

DEFAULT_WIDTH = 980
DEFAULT_HEIGHT = 400

BondIndices = Tuple[int, int]

MOLECULE_DEFAULT_REPS = [
    dict(type="licorice", params=dict(radius=0.25, multipleBond=True))
]


def in_ff(handler: str, param: ParameterType, ff: ForceField) -> bool:
    """
    True if param is in ff[handler], false if not
    """
    smirks = param.smirks
    try:
        ff_param = ff[handler][smirks]
    except ParameterLookupError:
        return False

    return ff_param.__dict__ == param.__dict__


def bond_in_ff(bond: BondType, ff: ForceField) -> bool:
    return in_ff("Bonds", bond, ff)


def blend_colors(
    a: tuple[float], b: tuple[float], t: float = 0.5, w: float = 0.0
) -> tuple[float]:
    """
    Blend two colors with white.

    By default, produces a color with an even mix of ``a`` and ``b`` and no
    white.

    Parameters
    ==========

    a
        The first color.
    b
        The second color.
    t
        The proportion of ``a`` in the final mix.
    w
        The proportion of white in the final mix.
    """
    a = np.asarray(a) / 255
    b = np.asarray(b) / 255
    blended = np.sqrt(w * WHITE**2 + t * a**2 + (1 - t - w) * b**2)
    return blended * 255


def color_to_floats(color: Color) -> tuple[float, float, float]:
    """The color as a 3-tuple of floats between 0 and 1 representing R, G, B values.

    Alpha is ignored."""
    r, g, b = color.as_rgb_tuple(alpha=False)
    return (r / 255.0, g / 255.0, b / 255.0)


class AtomHighlight(DefaultModel):
    """Color and radius of an atom's highlight."""

    color: Optional[Color] = None
    """The color of the highlight. Alpha is ignored."""
    radius: Optional[Quantity] = None
    """The radius of the highlight, in molecular-scale distance units."""

    @validator("radius")
    def radius_is_distance(cls, value: Optional[Quantity]) -> Optional[Quantity]:
        if value is not None and not value.is_compatible_with(unit.angstrom):
            raise ValueError("Radius should have dimensions of distance")
        return value

    @classmethod
    def _from_color_or_self(cls, obj: Union[Color, "AtomHighlight"]) -> "AtomHighlight":
        if isinstance(obj, cls):
            return obj.copy()
        return cls(color=obj)


def draw_molecule(
    molecule: Union[str, Molecule, Chem.Mol],
    image_width: int = DEFAULT_WIDTH,
    image_height: int = DEFAULT_HEIGHT,
    highlight_atoms: Optional[
        Union[list[int], dict[int, Union[Color, AtomHighlight]]]
    ] = None,
    highlight_bonds: Optional[
        Union[list[BondIndices], dict[BondIndices, Color]]
    ] = None,
    atom_notes: Optional[dict[int, str]] = None,
    bond_notes: Optional[dict[BondIndices, str]] = None,
    explicit_hydrogens: Optional[bool] = None,
    color_by_element: Optional[bool] = None,
) -> SVG:
    """Draw a molecule
    Parameters
    ==========
    molecule
        The molecule to draw.
    image_width
        The width of the resulting image in pixels.
    image_height
        The height of the resulting image in pixels.
    highlight_atoms
        A list of atom indices to highlight, or a map from indices to colors or
        ``AtomHighlights``. Colors are interpreted by Pydantic.
    highlight_bonds
        A list of pairs of atom indices indicating bonds to highlight, or a map
        from index pairs to colors. Colors are interpreted by Pydantic.
    atom_notes
        A map from atom indices to a string that should be printed near the
        atom.
    bond_notes
        A map from atom index pairs to a string that should be printed near the
        bond.
    explicit_hydrogens
        If ``False``, allow uncharged monovalent hydrogens to be hidden. If
        ``True``, make all hydrogens explicit. If ``None``
    color_by_element
        If True, color heteroatoms according to their element; if False, the
        image will be black and white. By default, uses black and white when
        highlight_atoms or highlight_bonds is provided, and color otherwise.
    Raises
    ======
    KeyError
        When an atom or bond in highlight_atoms or highlight_bonds is missing
        from the image, including when it is present in the molecule but hidden.
    """
    if isinstance(molecule, FrozenMolecule):
        rdmol = molecule.to_rdkit()
    elif isinstance(molecule, str):
        rdmol = Molecule.from_smiles(molecule).to_rdkit()
    else:
        rdmol = molecule

    if color_by_element is None:
        color_by_element = highlight_atoms is None and highlight_bonds is None

    if explicit_hydrogens is None:
        idx_map = {i: i for i in range(rdmol.GetNumAtoms())}
    elif explicit_hydrogens:
        idx_map = {i: i for i in range(rdmol.GetNumAtoms())}
        rdmol = Chem.AddHs(rdmol, explicitOnly=True)
    else:
        idx_map = {
            old: new
            for new, old in enumerate(
                a.GetIdx()
                for a in rdmol.GetAtoms()
                if a.GetAtomicNum() != 1 and a.GetMass() != 1
            )
        }
        rdmol = Chem.RemoveHs(rdmol, updateExplicitCount=True)

    if highlight_atoms is None:
        highlight_atoms = []
        highlight_atom_colors = None
        highlight_atom_radii = None
    elif isinstance(highlight_atoms, dict):
        highlight_atoms = {
            idx_map[i]: AtomHighlight._from_color_or_self(v)
            for i, v in highlight_atoms.items()
            if i in idx_map
        }
        highlight_atom_colors = {
            k: color_to_floats(v.color)
            for k, v in highlight_atoms.items()
            if v.color is not None
        }
        highlight_atom_radii = {
            k: v.radius.m_as(unit.angstrom)
            for k, v in highlight_atoms.items()
            if v.radius is not None
        }
        highlight_atoms = list(highlight_atoms.keys())
    else:
        highlight_atoms = [idx_map[i] for i in highlight_atoms if i in idx_map]
        highlight_atom_colors = None
        highlight_atom_radii = None

    if highlight_bonds is None:
        highlight_bonds = []
        highlight_bond_colors = None
    elif isinstance(highlight_bonds, dict):
        highlight_bond_colors = {
            rdmol.GetBondBetweenAtoms(
                idx_map[i_a], idx_map[i_b]
            ).GetIdx(): color_to_floats(v)
            for (i_a, i_b), v in highlight_bonds.items()
            if i_a in idx_map and i_b in idx_map
        }

        highlight_bonds = list(highlight_bond_colors.keys())
    else:
        highlight_bonds = [
            rdmol.GetBondBetweenAtoms(idx_map[i_a], idx_map[i_b])
            for i_a, i_b in highlight_bonds
            if i_a in idx_map and i_b in idx_map
        ]
        highlight_bond_colors = None

    if bond_notes is not None:
        for (i_a, i_b), note in bond_notes.items():
            if i_a not in idx_map or i_b not in idx_map:
                continue
            rdbond = rdmol.GetBondBetweenAtoms(idx_map[i_a], idx_map[i_b])
            rdbond.SetProp("bondNote", note)

    if atom_notes is not None:
        for i, note in atom_notes.items():
            if i not in idx_map:
                continue
            rdatom = rdmol.GetAtomWithIdx(idx_map[i])
            rdatom.SetProp("atomNote", note)

    Compute2DCoords(rdmol)

    drawer = Draw.MolDraw2DSVG(image_width, image_height)

    draw_options = drawer.drawOptions()
    if not color_by_element:
        draw_options.useBWAtomPalette()

    drawer.DrawMolecule(
        rdmol,
        highlightAtoms=highlight_atoms,
        highlightAtomColors=highlight_atom_colors,
        highlightAtomRadii=highlight_atom_radii,
        highlightBonds=highlight_bonds,
        highlightBondColors=highlight_bond_colors,
    )
    drawer.FinishDrawing()

    svg_contents = drawer.GetDrawingText()

    return SVG(svg_contents)


def depict_charges(
    offmol: Molecule,
    force_field: ForceField,
    image_width: int = DEFAULT_WIDTH,
    image_height: int = DEFAULT_HEIGHT,
) -> SVG:
    """
    Draw the library charges next to each atom
    """

    labels = force_field.label_molecules(offmol.to_topology())

    atom_notes = {}
    for indices, charge_type in labels[0]["LibraryCharges"].items():
        for charge_idx, atom_idx in enumerate(indices, start=1):
            charge = getattr(charge_type, f"charge{charge_idx}")
            atom_notes[atom_idx] = f"{charge.m:.1}"

    return draw_molecule(
        offmol,
        atom_notes=atom_notes,
        image_height=image_height,
        image_width=image_width,
    )


def depict_charge_source(
    offmol: Molecule,
    force_field: ForceField,
    image_width: int = DEFAULT_WIDTH,
    image_height: int = DEFAULT_HEIGHT,
) -> SVG:
    """
    Color each atom according to where it gets its charges from
    """
    labels = force_field.label_molecules(offmol.to_topology())

    highlight_atom_colors = {}
    for indices, charge_type in labels[0]["LibraryCharges"].items():
        # We don't assign IDs, but ff14sb does:
        if isinstance(charge_type.id, str):
            this_color = tuple(AMBER_COLOR)
        elif charge_type.id is None:
            this_color = tuple(NEW_PARAM_COLOR)
        else:
            continue

        for i in indices:
            highlight_atom_colors[i] = this_color

    return draw_molecule(
        offmol,
        highlight_atoms=highlight_atom_colors,
        image_height=image_height,
        image_width=image_width,
    )


def depict_bonds_source(
    offmol: Molecule,
    force_field: ForceField,
    image_width: int = DEFAULT_WIDTH,
    image_height: int = DEFAULT_HEIGHT,
    show_smirks: bool = False,
) -> SVG:
    """
    Color each bond according to where it gets its torsions from
    """

    labels = force_field.label_molecules(offmol.to_topology())
    sage = ForceField("openff-2.0.0.offxml")
    ff14sb = ForceField("ff14sb_off_impropers_0.0.3.offxml")

    highlight_bond_colors = {}
    bond_notes = {}
    for indices, bond_type in labels[0]["Bonds"].items():
        if bond_in_ff(bond_type, ff14sb):
            highlight_bond_colors[indices] = tuple(AMBER_COLOR)
        elif bond_in_ff(bond_type, sage):
            highlight_bond_colors[indices] = tuple(SAGE_COLOR)

        if show_smirks:
            bond_notes[indices] = bond_type.smirks

    return draw_molecule(
        offmol,
        highlight_bonds=highlight_bond_colors,
        bond_notes=bond_notes,
        image_height=image_height,
        image_width=image_width,
    )


def depict_torsions_source(
    offmol: Molecule,
    force_field: ForceField,
    image_width: int = DEFAULT_WIDTH,
    image_height: int = DEFAULT_HEIGHT,
    show_smirks: bool = False,
) -> SVG:
    """
    Color each torsion according to where it gets its torsions from
    """

    labels = force_field.label_molecules(offmol.to_topology())
    sage = ForceField("openff-2.0.0.offxml")
    ff14sb = ForceField("ff14sb_off_impropers_0.0.3.offxml")

    torsions_per_ff = {}
    bond_notes: Dict[Tuple[int, int], str] = {}
    for indices, torsion_type in labels[0]["ProperTorsions"].items():
        if in_ff("ProperTorsions", torsion_type, ff14sb):
            this_color = np.asarray([1, 0])
        elif in_ff("ProperTorsions", torsion_type, sage):
            this_color = np.asarray([0, 1])
        else:
            this_color = np.asarray([0, 0])

        bond_index = tuple(sorted(indices[1:3]))

        color = torsions_per_ff.setdefault(bond_index, np.zeros(2))
        color += this_color

        if show_smirks:
            bond_notes[bond_index] = torsion_type.smirks

    total_available_torsions = np.max([v[0] + v[1] for v in torsions_per_ff.values()])
    highlight_torsion_colors = {
        k: tuple(
            blend_colors(
                AMBER_COLOR,
                SAGE_COLOR,
                t=v[0] / total_available_torsions,
                w=1 - (v[0] + v[1]) / total_available_torsions,
            )
        )
        for k, v in torsions_per_ff.items()
    }

    return draw_molecule(
        offmol,
        image_width=image_width,
        image_height=image_height,
        highlight_bonds=highlight_torsion_colors,
        bond_notes=bond_notes,
        explicit_hydrogens=False,
    )


class OpenFFMoleculeTrajectory(Structure, Trajectory):
    """OpenFF Molecule adaptor.

    Example
    -------
    >>> import nglview as nv
    >>> import mdtraj as md
    >>> mol = Molecule.from_polymer_pdb(pdb_filename)
    >>> t = OpenFFMoleculeTrajectory(mol)
    >>> w = nv.NGLWidget(t)
    >>> w
    """

    def __init__(self, molecule: Molecule, ext: str = "MOL2"):
        if not molecule.conformers:
            raise ValueError(
                "Cannot visualize a molecule without conformers with NGLView"
            )
        self.molecule = molecule
        self.ext = ext.lower()
        self.params = {}
        self.id = str(uuid.uuid4())

    def get_coordinates(self, index: int):
        return self.molecule.conformers[index].m_as(unit.angstrom)

    @property
    def n_frames(self):
        return len(self.molecule.conformers)

    def get_structure_string(self):
        with StringIO() as f:
            self.molecule.to_file(f, file_format=self.ext)
            sdf_str = f.getvalue()
        return sdf_str


class OpenFFTopologyStructure(Structure):
    """OpenFF Molecule adaptor.

    Example
    -------
    >>> import nglview as nv
    >>> import mdtraj as md
    >>> mol = Molecule.from_polymer_pdb(pdb_filename)
    >>> t = OpenFFMoleculeTrajectory(mol)
    >>> w = nv.NGLWidget(t)
    >>> w
    """

    def __init__(
        self,
        topology: Topology,
        ext: str,
    ):
        self.topology = topology
        self.ext = ext.lower()
        self.params = {}
        self.id = str(uuid.uuid4())

    def get_structure_string(self):
        with StringIO() as f:
            self.topology.to_file(f, file_format=self.ext)
            sdf_str = f.getvalue()
        return sdf_str


def visualize(obj: Union[Molecule, Topology], ext=None, representations=None):
    """Visualize a topology or molecule with nglview"""
    if isinstance(obj, Molecule):
        if representations is None:
            representations = MOLECULE_DEFAULT_REPS
        if ext is None:
            ext = "MOL2"
        return nglview.NGLWidget(
            OpenFFMoleculeTrajectory(obj, ext=ext), representations=representations
        )
    elif isinstance(obj, Topology):
        if ext is None:
            ext = "PDB"
        return nglview.NGLWidget(
            OpenFFTopologyStructure(obj, ext=ext), representations=representations
        )
    else:
        if ext is not None:
            raise ValueError(
                "ext parameter only supported for topologies and molecules"
            )
        return nglview.NGLWidget(structure=obj, representations=representations)
