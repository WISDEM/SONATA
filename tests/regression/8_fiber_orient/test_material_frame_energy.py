"""Regression guard: material-frame (ply-local) stress/strain recovery must
stay energy-consistent with the global (section) frame recovery.

This closes a gap in the existing 6_box_beam regression suite: those tests
use isotropic steel, so the material frame and the global frame always
coincide and can never expose a frame-mismatch bug. Here we use a genuinely
off-axis (rotated) anisotropic ply -- the exact scenario in which `stressM`
was found to be silently paired downstream with the *global* (unrotated)
strain instead of a material-frame strain, corrupting any energy-based
post-processing (e.g. `strain_energy_eval`, modal-strain-energy damping)
for that material while leaving K/M/frequencies unaffected.
"""

import os

import numpy as np
import numpy.testing as npt

from SONATA.classBlade import Blade
from SONATA.utl.beam_struct_eval import beam_struct_eval, strain_energy_eval

run_dir = os.path.dirname(os.path.realpath(__file__))
EXAMPLE_YAML = os.path.join(run_dir, "..", "..", "..", "examples",
                           "8_fiber_orient", "8_box_beam.yaml")


def _run_single_cell_off_axis(loads_dict, fiber_deg=15.0):
    """Single-cell custom mesh with an off-axis CarbonUD ply.

    Mirrors the `8_sonata_orientation.py` example's custom-mesh shortcut so
    the test runs fast (no section meshing) while still exercising a
    nonzero, non-trivial ply rotation through the real b3_secfem /
    beam_struct_eval / recovery pipeline.
    """
    job_name = "Box_Beam"

    flags_dict = {
        "flag_wt_ontology": True,
        "flag_ref_axes_wt": True,
        "attribute_str": "MatID",
        "flag_plotDisplacement": False,
        "flag_plotTheta11": False,
        "flag_wf": False,
        "flag_lft": False,
        "flag_topo": False,
        "mesh_resolution": 400,
        "flag_recovery": True,
        "c2_axis": False,
    }

    # Two stations: strain_energy_eval needs >= 2 points for its trapezoid
    # span integration.
    radial_stations = [0.0, 1.0]
    job = Blade(name=job_name, filename=EXAMPLE_YAML, flags=flags_dict,
                stations=radial_stations)

    nodes = np.array([[-1.0, 0.5, 0.0],
                      [1.0, 0.5, 0.0],
                      [1.0, -0.5, 0.0],
                      [-1.0, -0.5, 0.0]])
    cells = np.array([[1, 2, 3, 0]])
    MatID = np.ones(1)
    job.blade_custom_mesh(nodes, cells, MatID, theta_3=fiber_deg)

    flags_dict["flag_csv_export"] = False
    flags_dict["flag_DeamDyn_def_transform"] = True
    flags_dict["flag_write_BeamDyn"] = True
    flags_dict["flag_write_BeamDyn_unit_convert"] = ""
    flags_dict["flag_output_zero_twist"] = False

    mu = np.zeros(6)
    beam_struct_eval(job_name, flags_dict, loads_dict, radial_stations, job,
                     run_dir, "8_box_beam.yaml", mu)

    return job


def _energy_density(stress, strain):
    """0.5*sigma:eps (Voigt, engineering shear) -- twice the strain energy
    density, but the factor of two is irrelevant for the invariance checks
    below.
    """
    return (stress.sigma11 * strain.epsilon11
            + stress.sigma22 * strain.epsilon22
            + stress.sigma33 * strain.epsilon33
            + stress.sigma23 * strain.gamma23
            + stress.sigma13 * strain.gamma13
            + stress.sigma12 * strain.gamma12)


def test_material_frame_energy_matches_global_frame_under_torque():
    """Energy is frame-invariant: sigma.eps must match whether evaluated in
    the global or the material (ply) frame, for every cell.
    """
    loads_dict = {"Forces": [0.0, 0.0, 0.0], "Moments": [1.0e3, 0.0, 0.0]}
    job = _run_single_cell_off_axis(loads_dict)

    checked = 0
    for _x, cs in job.sections:
        for c in cs.mesh:
            energy_global = _energy_density(c.stress, c.strain)
            energy_mat = _energy_density(c.stressM, c.strainM)
            npt.assert_allclose(
                energy_mat, energy_global, rtol=1e-6,
                atol=1e-6 * max(abs(energy_global), 1.0),
                err_msg="Material-frame energy is not frame-invariant.",
            )
            checked += 1
    assert checked > 0


def test_strainM_differs_from_global_strain_when_rotated():
    """Guard against silently re-aliasing `strainM = strain` (the actual
    historical bug): for a genuinely off-axis ply the two must differ.
    """
    loads_dict = {"Forces": [0.0, 0.0, 0.0], "Moments": [1.0e3, 0.0, 0.0]}
    job = _run_single_cell_off_axis(loads_dict)

    max_diff = 0.0
    max_scale = 1e-30
    for _x, cs in job.sections:
        for c in cs.mesh:
            for comp in ("epsilon11", "epsilon22", "gamma12", "gamma13"):
                max_diff = max(max_diff, abs(getattr(c.strainM, comp) - getattr(c.strain, comp)))
                max_scale = max(max_scale, abs(getattr(c.strain, comp)))
    assert max_diff > 1e-3 * max_scale


def test_strain_energy_eval_is_non_negative_under_torque():
    """`strain_energy_eval` integrates the *material-frame* stress/strain
    (per its own docstring). A frame mismatch shows up as catastrophic
    cancellation or a spurious sign here even though the section is under
    a purely elastic, energy-positive deformation.
    """
    loads_dict = {"Forces": [0.0, 0.0, 0.0], "Moments": [1.0e3, 0.0, 0.0]}
    job = _run_single_cell_off_axis(loads_dict)

    total_energy, directional_energy, _strain_all, _energy_all = \
        strain_energy_eval(job)

    assert total_energy > 0.0
    assert np.all(np.asarray(directional_energy) >= -1e-6 * abs(total_energy))
