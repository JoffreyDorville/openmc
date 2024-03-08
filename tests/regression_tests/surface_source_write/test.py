"""Test the 'surface_source_write' setting.

Results
-------

All results are generated using only 1 MPI process.

All results are generated using 1 thread except for "test_consistency_low_realization_number".
This specific test verifies that when the number of realization (i.e., point being candidate
to be stored) is higher than the capacity, results are repoducible even with multiple
threads (i.e., there is no potential thread competition that would produce different 
results).

All results are generated using the history-based mode except cases 21 to 23.

All results are visually verified using the '_visualize.py' script in the regression test folder.

OpenMC models
-------------

Three OpenMC models with CSG geometries only are used to cover the transmission, vacuum and
reflective Boundary Conditions (BC):

- model_1: complete model with a cylindrical core in 2 boxes (vacuum and transmission BC),
- model_2: simplified model with a cylindrical core in 1 box (vacuum BC),
- model_3: simplified model with a cylindrical core in 1 box (reflective BC).

Two models including DAGMC geometries are also used, based on the mesh file 'dagmc.h5m'
available from tests/regression_tests/dagmc/legacy:

- model_dagmc_1: model adapted from tests/regression_tests/dagmc/legacy,
- model_dagmc_2: model_dagmc_1 contained in two boxes to introduce multiple level of coordinates.

Test cases
----------

Test cases for CSG geometries only:

========  =======  =========  =========================  =====  ===================================
Folder    Model    Surface    Cell                       BC*    Expected particles
========  =======  =========  =========================  =====  ===================================
case-1    model_1  No         No                         T+V    Particles crossing any surface in
                                                                the model
case-2    model_1  1          No                         T      Particles crossing this surface
                                                                only
case-3    model_1  Multiple   No                         T      Particles crossing the declared
                                                                surfaces
case-4    model_1  Multiple   cell (lower universe)      T      Particles crossing the declared
                                                                surfaces that come from or are
                                                                coming to the cell
case-5    model_1  Multiple   cell (root universe)       T      Particles crossing the declared
                                                                surfaces that come from or are
                                                                coming to the cell
case-6    model_1  No         cell (lower universe)      T      Particles crossing any surface that
                                                                come from or are coming to the cell
case-7    model_1  No         cell (root universe)       T      Particles crossing any surface that
                                                                come from or are coming to the cell
case-8    model_1  No         cellfrom (lower universe)  T      Particles crossing any surface that
                                                                come from the cell
case-9    model_1  No         cellto (lower universe)    T      Particles crossing any surface that
                                                                are coming to the cell
case-10   model_1  No         cellfrom (root universe)   T      Particles crossing any surface that
                                                                come from the cell
case-11   model_1  No         cellto (root universe)     T      Particles crossing any surface that
                                                                are coming to the cell
case-12   model_2  Multiple   No                         V      Particles crossing the declared
                                                                surfaces
case-13   model_2  Multiple   cell (root universe)       V      Particles crossing any surface that
                                                                come from or are coming to the cell
case-14   model_2  Multiple   cellfrom (root universe)   V      Particles crossing any surface that
                                                                are coming to the cell
case-15   model_2  Multiple   cellto (root universe)     V      None
case-16   model_3  Multiple   No                         R      Particles crossing the declared
                                                                surfaces
case-17   model_3  Multiple   cell (root universe)       R      None
case-18   model_3  Multiple   cellfrom (root universe)   R      None
case-19   model_3  Multiple   cellto (root universe)     R      None
========  =======  =========  =========================  =====  ===================================

*: BC stands for Boundary Conditions, T for Transmission, R for Reflective, and V for Vacuum.

An additional case, called 'case-20', is used to check that the results are comparable when
the number of threads is set to 2 if the number of realization is lower than the capacity.

Cases 21 to 23 are the event-based cases corresponding to the history-based cases 4, 7 and 13,
respectively.

Test cases including DAGMC geometries:

========  =============  =========  =====================  =====  ===================================
Folder    Model          Surface    Cell                   BC*    Expected particles
========  =============  =========  =====================  =====  ===================================
case-24   model_dagmc_1  No         No                     T+V    Particles crossing any surface in
                                                                  the model
case-25   model_dagmc_1  1          No                     T      Particles crossing this surface
                                                                  only
case-26   model_dagmc_1  No         cell                   T      Particles crossing any surface that
                                                                  come from or are coming to the cell
case-27   model_dagmc_1  1          cell                   T      Particles crossing the declared
                                                                  surface that come from or are
                                                                  coming to the cell
case-28   model_dagmc_1  No         cellfrom               T      Particles crossing any surface that
                                                                  come from the cell
case-29   model_dagmc_1  No         cellto                 T      Particles crossing any surface that
                                                                  are coming to the cell
case-30   model_dagmc_2  Multiple   cell (lower universe)  T      Particles crossing the declared
                                                                  surfaces that come from or are
                                                                  coming to the cell
case-31   model_dagmc_2  Multiple   cell (root universe)   T      Particles crossing the declared
                                                                  surfaces that come from or are
                                                                  coming to the cell
========  =============  =========  =====================  =====  ===================================

*: BC stands for Boundary Conditions, T for Transmission, and V for Vacuum.

Notes:

- The test cases list is non-exhaustive compared to the number of possible combinations.
  Test cases have been selected based on use and internal code logic.
- Cases 8 to 11 are testing that the feature still works even if the level of coordinates
  before and after crossing a surface is different,
- Tests on boundary conditions are not performed on DAGMC models as the logic is shared
  with CSG only models,
- Cases that should return an error are tested in the 'test_exceptions' unit test
  from 'test_surf_source_write.py'.

TODO:

- Test with a lattice,
- Test with periodic boundary conditions.

"""

import os
import shutil
from pathlib import Path

import openmc
import openmc.lib

import h5py
import numpy as np
import pytest

from tests.testing_harness import PyAPITestHarness
from tests.regression_tests import config


@pytest.fixture(scope="function")
def single_thread(monkeypatch):
    """Set the number of OMP threads to 1 for the test."""
    monkeypatch.setenv("OMP_NUM_THREADS", "1")


@pytest.fixture(scope="function")
def two_threads(monkeypatch):
    """Set the number of OMP threads to 2 for the test."""
    monkeypatch.setenv("OMP_NUM_THREADS", "2")


@pytest.fixture(scope="function")
def single_process(monkeypatch):
    """Set the number of MPI process to 1 for the test."""
    monkeypatch.setitem(config, "mpi_np", "1")


@pytest.fixture(scope="module")
def model_1():
    """Cylindrical core contained in a first box which is contained in a larger box.
    A lower universe is used to describe the interior of the first box which
    contains the core and its surrounding space.

    """
    openmc.reset_auto_ids()
    model = openmc.model.Model()

    # =============================================================================
    # Materials
    # =============================================================================

    fuel = openmc.Material()
    fuel.add_element("U", 1.0, enrichment=5.0)
    fuel.add_nuclide("O16", 2.0)
    fuel.set_density("g/cm3", 11.0)

    water = openmc.Material()
    water.add_nuclide("H1", 2.0)
    water.add_nuclide("O16", 1.0)
    water.set_density("g/cm3", 1.0)

    # =============================================================================
    # Geometry
    # =============================================================================

    # -----------------------------------------------------------------------------
    # Cylindrical core
    # -----------------------------------------------------------------------------

    # Parameters
    core_radius = 2.0
    core_height = 4.0

    # Surfaces
    core_cylinder = openmc.ZCylinder(r=core_radius)
    core_lower_plane = openmc.ZPlane(z0=-core_height / 2.0)
    core_upper_plane = openmc.ZPlane(z0=core_height / 2.0)

    # Region
    core_region = -core_cylinder & +core_lower_plane & -core_upper_plane

    # Cells
    core = openmc.Cell(fill=fuel, region=core_region)
    outside_core_region = +core_cylinder | -core_lower_plane | +core_upper_plane
    outside_core = openmc.Cell(fill=water, region=outside_core_region)

    # Universe
    inside_box1_universe = openmc.Universe(cells=[core, outside_core])

    # -----------------------------------------------------------------------------
    # Box 1
    # -----------------------------------------------------------------------------

    # Parameters
    box1_size = 6.0

    # Surfaces
    box1_lower_plane = openmc.ZPlane(z0=-box1_size / 2.0)
    box1_upper_plane = openmc.ZPlane(z0=box1_size / 2.0)
    box1_left_plane = openmc.XPlane(x0=-box1_size / 2.0)
    box1_right_plane = openmc.XPlane(x0=box1_size / 2.0)
    box1_rear_plane = openmc.YPlane(y0=-box1_size / 2.0)
    box1_front_plane = openmc.YPlane(y0=box1_size / 2.0)

    # Region
    box1_region = (
        +box1_lower_plane
        & -box1_upper_plane
        & +box1_left_plane
        & -box1_right_plane
        & +box1_rear_plane
        & -box1_front_plane
    )

    # Cell
    box1 = openmc.Cell(fill=inside_box1_universe, region=box1_region)

    # -----------------------------------------------------------------------------
    # Box 2
    # -----------------------------------------------------------------------------

    # Parameters
    box2_size = 8

    # Surfaces
    box2_lower_plane = openmc.ZPlane(z0=-box2_size / 2.0, boundary_type="vacuum")
    box2_upper_plane = openmc.ZPlane(z0=box2_size / 2.0, boundary_type="vacuum")
    box2_left_plane = openmc.XPlane(x0=-box2_size / 2.0, boundary_type="vacuum")
    box2_right_plane = openmc.XPlane(x0=box2_size / 2.0, boundary_type="vacuum")
    box2_rear_plane = openmc.YPlane(y0=-box2_size / 2.0, boundary_type="vacuum")
    box2_front_plane = openmc.YPlane(y0=box2_size / 2.0, boundary_type="vacuum")

    # Region
    inside_box2 = (
        +box2_lower_plane
        & -box2_upper_plane
        & +box2_left_plane
        & -box2_right_plane
        & +box2_rear_plane
        & -box2_front_plane
    )
    outside_box1 = (
        -box1_lower_plane
        | +box1_upper_plane
        | -box1_left_plane
        | +box1_right_plane
        | -box1_rear_plane
        | +box1_front_plane
    )

    box2_region = inside_box2 & outside_box1

    # Cell
    box2 = openmc.Cell(fill=water, region=box2_region)

    # Root universe
    root = openmc.Universe(cells=[box1, box2])

    # Register geometry
    model.geometry = openmc.Geometry(root)

    # =============================================================================
    # Settings
    # =============================================================================

    model.settings = openmc.Settings()
    model.settings.run_mode = "eigenvalue"
    model.settings.particles = 1000
    model.settings.batches = 5
    model.settings.inactive = 1
    model.settings.seed = 1

    bounds = [
        -core_radius,
        -core_radius,
        -core_height / 2.0,
        core_radius,
        core_radius,
        core_height / 2.0,
    ]
    distribution = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
    model.settings.source = openmc.source.IndependentSource(space=distribution)

    return model


@pytest.fixture
def model_2():
    """Cylindrical core contained in a box.
    A lower universe is used to describe the interior of the box which
    contains the core and its surrounding space.

    The box is defined with vacuum boundary conditions.

    """
    openmc.reset_auto_ids()
    model = openmc.model.Model()

    # =============================================================================
    # Materials
    # =============================================================================

    fuel = openmc.Material()
    fuel.add_element("U", 1.0, enrichment=5.0)
    fuel.add_nuclide("O16", 2.0)
    fuel.set_density("g/cm3", 11.0)

    water = openmc.Material()
    water.add_nuclide("H1", 2.0)
    water.add_nuclide("O16", 1.0)
    water.set_density("g/cm3", 1.0)

    # =============================================================================
    # Geometry
    # =============================================================================

    # -----------------------------------------------------------------------------
    # Cylindrical core
    # -----------------------------------------------------------------------------

    # Parameters
    core_radius = 2.0
    core_height = 4.0

    # Surfaces
    core_cylinder = openmc.ZCylinder(r=core_radius)
    core_lower_plane = openmc.ZPlane(z0=-core_height / 2.0)
    core_upper_plane = openmc.ZPlane(z0=core_height / 2.0)

    # Region
    core_region = -core_cylinder & +core_lower_plane & -core_upper_plane

    # Cells
    core = openmc.Cell(fill=fuel, region=core_region)
    outside_core_region = +core_cylinder | -core_lower_plane | +core_upper_plane
    outside_core = openmc.Cell(fill=water, region=outside_core_region)

    # Universe
    inside_box1_universe = openmc.Universe(cells=[core, outside_core])

    # -----------------------------------------------------------------------------
    # Box 1
    # -----------------------------------------------------------------------------

    # Parameters
    box1_size = 6.0

    # Surfaces
    box1_lower_plane = openmc.ZPlane(z0=-box1_size / 2.0, boundary_type="vacuum")
    box1_upper_plane = openmc.ZPlane(z0=box1_size / 2.0, boundary_type="vacuum")
    box1_left_plane = openmc.XPlane(x0=-box1_size / 2.0, boundary_type="vacuum")
    box1_right_plane = openmc.XPlane(x0=box1_size / 2.0, boundary_type="vacuum")
    box1_rear_plane = openmc.YPlane(y0=-box1_size / 2.0, boundary_type="vacuum")
    box1_front_plane = openmc.YPlane(y0=box1_size / 2.0, boundary_type="vacuum")

    # Region
    box1_region = (
        +box1_lower_plane
        & -box1_upper_plane
        & +box1_left_plane
        & -box1_right_plane
        & +box1_rear_plane
        & -box1_front_plane
    )

    # Cell
    box1 = openmc.Cell(fill=inside_box1_universe, region=box1_region)

    # Root universe
    root = openmc.Universe(cells=[box1])

    # Register geometry
    model.geometry = openmc.Geometry(root)

    # =============================================================================
    # Settings
    # =============================================================================

    model.settings = openmc.Settings()
    model.settings.run_mode = "eigenvalue"
    model.settings.particles = 1000
    model.settings.batches = 5
    model.settings.inactive = 1
    model.settings.seed = 1

    bounds = [
        -core_radius,
        -core_radius,
        -core_height / 2.0,
        core_radius,
        core_radius,
        core_height / 2.0,
    ]
    distribution = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
    model.settings.source = openmc.source.IndependentSource(space=distribution)

    return model


@pytest.fixture
def model_3():
    """Cylindrical core contained in a box.
    A lower universe is used to describe the interior of the box which
    contains the core and its surrounding space.

    The box is defined with reflective boundary conditions.

    """
    openmc.reset_auto_ids()
    model = openmc.model.Model()

    # =============================================================================
    # Materials
    # =============================================================================

    fuel = openmc.Material()
    fuel.add_element("U", 1.0, enrichment=5.0)
    fuel.add_nuclide("O16", 2.0)
    fuel.set_density("g/cm3", 11.0)

    water = openmc.Material()
    water.add_nuclide("H1", 2.0)
    water.add_nuclide("O16", 1.0)
    water.set_density("g/cm3", 1.0)

    # =============================================================================
    # Geometry
    # =============================================================================

    # -----------------------------------------------------------------------------
    # Cylindrical core
    # -----------------------------------------------------------------------------

    # Parameters
    core_radius = 2.0
    core_height = 4.0

    # Surfaces
    core_cylinder = openmc.ZCylinder(r=core_radius)
    core_lower_plane = openmc.ZPlane(z0=-core_height / 2.0)
    core_upper_plane = openmc.ZPlane(z0=core_height / 2.0)

    # Region
    core_region = -core_cylinder & +core_lower_plane & -core_upper_plane

    # Cells
    core = openmc.Cell(fill=fuel, region=core_region)
    outside_core_region = +core_cylinder | -core_lower_plane | +core_upper_plane
    outside_core = openmc.Cell(fill=water, region=outside_core_region)

    # Universe
    inside_box1_universe = openmc.Universe(cells=[core, outside_core])

    # -----------------------------------------------------------------------------
    # Box 1
    # -----------------------------------------------------------------------------

    # Parameters
    box1_size = 6.0

    # Surfaces
    box1_lower_plane = openmc.ZPlane(z0=-box1_size / 2.0, boundary_type="reflective")
    box1_upper_plane = openmc.ZPlane(z0=box1_size / 2.0, boundary_type="reflective")
    box1_left_plane = openmc.XPlane(x0=-box1_size / 2.0, boundary_type="reflective")
    box1_right_plane = openmc.XPlane(x0=box1_size / 2.0, boundary_type="reflective")
    box1_rear_plane = openmc.YPlane(y0=-box1_size / 2.0, boundary_type="reflective")
    box1_front_plane = openmc.YPlane(y0=box1_size / 2.0, boundary_type="reflective")

    # Region
    box1_region = (
        +box1_lower_plane
        & -box1_upper_plane
        & +box1_left_plane
        & -box1_right_plane
        & +box1_rear_plane
        & -box1_front_plane
    )

    # Cell
    box1 = openmc.Cell(fill=inside_box1_universe, region=box1_region)

    # Root universe
    root = openmc.Universe(cells=[box1])

    # Register geometry
    model.geometry = openmc.Geometry(root)

    # =============================================================================
    # Settings
    # =============================================================================

    model.settings = openmc.Settings()
    model.settings.run_mode = "eigenvalue"
    model.settings.particles = 1000
    model.settings.batches = 5
    model.settings.inactive = 1
    model.settings.seed = 1

    bounds = [
        -core_radius,
        -core_radius,
        -core_height / 2.0,
        core_radius,
        core_radius,
        core_height / 2.0,
    ]
    distribution = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
    model.settings.source = openmc.source.IndependentSource(space=distribution)

    return model


def return_surface_source_data(filepath):
    """Read a surface source file and return a sorted array composed
    of flatten arrays of source data for each surface source point.

    TODO:

    - use read_source_file from source.py instead. Or a dedicated function
      to produce sorted list of source points for a given file.

    Parameters
    ----------
    filepath : str
        Path to the surface source file

    Returns
    -------
    data : np.array
        Sorted array composed of flatten arrays of source data for
        each surface source point

    """
    data = []
    keys = []

    # Read source file
    with h5py.File(filepath, "r") as f:
        for point in f["source_bank"]:
            r = point["r"]
            u = point["u"]
            e = point["E"]
            time = point["time"]
            wgt = point["wgt"]
            delayed_group = point["delayed_group"]
            surf_id = point["surf_id"]
            particle = point["particle"]

            key = (
                f"{r[0]:.10e} {r[1]:.10e} {r[2]:.10e} {u[0]:.10e} {u[1]:.10e} {u[2]:.10e}"
                f"{e:.10e} {time:.10e} {wgt:.10e} {delayed_group} {surf_id} {particle}"
            )

            keys.append(key)

            values = [*r, *u, e, time, wgt, delayed_group, surf_id, particle]
            assert len(values) == 12
            data.append(values)

    data = np.array(data)
    keys = np.array(keys)
    sorted_idx = np.argsort(keys)

    data = data[sorted_idx]

    return data


class SurfaceSourceWriteTestHarness(PyAPITestHarness):
    def __init__(self, statepoint_name, model=None, inputs_true=None, workdir=None):
        super().__init__(statepoint_name, model, inputs_true)
        self.workdir = workdir

    def _test_output_created(self):
        """Make sure surface_source.h5 has also been created."""
        super()._test_output_created()
        if self._model.settings.surf_source_write:
            assert os.path.exists(
                "surface_source.h5"
            ), "Surface source file has not been created."

    def _compare_output(self):
        """Compare surface_source.h5 files."""
        if self._model.settings.surf_source_write:
            source_true = return_surface_source_data("surface_source_true.h5")
            source_test = return_surface_source_data("surface_source.h5")
            np.testing.assert_allclose(source_true, source_test, rtol=1e-07)

    def main(self):
        """Accept commandline arguments and either run or update tests."""
        if config["build_inputs"]:
            self.build_inputs()
        elif config["update"]:
            self.update_results()
        else:
            self.execute_test()

    def build_inputs(self):
        """Build inputs."""
        base_dir = os.getcwd()
        try:
            os.chdir(self.workdir)
            self._build_inputs()
        finally:
            os.chdir(base_dir)

    def execute_test(self):
        """Build inputs, run OpenMC, and verify correct results."""
        base_dir = os.getcwd()
        try:
            os.chdir(self.workdir)
            self._build_inputs()
            inputs = self._get_inputs()
            self._write_inputs(inputs)
            self._compare_inputs()
            self._run_openmc()
            self._test_output_created()
            self._compare_output()
            results = self._get_results()
            self._write_results(results)
            self._compare_results()
        finally:
            self._cleanup()
            os.chdir(base_dir)

    def update_results(self):
        """Update results_true.dat and inputs_true.dat"""
        base_dir = os.getcwd()
        try:
            os.chdir(self.workdir)
            self._build_inputs()
            inputs = self._get_inputs()
            self._write_inputs(inputs)
            self._overwrite_inputs()
            self._run_openmc()
            self._test_output_created()
            results = self._get_results()
            self._write_results(results)
            self._overwrite_results()
        finally:
            self._cleanup()
            os.chdir(base_dir)

    def _overwrite_results(self):
        """Also add the 'surface_source.h5' file during overwriting."""
        super()._overwrite_results()
        if os.path.exists("surface_source.h5"):
            shutil.copyfile("surface_source.h5", "surface_source_true.h5")

    def _cleanup(self):
        """Also remove the 'surface_source.h5' file while cleaning."""
        super()._cleanup()
        fs = "surface_source.h5"
        if os.path.exists(fs):
            os.remove(fs)


@pytest.mark.skipif(config["event"] is True, reason="Results from history-based mode.")
@pytest.mark.parametrize(
    "folder, model_name, parameter",
    [
        ("case-1", "model_1", {"max_particles": 3000}),
        ("case-2", "model_1", {"max_particles": 3000, "surface_ids": [4]}),
        (
            "case-3",
            "model_1",
            {"max_particles": 3000, "surface_ids": [4, 5, 6, 7, 8, 9]},
        ),
        (
            "case-4",
            "model_1",
            {"max_particles": 3000, "surface_ids": [4, 5, 6, 7, 8, 9], "cell": 2},
        ),
        (
            "case-5",
            "model_1",
            {"max_particles": 3000, "surface_ids": [4, 5, 6, 7, 8, 9], "cell": 3},
        ),
        ("case-6", "model_1", {"max_particles": 3000, "cell": 2}),
        ("case-7", "model_1", {"max_particles": 3000, "cell": 3}),
        ("case-8", "model_1", {"max_particles": 3000, "cellfrom": 2}),
        ("case-9", "model_1", {"max_particles": 3000, "cellto": 2}),
        ("case-10", "model_1", {"max_particles": 3000, "cellfrom": 3}),
        ("case-11", "model_1", {"max_particles": 3000, "cellto": 3}),
        (
            "case-12",
            "model_2",
            {"max_particles": 3000, "surface_ids": [4, 5, 6, 7, 8, 9]},
        ),
        (
            "case-13",
            "model_2",
            {"max_particles": 3000, "surface_ids": [4, 5, 6, 7, 8, 9], "cell": 3},
        ),
        (
            "case-14",
            "model_2",
            {"max_particles": 3000, "surface_ids": [4, 5, 6, 7, 8, 9], "cellfrom": 3},
        ),
        (
            "case-15",
            "model_2",
            {"max_particles": 3000, "surface_ids": [4, 5, 6, 7, 8, 9], "cellto": 3},
        ),
        (
            "case-16",
            "model_3",
            {"max_particles": 3000, "surface_ids": [4, 5, 6, 7, 8, 9]},
        ),
        (
            "case-17",
            "model_3",
            {"max_particles": 3000, "surface_ids": [4, 5, 6, 7, 8, 9], "cell": 3},
        ),
        (
            "case-18",
            "model_3",
            {"max_particles": 3000, "surface_ids": [4, 5, 6, 7, 8, 9], "cellfrom": 3},
        ),
        (
            "case-19",
            "model_3",
            {"max_particles": 3000, "surface_ids": [4, 5, 6, 7, 8, 9], "cellto": 3},
        ),
    ],
)
def test_surface_source_cell_history_based(
    folder, model_name, parameter, single_thread, single_process, request
):
    """Test on a generic model with vacuum and transmission boundary conditions."""
    assert os.environ["OMP_NUM_THREADS"] == "1"
    assert config["mpi_np"] == "1"
    model = request.getfixturevalue(model_name)
    model.settings.surf_source_write = parameter
    harness = SurfaceSourceWriteTestHarness(
        "statepoint.5.h5", model=model, workdir=folder
    )
    harness.main()


@pytest.mark.skipif(config["event"] is True, reason="Results from history-based mode.")
def test_consistency_low_realization_number(model_1, two_threads, single_process):
    """The objective is to test that the results produced, in a case where
    the number of potential realization (particle storage) is low
    compared to the capacity of storage, are still consistent.

    This configuration ensures that the competition between threads does not
    occur and that the content of the source file created can be compared.

    """
    assert os.environ["OMP_NUM_THREADS"] == "2"
    assert config["mpi_np"] == "1"
    model_1.settings.surf_source_write = {
        "max_particles": 200,
        "surface_ids": [2],
        "cellfrom": 2,
    }
    harness = SurfaceSourceWriteTestHarness(
        "statepoint.5.h5", model=model_1, workdir="case-20"
    )
    harness.main()


@pytest.mark.skipif(config["event"] is False, reason="Results from event-based mode.")
@pytest.mark.parametrize(
    "folder, model_name, parameter",
    [
        (
            "case-21",
            "model_1",
            {"max_particles": 3000, "surface_ids": [4, 5, 6, 7, 8, 9], "cell": 2},
        ),
        ("case-22", "model_1", {"max_particles": 3000, "cell": 3}),
        (
            "case-23",
            "model_2",
            {"max_particles": 3000, "surface_ids": [4, 5, 6, 7, 8, 9], "cell": 3},
        ),
    ],
)
def test_surface_source_cell_event_based(
    folder, model_name, parameter, single_thread, single_process, request
):
    """Test on event-based results."""
    assert os.environ["OMP_NUM_THREADS"] == "1"
    assert config["mpi_np"] == "1"
    model = request.getfixturevalue(model_name)
    model.settings.surf_source_write = parameter
    harness = SurfaceSourceWriteTestHarness(
        "statepoint.5.h5", model=model, workdir=folder
    )
    harness.main()


@pytest.fixture(scope="module")
def model_dagmc_1():
    """Model based on the mesh file 'dagmc.h5m' available from
    tests/regression_tests/dagmc/legacy.

    """
    openmc.reset_auto_ids()
    model = openmc.Model()

    # =============================================================================
    # Materials
    # =============================================================================

    u235 = openmc.Material(name="no-void fuel")
    u235.add_nuclide("U235", 1.0, "ao")
    u235.set_density("g/cc", 11)
    u235.id = 40

    water = openmc.Material(name="water")
    water.add_nuclide("H1", 2.0, "ao")
    water.add_nuclide("O16", 1.0, "ao")
    water.set_density("g/cc", 1.0)
    water.add_s_alpha_beta("c_H_in_H2O")
    water.id = 41

    materials = openmc.Materials([u235, water])
    model.materials = materials

    # =============================================================================
    # Geometry
    # =============================================================================

    dagmc_univ = openmc.DAGMCUniverse(Path("../dagmc.h5m"))
    model.geometry = openmc.Geometry(dagmc_univ)

    # =============================================================================
    # Settings
    # =============================================================================

    model.settings = openmc.Settings()
    model.settings.run_mode = "eigenvalue"
    model.settings.particles = 100
    model.settings.batches = 5
    model.settings.inactive = 1
    model.settings.seed = 1

    source_box = openmc.stats.Box([-4, -4, -20], [4, 4, 20], only_fissionable=True)
    model.settings.source = openmc.IndependentSource(space=source_box)

    return model


@pytest.fixture(scope="module")
def model_dagmc_2():
    """Model based on the mesh file 'dagmc.h5m' available from
    tests/regression_tests/dagmc/legacy.

    This model corresponds to the model_dagmc_1 contained in two boxes to introduce
    multiple level of coordinates from CSG geometry.

    """
    openmc.reset_auto_ids()
    model = openmc.Model()

    # =============================================================================
    # Materials
    # =============================================================================

    u235 = openmc.Material(name="no-void fuel")
    u235.add_nuclide("U235", 1.0, "ao")
    u235.set_density("g/cc", 11)
    u235.id = 40

    water = openmc.Material(name="water")
    water.add_nuclide("H1", 2.0, "ao")
    water.add_nuclide("O16", 1.0, "ao")
    water.set_density("g/cc", 1.0)
    water.add_s_alpha_beta("c_H_in_H2O")
    water.id = 41

    materials = openmc.Materials([u235, water])
    model.materials = materials

    # =============================================================================
    # Geometry
    # =============================================================================

    dagmc_univ = openmc.DAGMCUniverse(Path("../dagmc.h5m"))

    # -----------------------------------------------------------------------------
    # Box 1
    # -----------------------------------------------------------------------------

    # Parameters
    box1_size = 44

    # Surfaces
    box1_lower_plane = openmc.ZPlane(z0=-box1_size / 2.0, surface_id=101)
    box1_upper_plane = openmc.ZPlane(z0=box1_size / 2.0, surface_id=102)
    box1_left_plane = openmc.XPlane(x0=-box1_size / 2.0, surface_id=103)
    box1_right_plane = openmc.XPlane(x0=box1_size / 2.0, surface_id=104)
    box1_rear_plane = openmc.YPlane(y0=-box1_size / 2.0, surface_id=105)
    box1_front_plane = openmc.YPlane(y0=box1_size / 2.0, surface_id=106)

    # Region
    box1_region = (
        +box1_lower_plane
        & -box1_upper_plane
        & +box1_left_plane
        & -box1_right_plane
        & +box1_rear_plane
        & -box1_front_plane
    )

    # Cell
    box1 = openmc.Cell(fill=dagmc_univ, region=box1_region, cell_id=8)

    # -----------------------------------------------------------------------------
    # Box 2
    # -----------------------------------------------------------------------------

    # Parameters
    box2_size = 48

    # Surfaces
    box2_lower_plane = openmc.ZPlane(
        z0=-box2_size / 2.0, boundary_type="vacuum", surface_id=107
    )
    box2_upper_plane = openmc.ZPlane(
        z0=box2_size / 2.0, boundary_type="vacuum", surface_id=108
    )
    box2_left_plane = openmc.XPlane(
        x0=-box2_size / 2.0, boundary_type="vacuum", surface_id=109
    )
    box2_right_plane = openmc.XPlane(
        x0=box2_size / 2.0, boundary_type="vacuum", surface_id=110
    )
    box2_rear_plane = openmc.YPlane(
        y0=-box2_size / 2.0, boundary_type="vacuum", surface_id=111
    )
    box2_front_plane = openmc.YPlane(
        y0=box2_size / 2.0, boundary_type="vacuum", surface_id=112
    )

    # Region
    inside_box2 = (
        +box2_lower_plane
        & -box2_upper_plane
        & +box2_left_plane
        & -box2_right_plane
        & +box2_rear_plane
        & -box2_front_plane
    )
    outside_box1 = (
        -box1_lower_plane
        | +box1_upper_plane
        | -box1_left_plane
        | +box1_right_plane
        | -box1_rear_plane
        | +box1_front_plane
    )

    box2_region = inside_box2 & outside_box1

    # Cell
    box2 = openmc.Cell(fill=water, region=box2_region, cell_id=9)

    # Root universe
    root = openmc.Universe(cells=[box1, box2])

    # Register geometry
    model.geometry = openmc.Geometry(root)

    # =============================================================================
    # Settings
    # =============================================================================

    model.settings = openmc.Settings()
    model.settings.run_mode = "eigenvalue"
    model.settings.particles = 100
    model.settings.batches = 5
    model.settings.inactive = 1
    model.settings.seed = 1

    source_box = openmc.stats.Box([-4, -4, -20], [4, 4, 20], only_fissionable=True)
    model.settings.source = openmc.IndependentSource(space=source_box)

    return model


@pytest.mark.skipif(
    not openmc.lib._dagmc_enabled(), reason="DAGMC CAD geometry is not enabled."
)
@pytest.mark.skipif(config["event"] is True, reason="Results from history-based mode.")
@pytest.mark.parametrize(
    "folder, model_name, parameter",
    [
        ("case-24", "model_dagmc_1", {"max_particles": 300}),
        ("case-25", "model_dagmc_1", {"max_particles": 300, "surface_ids": [1]}),
        ("case-26", "model_dagmc_1", {"max_particles": 300, "cell": 2}),
        (
            "case-27",
            "model_dagmc_1",
            {"max_particles": 300, "surface_ids": [1], "cell": 2},
        ),
        ("case-28", "model_dagmc_1", {"max_particles": 300, "cellfrom": 2}),
        ("case-29", "model_dagmc_1", {"max_particles": 300, "cellto": 2}),
        (
            "case-30",
            "model_dagmc_2",
            {
                "max_particles": 300,
                "surface_ids": [101, 102, 103, 104, 105, 106],
                "cell": 7,
            },
        ),
        (
            "case-31",
            "model_dagmc_2",
            {
                "max_particles": 300,
                "surface_ids": [101, 102, 103, 104, 105, 106],
                "cell": 8,
            },
        ),
    ],
)
def test_surface_source_cell_dagmc(
    folder, model_name, parameter, single_thread, single_process, request
):
    """Test on a generic model with vacuum and transmission boundary conditions."""
    assert os.environ["OMP_NUM_THREADS"] == "1"
    assert config["mpi_np"] == "1"
    model = request.getfixturevalue(model_name)
    model.settings.surf_source_write = parameter
    harness = SurfaceSourceWriteTestHarness(
        "statepoint.5.h5", model=model, workdir=folder
    )
    harness.main()
