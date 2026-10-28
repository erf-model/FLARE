FLARE: Fire-Line Advancement and Reaction Engine
================================================

FLARE is a wildland fire, smoke and dust model built on the
`Energy Research and Forecasting (ERF) <https://github.com/erf-model/ERF>`__
atmospheric code. It is the fire engine behind ERF's fire-atmosphere coupling
layer (``Source/FireCoupling``). ERF supplies the wind, temperature and humidity
near the ground. FLARE spreads the fire on its own finer surface grid and hands
back the heat, moisture and smoke the fire releases. The same surface-grid
design carries a dust emission model for mining and industrial sites. Everything in ERF (compressible and anelastic dynamics, terrain-fitted
meshes, LES and PBL closures, surface layer, radiation, microphysics and land
surface models) is available to a FLARE run.

This repository was the ``ERF-Fire`` branch of ``hgopalan/ERF`` until October
2026. Its ``development`` branch takes upstream ERF ``development`` merges
regularly, so the atmospheric core, build system and inputs stay the same as
ERF's.

=================  ================
Regression Tests    |firetests|
=================  ================

.. |firetests| image:: https://github.com/erf-model/FLARE/actions/workflows/ci.yml/badge.svg?branch=development

What FLARE does
---------------

**Surface fire spread** (``FireEngine/``). The fire lives on a two-dimensional
grid that refines one level of the atmosphere's mesh by an integer factor. It is
built from pieces that are selected independently:

* seven rate-of-spread models under one key, ``erf.fire.ros_model``: Rothermel
  (the default), the BEHAVE multi-class Rothermel, MacArthur, Cheney-Gould and
  its simple grass fit, the Canadian FBP system and Balbi, plus a per-cell hybrid
  of two of them and a prescribed field;
* two front propagation methods, a cell-based Huygens scheme in the spirit of
  FARSITE and a level-set solver;
* fuel as one of Anderson's 13 standard fuel models everywhere, or as a fuel-map
  raster (ESRI ASCII text or a FARSITE ``.lcp`` landscape file), with firebreaks
  and a non-burnable mask from building footprints and listed fuel codes;
* dead and live fuel moisture driven by the atmosphere;
* ignition from a disc, a polygon or polyline perimeter, or a timed schedule;
* startup acceleration, ember spotting and crown fire;
* suppression during the run: fire lines built at a rate, retardant drops, a
  hold test for lines under a hot enough fire;
* two-way coupling of heat, moisture and smoke to the atmosphere, with the fire
  grid placed on any refinement level (a static one above level 0);
* diagnostics: fireline intensity, flame length, flame temperature and tilt,
  fire plotfiles, a statistics CSV, arrival-time probes and full checkpoint and
  restart of the fire state.

**Wildland-urban interface (WUI).** Buildings are obstacles the fire goes around
and targets the exposure diagnostics track. With structure ignition on, each
building carries a state: unignited, burning or burned out. It ignites when the
heat load, the ember count or the fireline intensity at its walls crosses a
threshold (the thresholds are placeholders awaiting calibration). Once burning
it is itself a heat, ember and radiation source for its neighbours, so
house-to-house spread needs no rule of its own.

**Dust** (``Source/Dust``, compiled in with ``ERF_ENABLE_DUST=ON``). Wind-blown
and activity-driven particulate emission from bare mineral surfaces such as
tailings impoundments, open pits, haul roads and evaporation ponds, with the
surface state (silt, crust, salt crusts, moisture, suppression agents) read
from rasters and from the output tables of the geochemistry code PHREEQC. The
dust is a passive scalar of the atmosphere that settles and deposits. The model
reports the regulatory quantities built on it: PM2.5 and PM10 averages against
the EPA air-quality standards (NAAQS), mine-safety (MSHA) worker exposure,
visibility, respirable silica, short-term exposure limits and the emission
budgets of critical minerals. The fire-dust coupling
(``Source/FireDust``) feeds the fire's burned area, outflow wind and heat flux
into the dust emission.

**Fire weather over real terrain.** A one-equation k RANS closure
(``erf.rans_type = kEqn``) for coarse grids on which the boundary layer is not
resolved (horizontal spacings of tens of metres and up), terrain-following
inflow profiles, and the ERF atmosphere's
terrain-fitted meshes, surface layer and PBL schemes.

**Cases.** Every capability has a case folder with its inputs and a README; the
fire and RANS cases also carry a Python check script:

* ``Exec/CanonicalTests/Fire``: fire behaviour, fire-atmosphere coupling,
  heat-flux diagnostics, a WUI subdivision, analytic verification cases
  (level-set advection, slope without wind, merging fires, junction fire,
  obstacle shadow, the UNR WRF-Fire comparison) and the Marshall and Palisades
  fires on real terrain;
* ``Exec/CanonicalTests/Dust`` and ``Exec/CanonicalTests/Hazard``: the dust
  sources, thresholds, deposition and health diagnostics, and dust over hills
  and pits;
* ``Exec/CanonicalTests/Canonical_RANS``: neutral, stable and convective
  boundary layers and flow over hills with the k closure;
* ``Exec/RegTests``: the short regression decks CTest runs.

``Exec/CanonicalTests/Fire/inputs_fire_master_reference`` and its dust and hazard
siblings list every input key with its default; they are documentation, not
decks to run.

Building and running
--------------------

FLARE builds with CMake exactly like ERF; the ERF
`Getting Started <https://erf.readthedocs.io/en/latest/GettingStarted.html>`_
page gives the system requirements. Clone with the submodules::

   git clone --recursive https://github.com/erf-model/FLARE.git FLARE
   cd FLARE
   cmake -S . -B cmake-build -DERF_ENABLE_MPI=ON
   cmake --build cmake-build -j

The fire engine is built in by default (``ERF_ENABLE_FIRE=ON``); dust is opt-in
(``ERF_ENABLE_DUST=ON``). A deck turns the fire on with ``erf.fire.enable = true``
and selects the engine with ``erf.fire.engine = flare``, the default; the
``erf.fire.*`` keys are listed in the inputs documentation. Each case folder
holds the inputs it needs and a README with the expected behaviour.

Python tools for fuel maps from LANDFIRE, HRRR point winds, input-deck checking
and fire-map animations are in ``Exec/Tools``; each real-terrain case keeps its
own terrain generator next to its deck.

Documentation
-------------

The theory and inputs of the fire, WUI, dust and RANS capabilities are
documented in the Sphinx sources under ``Docs/sphinx_doc`` (``theory/Fire.rst``,
``theory/DustModule.rst``, ``theory/RANS.rst`` and the pages they link to, and
the ``erf.fire`` and ``erf.dust`` tables in ``Inputs.rst``), alongside the ERF
documentation. There is no published site yet; build them locally with::

   pip install -r Docs/sphinx_doc/requirements.txt
   sphinx-build -b html Docs/sphinx_doc cmake-build/docs

Development model
-----------------

Development happens on topic branches off ``development`` with pull requests
back into ``development``. New physics is added as an opt-in option with the
existing behaviour kept as the default, and every change to the fire or dust
model comes with a regression or canonical case. Fixes to the shared
atmospheric core are contributed to upstream ERF and reach this repository
through the periodic upstream merges. The ERF coding conventions in
`CONTRIBUTING.md <CONTRIBUTING.md>`_ apply here as well.

Built on ERF
------------

ERF, the Energy Research and Forecasting code (`erf-model/ERF
<https://github.com/erf-model/ERF>`_), solves the atmospheric equations in
compressible or anelastic form on the `AMReX <https://amrex-codes.github.io/amrex/>`_
block-structured adaptive mesh framework, from laptops to GPU systems. Its
`ERF documentation <https://erf.readthedocs.io/en/latest/>`_ and
`Doxygen pages <https://erf-model.github.io/docs/index.html>`_ describe the
atmosphere FLARE runs in; its Python pre- and post-processing tools are in the
`erftools repository <https://github.com/erf-model/erftools/>`_.

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.8102984.svg
   :target: https://doi.org/10.5281/zenodo.8102984

The development of ERF is funded by the U.S. Department of Energy (DOE). When
using FLARE, please cite the ERF publications listed at that DOI and the ERF
model paper in JAMES (`doi:10.1029/2024MS004884
<https://doi.org/10.1029/2024MS004884>`_).

FLARE is distributed under the ERF license; see `LICENSE.md <LICENSE.md>`_.
