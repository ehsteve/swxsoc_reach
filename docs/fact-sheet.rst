================
REACH Fact Sheet
================

The Responsive Environmental Assessment Commercially Hosted (REACH)
dosimeter network provides global, low-latency monitoring of space weather
hazards. The ``swxsoc_reach`` Python package gets, processes, and analyzes
REACH data.

At a Glance
===========

* **Constellation:** 64 dosimeters on 32 Iridium NEXT spacecraft.
* **Measurements:** Electron and proton populations, total dose, internal
  charging, and single-event effects.
* **Coverage:** High-cadence, global measurements across six sensor types.
* **Cadence:** 1 Hz sampling with less than 20-minute constellation revisit.
* **Users:** Space weather researchers and analysts working with REACH data.

Mission Context
===============

REACH is an Air Force-led demonstration of radiation-sensing payloads across
six Iridium NEXT orbital planes. The network helps characterize radiation
environments that can contribute to spacecraft anomalies and separates natural
space weather effects from other disturbances.

Dosimeter Flavors
=================

Each REACH spacecraft carries two channels selected from six flavors.

.. list-table::
   :header-rows: 1
   :widths: 12 68 20

   * - Flavor
     - Energy threshold
     - Sensors
   * - U
     - >= 5.0 MeV electrons; >= 57 MeV protons
     - 5
   * - V
     - >= 3.4 MeV electrons; >= 47 MeV protons
     - 7
   * - W
     - >= 12 MeV protons
     - 14
   * - X
     - >= 360 keV electrons; >= 12 MeV protons
     - 20
   * - Y
     - >= 1.6 MeV electrons; >= 31 MeV protons
     - 12
   * - Z
     - >= 50 keV electrons; >= 200 keV protons
     - 6

What It Does
============

``swxsoc_reach`` supports the complete REACH data workflow:

* Download and process REACH telemetry.
* Apply calibration and coordinate transformations.
* Analyze data by track and geographic region.
* Create visualizations and work with CDF metadata.

Where to Get Data
=================

Public REACH Level 1C data are available from the
`NASA SPDF REACH archive <https://spdf.gsfc.nasa.gov/pub/data/reach/dosimeter/l1c/all_satellites/prelim/>`_.

REACH archive: https://spdf.gsfc.nasa.gov/pub/data/reach/dosimeter/l1c/all_satellites/prelim/

New data are published daily with an approximate two-day latency. For
authorized historical retrieval from the Unified Data Library (UDL), use the
:ref:`historical-download` command-line workflow. UDL downloads require a
``BASICAUTH`` credential or an AWS Secrets Manager configuration.

Typical Workflow
================

#. Obtain REACH data from NASA SPDF or the authorized UDL service.
#. Load and validate the data with ``swxsoc_reach``.
#. Apply the appropriate calibration and transformations.
#. Analyze, visualize, or export the resulting measurements.

Open a File and Generate a Map
==============================

.. code-block:: python

  from swxsoc_reach.track import REACHTrack
  from swxsoc_reach.util import Flavor

  track = REACHTrack.load("reach_l1c.cdf")
  geomap = track.to_geomap()
  geomap.plot(Flavor.X, statistic="median")

Links
=====

* Source repository: https://github.com/swxsoc/swxsoc_reach/
* Documentation: https://swxsoc-reach.readthedocs.io/en/latest/
* User guide: https://swxsoc-reach.readthedocs.io/en/latest/user-guide/overview.html
* Historical UDL download guide:
  https://swxsoc-reach.readthedocs.io/en/latest/user-guide/historical-download.html
