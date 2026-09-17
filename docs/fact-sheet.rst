==========================
Aerospace REACH Fact Sheet
==========================

The Aerospace Responsive Environmental Assessment Commercially Hosted (REACH)
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
six Iridium NEXT orbital planes on 32 spacecraft. Each spacecraft carries two REACH dosimeters
built and managed by the Aerospace Corporation.
The goal of the network is to help characterize radiation
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
     - Number of Sensors in Constellation
   * - U
     - >= 2.15 MeV electrons; >= 51.5 MeV protons
     - 5
   * - V
     - >= 2.20 MeV electrons; >= 43.9 MeV protons
     - 7
   * - W
     - >=1.43 MeV electrons, >= 10.5 MeV protons
     - 14
   * - X
     - >= 0.798 keV electrons; >= 12.1 MeV protons
     - 20
   * - Y
     - >= 2.47 MeV electrons; >= 30.3 MeV protons
     - 12
   * - Z
     - >= 91.6 keV electrons; >= 1.29 MeV protons
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

New data are published daily to NASA SPDF REACH Archive which can be found at 
`https://spdf.gsfc.nasa.gov/pub/data/reach/dosimeter/l1c/all_satellites/prelim/ <https://spdf.gsfc.nasa.gov/pub/data/reach/dosimeter/l1c/all_satellites/prelim/>`_
The data are made available with an approximate two-day latency.
To programmatically search for and download REACH data use the **SPDF REACH Client**
(Fido integration), which provides processed Level 1C CDF data without requiring
credentials (see below).
Historical data can be found at `Zenodo REACH Historical Data <https://zenodo.org/records/7038285>`_.

Downloading Data
================

Download the last few days of REACH data using the SPDF REACH Client:

.. code-block:: python

  from sunpy.net import Fido, attrs as a
  from swxsoc_reach.net import REACHClient, DataType, Vehicle

  results = Fido.search(
    a.Time("2025-05-01", "2025-05-05")
    & a.Source.reach
    & a.Instrument.dosimeter
    & Vehicle.all_satellites
  )
  files = Fido.fetch(results, path="./reach_data/")

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
* Software Documentation: https://swxsoc-reach.readthedocs.io/en/latest/
* REACH Historical Data and Reference Documentation: https://zenodo.org/records/7038285 (DOI 10.5281/zenodo.5988170)
* REACH Constellation Info: https://zenodo.org/records/7038285?preview_file=TOR-2019-02650.pdf
* ISWA Gallery: https://iswa.ccmc.gsfc.nasa.gov/app/info?cygnetId=799&dataId=3621