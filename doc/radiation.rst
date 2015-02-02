Radiation correction on a digital elevation model
=================================================


Introduction
------------
The radiation module can be used to adjust radiation estimates
using a (high resolution) digital elevation model. In turn, these results can be used
by the evaporation module to bettes estimate local evaporation.

In general evaporation amounts are determined for about 90% by
radiation input. Radiation at the earth’s surface is determined by
the potential solar radiation at the edge of the earths atmosphere
and the filtering within the atmosphere. The first component can be
easily determined from equations. The reduction due to clouds etc
can be estimated by incorporation short wave measurements but if
these are not available cloud cover estimates can also be used. By
combining this with a DEM radiation at the earths surface can be
determined including the effects of aspect and shading.


Description
-----------

Adjusted after van Dam 2000

This section gives a short description. Another description can
be found at http://re.jrc.cec.eu.int/pvgis/pv/solres/solres.htm.

Potential solar radiation is the radiation of an unobstructed or
cloudless sky. The magnitude of the potential solar radiation depends
on the position of the sun the solar altitude or solar angleduring the day,
the inclination of the solar rays with the earth’s surface, the amount of
radiation at the outer layer of the earth’s atmosphere, the transmissivity
of the sky and the altitude of the earth’s surface.

Solar declination is the annual fluctuation of the sun between the two
tropics and varies between –23 and +23 degrees latitude. Solar declination is
calculated per Day (Julian day number):


.. math::

    \delta = -23.4 cos(360 (Day + 10 )/365)

The hour angle describes the movement of the sun around the earth in 24 hours,
which  equals 15 degrees longitude per hour (360 :math:`deg` /24h). The hour angle n
is calculated for each Hour (whole hour of the day):

.. math::

    n = ( 15 \dot Hour - 12)


The position or height of the sun above the horizon is called the solar altitude
or solar angle. Solar altitude :math:`\alpha` (deg) is calculated for each location,
determined by the location’s latitude :math:`\psi` (deg), declination and hour angle:


.. math::

    sin(\alpha) = sin(\psi) sin(\delta) + cos(\psi) cos(\delta) cos(n)

Solar azimuth is the angle between the solar rays and the North-South axis of the
earth. Solar azimuth :math:`{\beta}_s` (deg) is calculated by:

.. math::

    &cos({\beta}_s) = (sin(\delta) cos(\psi) - cos(\delta) sin(\psi) cos(n))/cos(\alpha) \\
    &for Hour \le 12: {\beta}_s = {\beta}_s \\
    &for Hour > 12: {\beta}_s = 360 - {\beta}_s

Surface azimuth or aspect :math:`{\beta}_1` (deg) is the orientation of the land
surface or slope to the North-South axis of the sun. Slope :math:`\varphi` (deg) is
the maximum rate of change in elevation.

The angle of incidence is the angle between the perpendicular plane of the
incoming solar rays and the surface on which they are projected, defined by the
aspect and slope of that surface. The angle of incidence :math:`\vartheta` (deg) is
calculated with the solar angle :math:`\alpha` (deg), the slope of the land
surface :math:`\varphi` (deg), the azimuth of the sun :math:`{\beta}_s` (deg) and
azimuth of the land surface  :math:`{\beta}_1` (deg):

.. math::

    cos(\vartheta) = cos(\alpha) sin(\varphi) cos({\beta}_s - {\beta}_1) + sin(\alpha) cos(\varphi)

The second section of the radiation module calculates the potential solar energy. The
amount of solar radiation that reaches the outer atmosphere is decreased by the
travelling distance of the solar rays through the sky to the surface, the transmissivity
of the sky and the cloud factor.

Solar energy at the outer layer of the atmosphere :math:`Sout (Wm^2)` is
calculated by (Kreider & Kreith 1975):

.. math::

    S_{out} = S_c (1 + 0.034 cos(360 Day/365))

where :math:`S_c  (Wm^2)` is the solar constant of 1367 :math:`Wm^2` (Duffie & Beckman 1991).
The solar ‘constant’ is subject to much discussion. Gates (1980) gives a value
of 1360 :math:`Wm^2`. The NASA reports a value of 1353 :math:`Wm^2` (Jansen 1985),
while Duncan et al. (1982) give a value of 1367 :math:`Wm^2`. Monteith and Unsworth (1990)
measured the highest value of 1373 W.m-2. The World Radiation Centre uses a
value of 1367 :math:`Wm^2` (Duffie & Beckman 1991) and this value is also used in this study.

The solar radiation energy that reaches the earth’s surface is decreased due to the
length of the air mass it has to pass through and the transmissivity :math:`\tau`
(% or fraction) of the sky. The radiation flux through a hypothetical plane
normal to the beam  (:math:`S_{nor} Wm^2`) is given by (Gates 1980):

.. math::

    S_{nor} = S_{out} \tau^{Mh}

in which Mh (% or fraction) is the relative path length of the optical air mass at
altitude h (m). Transmissivity (:math:`\tau`) is usually between 0.5 and 0.8, but can be as low
as 0.4 in the tropics (Whitmore et al. 1993), but mostly a value of 0.6 is used
(Gates 1980). To calculate the relative path length of an optical air mass at
altitude h (m), the relative path length of an optical air mass at sea level M0
(% or fraction) is corrected for the atmospheric pressure at altitude h. Mh
(% or fraction) is calculated using (Kreider & Kreith 1975):

.. math::

    Mh = M_0 P_h/P_0

in which :math:`P_h / P_0` (mbar.mbar-1) is an atmospheric pressure correction. The relative path length
of the optical air mass at sea level M0 is obtained by (Kreider & Kreith 1975):

.. math::

    M_0 = \sqrt(1299 + (614 sin(\alpha))^2) - 614 sin(\alpha)

The atmospheric pressure correction :math:`P_h / P_0` is written as (List 1984):


.. math::

    P_h / P_0 = ((288 0.0065h) / 288)^5.256


The incoming radiation normal to the beam Snor must be corrected by the orientation and
slope of the surface, defined by the angle of incidence :math:`\vartheta`, to calculate the incoming radiation
Sdir (:math:`Wm^2`) on the earth’s surface:

.. math::

    S_dir = S_nor cos(\vartheta)

Direct light is scattered in the atmosphere. This daylight scattering or diffuse radiation is
approximately 15% of direct radiation (Gates 1980). A more accurate empirical estimation for diffuse
radiation Sdif (:math:`Wm^2`) in a clear not dust-free sky reads as (Liu and Jordan in Gates 1980):


.. math::

    S_dif = S_out (0.271 - 0.294 \vartheta^{Mh} sin(\alpha)

During daylight when the sun is above the horizon, it is assumed that all cells receive the same amount
of diffuse radiation. Total incoming radiation Sin (:math:`Wm^2`) is the sum of direct and diffuse radiation:

.. math::

    S_in = S_dir + S_dif


Total incoming radiation Sin as calculated with the above is actually a radiation flux for that moment. In the procedure
given above, radiation is calculated per time step. If this amount of radiation is used in a water balance model, the amount
of radiation and therewith the amount of evapotranspiration will be overestimated or under estimated, depending on the
time of the day and the position of the sun.



Most of the work done for the shading is implemented in the pcraster horizontan function.


How to use the maps generated to correct model incoming radiation from models or measurements
=============================================================================================

The paragraph below is adapted from the r.sun grass manual:

The real-sky irradiance/irradiation are calculated from clear-sky raster maps by the
application of a factor parameterizing the attenuation of cloud cover. Examples of explicit
calculations of this parameter can be found in Becker (2001), Kitler and Mikler (1986). However, the cloudiness
observation by a meteorological service routine is usually prone to subjective errors and does
not describe sufficiently the physical nature and dynamic spatial-temporal pattern of different
types of cloud cover. Therefore, a simpler parameter has to be used. The solutions for horizontal
and inclined surfaces are slightly different. For the assessment of global irradiance/irradiation
on a horizontal surface under overcast conditions Gh the clear-sky values Ghc are multiplied by
clear-sky index kc (Beyer et al 1996, Hammer et al 1998, Rigollier et al. 2001):

.. math::

	Gh = Ghc kc

The index kc represents the atmospheric transmission expressed as a ratio between horizontal
global radiation under overcast and clear-sky conditions. For a set of ground meteorological
stations the clear-sky index can be calculated from measured global radiation Ghs and
computed values of clear-sky global radiation Ghc:

.. math::

	kc = Ghs/Ghc

As an alternative the kc can be derived also from other climatologic data
(e.g. cloudiness, cf. Kasten and Czeplak 1980). The raster maps of kc must be
then derived by spatial interpolation. The kc can be calculated directly as a raster map from
short-wave surface irradiance measured by satellites. This method is based on the complementarity
between the planetary albedo recorded by the radiometer and the surface radiant flux
(Cano et al 1986, Beyer et al 1996, Hammer et al 1998).
To compute the overcast global irradiance/irradiation for inclined surfaces, Gi
the diffuse Dh and beam Bh components of overcast global radiation and of the clear-sky index kc
have to be treated separately as follows from the following equations:

.. math::

	Dh = Dhc kdc

	Bh = Bhc kbc

The ratio of diffuse to the global radiation Dh/Gh for clear and overcast skies changes
according to the cloudiness. In Europe the Dh/Gh values are typically in interval 0.3-1.0
(Kasten and Czeplak 1980). The underlying physical processes are quite complicated and computationally
represented only by empirical equations (cf. Scharmer and Greif, 2000, Kasten and Czeplak 1980, Hrvoľ 1991).
However, for many meteorological stations, besides the global horizontal radiation Ghs, the diffuse component
Dhs is either measured or calculated from cloudiness, sunshine or other climatologic data.
The raster map of Dhs/Ghs can be derived from the point values by spatial interpolation.
Consecutively, the raster maps of diffuse and beam components of the clear sky index can be computed:

.. math::

	Dh = Gh Dhs/Ghs

	Bh = Gh – Dh

	kdc = Dh/Dh

	kbc = Bh/Bhc


where subscript s is meant to distinguish data measured on meteorological stations Bhs
nd Dhs from the estimated values Bh, and Dh.


Example
=======


.. raw:: html

        <object width="480" height="385"><param name="movie"
        value="http://youtu.be/CEgsx6MPqEM"></param><param
        name="allowFullScreen" value="true"></param><param
        name="allowscriptaccess" value="always"></param><embed
        src="http://youtu.be/CEgsx6MPqEM"
        type="application/x-shockwave-flash" allowscriptaccess="always"
        allowfullscreen="true" width="480"
        height="385"></embed></object>


Implementation
==============

.. automodule:: e2o_dstools.e2o_radiation
    :members:

