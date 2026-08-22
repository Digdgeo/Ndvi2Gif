"""
hydroperiod.py - GEE-native Hydroperiod Analysis
=================================================

Computes hydroperiod (flood days per hydrological year) entirely in Google
Earth Engine using the midpoint temporal weighting method.

The hydroperiod represents the number of days a pixel remains flooded during
a hydrological cycle (default: September 1st – August 31st).

Methodology
-----------
The midpoint weighting method distributes 365 days among available scenes:
each scene is assigned a temporal "territory" that extends from the midpoint
with its predecessor to the midpoint with its successor (boundaries at day 0
and day 365 for the first and last scenes). Weighted flood and valid days are
then accumulated per pixel and optionally normalised to a full year.

References
----------
Bustamante, J., Aragonés, D., Afán, I., Luque, C.J., Pérez-Vázquez, A.,
Castellanos, E.M., Díaz-Delgado, R. (2016). Predictive models of floristic
composition based on flooding frequency in Mediterranean wetlands.
Remote Sensing, 8(9), 776.

Author: Diego García Díaz
"""

import ee
import numpy as np
from datetime import date


class HydroperiodAnalyzer:
    """GEE-native hydroperiod computation from satellite water indices.

    Computes hydroperiod entirely in Google Earth Engine using the midpoint
    temporal weighting method. Each scene in the satellite collection receives
    a proportional number of days based on its position within the hydrological
    cycle, then weighted flood and valid days are accumulated per pixel.

    Parameters
    ----------
    ndvi_seasonality : NdviSeasonality
        Configured :class:`NdviSeasonality` instance. The ROI, satellite,
        date range, cloud filtering, and band standardisation are reused
        directly — no duplicate configuration needed.
    hydrological_year_start : tuple of (int, int), optional
        ``(month, day)`` marking the start of the hydrological cycle.
        Default ``(9, 1)`` = September 1st.

    Examples
    --------
    ::

        import ee
        from ndvi2gif import NdviSeasonality, HydroperiodAnalyzer

        ee.Initialize()

        ns = NdviSeasonality(roi, sat='S2', start_year=2023, end_year=2024)
        ha = HydroperiodAnalyzer(ns)

        # Full GEE-native computation
        result = ha.compute_hydroperiod(index='mndwi', threshold=0.1)

        # result is an ee.Image with bands:
        #   hydroperiod, valid_days, normalized,
        #   first_flood_doy, last_flood_doy

        # Export to Drive
        ha.export_to_drive(folder='my_wetland')

        # Same, plus the per-date water masks as a separate uint8 file
        # (0 = dry, 1 = water, 2 = cloud-masked, 255 = no acquisition)
        ha.export_to_drive(folder='my_wetland', include_masks=True)

        # Per-pixel IRT (temporal representativity)
        irt = ha.compute_irt_image()

    Notes
    -----
    Supported water indices (optical sensors):
    ``'ndwi'``, ``'mndwi'``, ``'awei'``, ``'aweinsh'``, ``'wi2015'``

    The per-pixel IRT is computed as Simpson's diversity index across monthly
    periods, which is fully computable per-pixel in GEE. The global IRT uses
    the Gini-based formula identical to the standalone *phydroperiod* library
    (computed client-side from scene dates).
    """

    #: Water indices supported for hydroperiod computation.
    WATER_INDICES = frozenset({'ndwi', 'mndwi', 'awei', 'aweinsh', 'wi2015'})

    def __init__(self, ndvi_seasonality, hydrological_year_start=(9, 1)):
        self.ns = ndvi_seasonality
        self.hyd_month, self.hyd_day = hydrological_year_start
        self._water_masks = None
        self._hyd_year = None
        self._d0 = None
        self._results = None
        self._last_index = None
        self._last_threshold = None

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _hyd_year_bounds(self, hyd_year):
        """Return ``(start_str, end_str)`` for the given hydrological year.

        The end date is the exclusive upper bound (start of the *next* cycle),
        matching GEE's ``filterDate`` convention.
        """
        start = f'{hyd_year}-{self.hyd_month:02d}-{self.hyd_day:02d}'
        end = f'{hyd_year + 1}-{self.hyd_month:02d}-{self.hyd_day:02d}'
        return start, end

    def _mosaic_by_day(self, collection):
        """Mosaic all scenes acquired on the same calendar day.

        Sentinel-2 (and other sensors) can produce multiple scenes per day
        when adjacent tiles overlap the ROI. Without this step, same-day
        scenes get a midpoint distance of 0 and therefore weight = 0 days,
        silently discarding valid observations.

        For binary water masks (1 = water, 0 = dry) we use ``.max()``:
        if *any* tile shows water for a pixel on that day, the pixel is
        classified as water. The resulting image keeps ``system:time_start``
        set to midnight of that day.
        """
        # Get unique day-level timestamps (midnight UTC)
        def _floor_to_day(t):
            return ee.Date(t).update(hour=0, minute=0, second=0).millis()

        times = collection.aggregate_array('system:time_start')
        unique_days = times.map(_floor_to_day).distinct()

        def _mosaic_one_day(day_millis):
            d = ee.Date(day_millis)
            day_col = collection.filterDate(d, d.advance(1, 'day'))
            return (
                day_col.max()                        # water wins over dry
                .set('system:time_start', day_millis)
            )

        mosaics = unique_days.map(_mosaic_one_day)
        return ee.ImageCollection(mosaics)

    def _add_scene_weights(self, collection, d0):
        """Attach temporal weight properties to each image (server-side).

        Implements the midpoint weighting method:

        1. Sort images by acquisition date.
        2. Compute the midpoint (in days) between each pair of consecutive
           scenes.
        3. Each scene's territory spans ``[midpoint_prev, midpoint_next]``,
           with 0 and 365 as the outer boundaries for the first and last
           scenes.
        4. Set image properties: ``weight`` (days), ``start_doy``, ``end_doy``.

        All values are expressed as days elapsed since ``d0``.
        """
        collection = collection.sort('system:time_start')
        image_list = collection.toList(collection.size())
        n = collection.size()

        # Days from hydrological year start for each scene
        times = collection.aggregate_array('system:time_start')
        days = times.map(lambda t: ee.Date(t).difference(d0, 'day'))

        # Midpoints between consecutive scenes
        n_imgs = days.size()
        days_left = days.slice(0, n_imgs.subtract(1))
        days_right = days.slice(1)
        midpoints = days_left.zip(days_right).map(
            lambda pair: ee.List(pair).getNumber(0)
            .add(ee.List(pair).getNumber(1))
            .divide(2)
        )

        # Full boundary array: [0, mid_01, mid_12, …, mid_(n-2,n-1), 365]
        boundaries = ee.List([0]).cat(midpoints).cat([365])

        def _set_props(i):
            i = ee.Number(i)
            img = ee.Image(image_list.get(i))
            start_doy = boundaries.getNumber(i)
            end_doy = boundaries.getNumber(i.add(1))
            return img.set({
                'weight': end_doy.subtract(start_doy),
                'start_doy': start_doy,
                'end_doy': end_doy,
            })

        weighted = ee.List.sequence(0, n.subtract(1)).map(_set_props)
        return ee.ImageCollection(weighted)

    def _masks_params_for(self, image):
        """Return the (index, threshold, hyd_year) matching an exported image.

        :meth:`compute_hydroperiod` stamps ``index``, ``threshold`` and
        ``hyd_year_start`` on every result, so the masks that accompany an
        export can be rebuilt from the image itself rather than from whatever
        the last call happened to leave behind. Falls back to the cached values
        for images that carry no such metadata (a mean or anomaly composite).
        """
        info = image.toDictionary(
            ['index', 'threshold', 'hyd_year_start']
        ).getInfo() or {}

        index = info.get('index', self._last_index)
        threshold = info.get('threshold', self._last_threshold)
        hyd_year = info.get('hyd_year_start', self._hyd_year)

        if index is None:
            raise ValueError(
                "include_masks=True needs to know which index and threshold "
                "to use. The image carries no 'index'/'threshold' property "
                "and no compute_hydroperiod() has run on this analyzer. "
                "Call compute_hydroperiod() first, or export the masks "
                "explicitly with export_water_masks_to_drive(index=..., "
                "threshold=...)."
            )
        return index, threshold, hyd_year

    def _validate_index(self, index):
        """Raise informative errors for unsupported indices."""
        if index not in self.WATER_INDICES:
            raise ValueError(
                f"'{index}' is not a supported water index. "
                f"Choose from: {sorted(self.WATER_INDICES)}"
            )
        available = self.ns.sensor_indices.get(self.ns.sat, set())
        if index not in available:
            valid_water = sorted(self.WATER_INDICES & available)
            raise ValueError(
                f"Index '{index}' is not available for sensor '{self.ns.sat}'. "
                f"Supported water indices for {self.ns.sat}: {valid_water}"
            )

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def get_water_masks(self, index='mndwi', threshold=0.0, hyd_year=None,
                        add_footprint=False):
        """Generate weighted binary water masks from the satellite collection.

        For each scene in the hydrological year:

        1. Compute the selected water index.
        2. Classify pixels: ``index > threshold`` → water (1), else → dry (0).
        3. Preserve the cloud/nodata mask already applied by ``NdviSeasonality``
           (SCL-based for S2 when ``scl_mask=True``, QA60 otherwise).
        4. Mosaic same-day scenes (``max`` — water wins over dry).
        5. Attach temporal weight properties via the midpoint method.

        Parameters
        ----------
        index : str
            Water index to compute. One of:
            ``'ndwi'``, ``'mndwi'``, ``'awei'``, ``'aweinsh'``, ``'wi2015'``.
        threshold : float
            Classification threshold. Default ``0.0`` works well for most
            water indices (positive = water). Increase for stricter mapping.
        hyd_year : int, optional
            Starting calendar year of the hydrological cycle
            (e.g. ``2023`` for the Sep 2023 – Aug 2024 cycle).
            Defaults to ``ns.start_year``.
        add_footprint : bool
            Add a second band ``'footprint'`` (1 inside the acquisition
            footprint of that date, masked outside), which tells apart *"the
            satellite looked here but the pixel was masked as cloud"* from
            *"no scene covered this pixel that day"*. Used by
            :meth:`get_water_masks_stack`. Default ``False``.

        Returns
        -------
        ee.ImageCollection
            Binary water masks (band ``'water'``: 1 = water, 0 = dry,
            masked = nodata/cloud) with properties:
            ``weight`` (days), ``start_doy``, ``end_doy``.
            A ``'footprint'`` band is added when ``add_footprint=True``.

        Notes
        -----
        Pixel-level cloud/shadow masking is controlled by the ``scl_mask``
        parameter of :class:`NdviSeasonality` (default ``True`` for S2).
        """
        self._validate_index(index)

        if hyd_year is None:
            hyd_year = self.ns.start_year

        start_str, end_str = self._hyd_year_bounds(hyd_year)
        d0 = ee.Date(start_str)

        col = self.ns.ndvi_col.filterDate(start_str, end_str)
        index_fn = self.ns.d[index]

        def _to_mask(image):
            idx = index_fn(image).select(0)          # single-band, any name
            water = idx.gt(threshold).rename('water')
            water = water.updateMask(idx.mask())     # preserve nodata/cloud mask
            # The scene footprint survives updateMask, so it still describes
            # where the sensor actually acquired data — including the slanted
            # edges of Landsat scenes
            footprint = ee.Image(1).clip(image.geometry()).rename('footprint')
            return (
                water.addBands(footprint)
                .copyProperties(image, ['system:time_start', 'system:time_end'])
            )

        water_col = col.map(_to_mask)
        water_col = self._mosaic_by_day(water_col)   # deduplicate same-day tiles
        water_col = self._add_scene_weights(water_col, d0)

        if not add_footprint:
            water_col = water_col.select('water')

        self._water_masks = water_col
        self._hyd_year = hyd_year
        self._d0 = d0
        self._last_index = index
        self._last_threshold = threshold

        return water_col

    def get_water_masks_stack(
        self,
        index='mndwi',
        threshold=0.0,
        hyd_year=None,
        masked_value=2,
        nodata=255,
    ):
        """Stack the per-date water masks into a single 8-bit image.

        Turns the collection returned by :meth:`get_water_masks` into one
        multi-band ``ee.Image``, with **one band per acquisition date** named
        ``water_YYYY_MM_DD``. Values are stored as ``uint8``, which is four
        times smaller than the ``int16`` hydroperiod bands and is the whole
        point of exporting the masks separately: a GeoTIFF has a single data
        type for the entire file, so stacking the masks together with the
        hydroperiod bands would promote them to ``int16`` and waste the saving.

        Two different kinds of "no value" are distinguished, because they mean
        different things when interpreting a time series:

        - ``masked_value`` — the sensor **did** acquire the pixel that day, but
          it was discarded as cloud, shadow or another invalid class.
        - ``nodata`` — no scene covered the pixel at all: outside the
          acquisition footprint (the slanted edges of Landsat scenes, the gaps
          between Sentinel-2 orbits) or outside the ROI.

        Parameters
        ----------
        index : str
            Water index used to build the masks. Default ``'mndwi'``.
        threshold : float
            Water classification threshold. Default ``0.0``.
        hyd_year : int, optional
            Starting calendar year of the hydrological cycle.
            Defaults to ``ns.start_year``.
        masked_value : int
            Value for pixels observed but masked as cloud/shadow. Default ``2``.
        nodata : int
            Value for pixels no scene covered, and for everything outside the
            ROI. Default ``255``.

        Both must be integers in ``[0, 255]``, different from 0, from 1 and
        from each other.

        Returns
        -------
        ee.Image
            ``uint8`` image clipped to the ROI, one band per date:

            - ``0``            : dry
            - ``1``            : water
            - ``masked_value`` : observed but masked as cloud/shadow (2)
            - ``nodata``       : no acquisition, or outside the ROI (255)

            Image properties: ``'hyd_year_start'``, ``'hyd_year_end'``,
            ``'index'``, ``'threshold'``, ``'masked_value'``, ``'nodata'``.

        Examples
        --------
        ::

            stack = ha.get_water_masks_stack(index='mndwi', threshold=0.1)
            print(stack.bandNames().getInfo()[:3])
            # ['water_2022_09_04', 'water_2022_09_09', 'water_2022_09_14']

        See Also
        --------
        export_water_masks_to_drive : Send this stack to Google Drive.
        export_water_masks_to_asset : Send this stack to an Earth Engine asset.
        """
        for name, value in (('masked_value', masked_value), ('nodata', nodata)):
            if not isinstance(value, int) or isinstance(value, bool) \
                    or not 0 <= value <= 255:
                raise ValueError(
                    f"{name} must be an integer in [0, 255], got {value!r}."
                )
            if value in (0, 1):
                raise ValueError(
                    f"{name}={value} collides with the mask values "
                    f"(0 = dry, 1 = water). Use a different value."
                )
        if masked_value == nodata:
            raise ValueError(
                f"masked_value and nodata must differ, both are {nodata}. "
                f"They encode different things: masked_value means the pixel "
                f"was observed but discarded as cloud, nodata means no scene "
                f"covered it."
            )

        masks = self.get_water_masks(
            index=index, threshold=threshold, hyd_year=hyd_year,
            add_footprint=True,
        )

        def _encode(image):
            # 1 where any scene of that date covered the pixel, 0 elsewhere
            covered = image.select('footprint').unmask(0)
            # 0 dry / 1 water, then masked_value wherever the observation was
            # discarded, then nodata wherever nothing was acquired at all
            encoded = image.select('water').unmask(masked_value)
            return encoded.where(covered.Not(), nodata).toUint8()

        encoded_col = masks.map(_encode)

        # One band per date, named from system:time_start. aggregate_array
        # preserves collection order, and the collection is already sorted by
        # date and deduplicated by day, so names line up with toBands() output
        band_names = masks.aggregate_array('system:time_start').map(
            lambda t: ee.String('water_').cat(ee.Date(t).format('YYYY_MM_dd'))
        )

        stack = encoded_col.toBands().rename(band_names)

        # Clip first, then unmask: everything outside the ROI ends up carrying
        # the nodata value explicitly
        stack = stack.clip(self.ns.roi).unmask(nodata).toUint8()

        return stack.set({
            'hyd_year_start': self._hyd_year,
            'hyd_year_end': self._hyd_year + 1,
            'index': index,
            'threshold': threshold,
            'masked_value': masked_value,
            'nodata': nodata,
        })

    def compute_hydroperiod(
        self,
        index='mndwi',
        threshold=0.0,
        hyd_year=None,
        normalize=True,
        compute_first_last=True,
        min_flood_days=3,
        permanent_threshold=0.95,
    ):
        """Compute hydroperiod entirely in Google Earth Engine.

        Implements the midpoint weighting algorithm:

        1. Generate binary water masks (via :meth:`get_water_masks`).
        2. For each scene, compute weighted flood and valid-day contributions:

           - ``flood_contrib = water * weight``  (weight where flooded, 0 where dry)
           - ``valid_contrib = valid_mask * weight``  (weight where any valid obs)

        3. Sum contributions per pixel across all scenes →
           ``hydroperiod`` and ``valid_days``.
        4. Optionally normalise: ``normalized = (hydroperiod / valid_days) * 365``.
        5. Optionally compute first / last flood day (day of hydrological year).

        Parameters
        ----------
        index : str
            Water index for mask generation. Default ``'mndwi'``.
        threshold : float
            Water classification threshold. Default ``0.0``.
        hyd_year : int, optional
            Starting year of the hydrological cycle. Default ``ns.start_year``.
        normalize : bool
            Add a ``'normalized'`` band corrected for uneven temporal coverage.
            Default ``True``.
        compute_first_last : bool
            Add ``'first_flood_doy'`` and ``'last_flood_doy'`` bands (days since
            hydrological year start, 0-indexed). Default ``True``.
        min_flood_days : int
            Minimum accumulated flood days required for a pixel to receive a
            valid ``first_flood_doy`` / ``last_flood_doy`` value. Pixels below
            this threshold are masked (likely noise or isolated cloud artefacts).
            Default ``3``.
        permanent_threshold : float
            Fraction of valid days that must be flooded for a pixel to be
            considered permanently inundated (e.g. open sea, permanent lagoons).
            These pixels get ``first_flood_doy = 0`` and ``last_flood_doy = 365``
            because a meaningful onset/recession date cannot be determined.
            Range [0, 1]. Default ``0.95``.

        Returns
        -------
        ee.Image
            Multi-band image clipped to the ROI with bands:

            - ``hydroperiod``    : accumulated weighted flood days (float32)
            - ``valid_days``     : accumulated weighted valid days (float32)
            - ``normalized``     : flood days normalised to 365-day year (int16, optional)
            - ``first_flood_doy``: first flood day of cycle (int16, optional)
            - ``last_flood_doy`` : last flood day of cycle (int16, optional)

            Image properties: ``'hyd_year_start'``, ``'hyd_year_end'``,
            ``'index'``, ``'threshold'``.
        """
        water_masks = self.get_water_masks(
            index=index, threshold=threshold, hyd_year=hyd_year
        )

        # --- Weighted flood / valid contributions ---------------------------------
        def _contribution(image):
            weight = ee.Number(image.get('weight'))
            water = image.select('water')
            # valid: any pixel that has a real observation (not cloud/nodata)
            valid = water.mask().toFloat()

            flood_c = water.multiply(weight).toFloat().rename('flood')
            valid_c = valid.multiply(weight).rename('valid')
            return flood_c.addBands(valid_c)

        contributions = water_masks.map(_contribution)

        # Keep float sums for accurate division in normalization
        flood_sum = contributions.select('flood').sum()
        valid_sum = contributions.select('valid').sum()

        hydroperiod = flood_sum.round().toInt16().rename('hydroperiod')
        valid_days = valid_sum.round().toInt16().rename('valid_days')
        result = hydroperiod.addBands(valid_days)

        # --- Normalisation --------------------------------------------------------
        if normalize:
            normalized = (
                flood_sum
                .divide(valid_sum)
                .multiply(365)
                .round()
                .toInt16()
                .rename('normalized')
            )
            result = result.addBands(normalized)

        # --- First / last flood day -----------------------------------------------
        if compute_first_last:
            def _first_doy(image):
                doy = ee.Image.constant(
                    ee.Number(image.get('start_doy'))
                ).toFloat()
                return doy.updateMask(image.select('water')).rename('first_flood_doy')

            def _last_doy(image):
                doy = ee.Image.constant(
                    ee.Number(image.get('end_doy'))
                ).toFloat()
                return doy.updateMask(image.select('water')).rename('last_flood_doy')

            # Raw min/max DOY across all scenes where water was detected
            first_raw = water_masks.map(_first_doy).min()
            last_raw = water_masks.map(_last_doy).max()

            # Mask pixels with too few flood days (likely noise / cloud artefacts)
            enough_flood = flood_sum.gte(min_flood_days)

            # Permanently flooded pixels (sea, permanent lagoons): ratio ≈ 1
            # → assign fixed boundaries 0 / 365 (onset/recession not meaningful)
            permanent = flood_sum.divide(valid_sum).gte(permanent_threshold)

            first_flood = (
                first_raw
                .updateMask(enough_flood)
                .where(permanent, ee.Image.constant(0))
                .toInt16()
                .rename('first_flood_doy')
            )
            last_flood = (
                last_raw
                .updateMask(enough_flood)
                .where(permanent, ee.Image.constant(365))
                .toInt16()
                .rename('last_flood_doy')
            )
            result = result.addBands(first_flood).addBands(last_flood)

        # --- Metadata & clip ------------------------------------------------------
        hyd_year_val = self._hyd_year
        result = result.set({
            'hyd_year_start': hyd_year_val,
            'hyd_year_end': hyd_year_val + 1,
            'index': index,
            'threshold': threshold,
        }).clip(self.ns.roi)

        self._results = result
        return result

    def compute_irt_global(self, hyd_year=None):
        """Compute the global Temporal Representativity Index (IRT).

        Measures how uniformly the available scenes are distributed across
        the hydrological year. A value of 1 means perfect temporal coverage;
        0 means all scenes cluster in a single period.

        This uses the Gini-based formula identical to *phydroperiod*'s
        ``calculate_temporal_representativity()``, computed client-side
        from the scene date metadata.

        Parameters
        ----------
        hyd_year : int, optional
            Starting year of the hydrological cycle.

        Returns
        -------
        float
            IRT value between 0 (worst) and 1 (best). Requires one
            ``.getInfo()`` call to retrieve scene dates.
        """
        if hyd_year is None:
            hyd_year = self.ns.start_year

        start_str, end_str = self._hyd_year_bounds(hyd_year)
        times = (
            self.ns.ndvi_col
            .filterDate(start_str, end_str)
            .aggregate_array('system:time_start')
            .getInfo()
        )

        if not times:
            return 0.0

        n_periods = 12
        period_length = 365.0 / n_periods
        d0 = date(hyd_year, self.hyd_month, self.hyd_day)

        period_counts = np.zeros(n_periods)
        for t_ms in times:
            d1 = date.fromtimestamp(t_ms / 1000)
            doy = (d1 - d0).days % 365
            idx = min(int(doy / period_length), n_periods - 1)
            period_counts[idx] += 1

        # Gini coefficient (same formula as phydroperiod)
        sorted_counts = np.sort(period_counts)
        n = len(sorted_counts)
        cumsum = np.cumsum(sorted_counts)
        total = float(cumsum[-1])
        if total == 0:
            return 0.0
        gini = 1 - 2 * np.sum(cumsum) / (n * total) + 1 / n
        return round(float(1 - gini), 4)

    def compute_irt_image(self, hyd_year=None, n_periods=12):
        """Compute per-pixel Temporal Representativity Index as an ``ee.Image``.

        For each pixel, estimates how uniformly its valid (non-cloud)
        observations are distributed across the hydrological year. This is
        the GEE-native counterpart of *phydroperiod*'s
        ``calculate_pixel_irt()``.

        Uses Simpson's diversity index as a per-pixel-computable proxy:

        .. math::

            \\text{IRT}_{pixel} = \\frac{N^2}{n_{periods} \\cdot \\sum n_i^2}

        where :math:`N` is the total number of valid observations for that
        pixel and :math:`n_i` is the count in period *i*. This ranges from
        :math:`1/n_{periods}` (all observations in one period) to 1 (perfectly
        uniform).

        Parameters
        ----------
        hyd_year : int, optional
            Starting year of the hydrological cycle.
        n_periods : int
            Number of periods the year is divided into. Default 12 (monthly).

        Returns
        -------
        ee.Image
            Single-band IRT image (band ``'irt'``), values 0–1, clipped to
            ROI. Masked where no valid observations exist.

        Notes
        -----
        Call :meth:`get_water_masks` first (or let this method call it
        automatically via the ``hyd_year`` argument).
        """
        if self._water_masks is None:
            self.get_water_masks(hyd_year=hyd_year)

        # Select explicitly: the cached collection carries an extra 'footprint'
        # band when it was last built by get_water_masks_stack()
        water_masks = self._water_masks.select('water')
        period_length = 365.0 / n_periods

        # For each period: count valid (unmasked) observations per pixel
        period_images = []
        for p in range(n_periods):
            p_start = p * period_length
            p_end = (p + 1) * period_length

            period_col = water_masks.filter(
                ee.Filter.And(
                    ee.Filter.gte('start_doy', p_start),
                    ee.Filter.lt('start_doy', p_end),
                )
            )

            valid_count = (
                period_col
                .map(lambda img: img.mask().toFloat())
                .sum()
                .rename(f'p{p}')
            )
            period_images.append(valid_count)

        period_stack = ee.Image.cat(period_images)

        # Simpson's diversity: IRT = N² / (n_periods × Σ nᵢ²)
        total_obs = period_stack.reduce(ee.Reducer.sum())
        sum_sq = period_stack.pow(2).reduce(ee.Reducer.sum())

        irt = (
            total_obs.pow(2)
            .divide(sum_sq.multiply(n_periods))
            .rename('irt')
            .updateMask(total_obs.gt(0))
            .clip(self.ns.roi)
        )

        return irt

    def compute_all_cycles(
        self,
        index='mndwi',
        threshold=0.0,
        normalize=True,
        compute_first_last=True,
    ):
        """Compute hydroperiod for every hydrological cycle in the date range.

        Iterates from ``ns.start_year`` to ``ns.end_year`` (inclusive),
        computing one hydroperiod image per cycle. Each cycle spans from
        ``hydrological_year_start`` in year *Y* to the same date in year *Y+1*.

        Parameters
        ----------
        index : str
            Water index for mask generation. Default ``'mndwi'``.
        threshold : float
            Water classification threshold. Default ``0.0``.
        normalize : bool
            Include ``'normalized'`` band in each cycle image. Default ``True``.
        compute_first_last : bool
            Include ``'first_flood_doy'`` / ``'last_flood_doy'`` bands.
            Default ``True``.

        Returns
        -------
        dict
            Mapping ``{hyd_year: ee.Image}`` where *hyd_year* is the starting
            calendar year of each cycle (e.g. ``2020`` for 2020-2021).
            Each image has the same bands as :meth:`compute_hydroperiod` plus
            a ``'hyd_year'`` property for easy identification.

        Examples
        --------
        ::

            ns = NdviSeasonality(roi, sat='S2', start_year=2019, end_year=2023)
            ha = HydroperiodAnalyzer(ns)

            cycles = ha.compute_all_cycles(index='mndwi', threshold=0.1)
            # cycles = {2019: ee.Image, 2020: ee.Image, ..., 2023: ee.Image}

            # Export all cycles to Drive
            for yr, img in cycles.items():
                ha.export_to_drive(image=img, description=f'hydroperiod_{yr}_{yr+1}')
        """
        results = {}
        for yr in range(self.ns.start_year, self.ns.end_year + 1):
            print(f"Computing hydroperiod cycle {yr}/{yr + 1}…")
            img = self.compute_hydroperiod(
                index=index,
                threshold=threshold,
                hyd_year=yr,
                normalize=normalize,
                compute_first_last=compute_first_last,
            )
            results[yr] = img
        return results

    #: First available hydrological year for each sensor's full archive.
    HISTORICAL_START = {
        'S2': 2017,       # S2_SR_HARMONIZED quality data from 2017
        'Landsat': 1984,  # Landsat 5 TM
        'MODIS': 2000,    # MOD09A1 Terra
        'S1': 2014,       # Sentinel-1A launch
    }

    def _last_complete_hyd_year(self):
        """Return the last hydrological year whose cycle is fully complete."""
        from datetime import date as _date
        today = _date.today()
        if (today.month, today.day) > (self.hyd_month, self.hyd_day):
            return today.year - 1
        return today.year - 2

    def compute_anomalies(self, cycles=None, index='mndwi', threshold=0.0,
                          reference='period'):
        """Compute mean hydroperiod and per-cycle anomalies.

        Calculates the mean normalised hydroperiod across a reference period
        and the anomaly of each cycle relative to that mean:

        ``anomaly = cycle_normalized - mean_normalized``

        Positive values indicate a wetter-than-average year; negative values
        indicate a drier-than-average year.

        Parameters
        ----------
        cycles : dict, optional
            Output of :meth:`compute_all_cycles`: ``{hyd_year: ee.Image}``.
            If ``None``, runs :meth:`compute_all_cycles` automatically.
        index : str
            Water index, used only when ``cycles`` is ``None``.
        threshold : float
            Classification threshold, used only when ``cycles`` is ``None``.
        reference : str
            Reference period for computing the mean:

            - ``'period'`` *(default)*: mean of the cycles in the current
              analysis range (``ns.start_year`` → ``ns.end_year``).
            - ``'historical'``: mean of the full satellite archive, from
              the sensor's first available year to the last complete cycle.
              Gives a climatological baseline independent of the chosen
              analysis window. Note: GEE builds the full graph lazily so
              this call is fast, but rendering/export will be slower.

        Returns
        -------
        dict
            ``{'mean': ee.Image, 'anomalies': {hyd_year: ee.Image}}``

            - ``'mean'``: mean normalised hydroperiod (band
              ``'mean_hydroperiod'``, int16). Properties: ``'reference'``,
              ``'ref_start_year'``, ``'ref_end_year'``, ``'ref_n_cycles'``.
            - ``'anomalies'``: per-cycle anomaly images (band ``'anomaly'``,
              int16). Same metadata properties set on each image.

        Examples
        --------
        ::

            # Mean of the analysis period
            result = ha.compute_anomalies(cycles, reference='period')

            # Climatological baseline (full satellite archive)
            result = ha.compute_anomalies(cycles, reference='historical')

            Map.addLayer(result['mean'], vis, 'Mean hydroperiod')
            for yr, anom in result['anomalies'].items():
                Map.addLayer(anom, anom_vis, f'Anomaly {yr}/{yr+1}')
        """
        if cycles is None:
            cycles = self.compute_all_cycles(index=index, threshold=threshold)

        if reference == 'period':
            ref_images = [img.select('normalized') for img in cycles.values()]
            ref_years = sorted(cycles.keys())

        elif reference == 'historical':
            hist_start = self.HISTORICAL_START.get(self.ns.sat)
            if hist_start is None:
                raise ValueError(
                    f"No historical start defined for sensor '{self.ns.sat}'. "
                    f"Available sensors: {list(self.HISTORICAL_START)}"
                )
            last_complete = self._last_complete_hyd_year()
            ref_years = list(range(hist_start, last_complete + 1))

            print(
                f"Building historical reference: {len(ref_years)} cycles "
                f"({hist_start}/{hist_start+1} → {last_complete}/{last_complete+1}). "
                f"GEE evaluates lazily — computation happens on render/export."
            )

            # Reuse already-computed cycles to avoid redundant graph nodes
            hist_cycles = dict(cycles)
            for yr in ref_years:
                if yr not in hist_cycles:
                    hist_cycles[yr] = self.compute_hydroperiod(
                        index=index, threshold=threshold, hyd_year=yr
                    )
            ref_images = [hist_cycles[yr].select('normalized') for yr in ref_years]

        else:
            raise ValueError(
                f"reference must be 'period' or 'historical', got '{reference!r}'"
            )

        mean_hyd = (
            ee.ImageCollection(ref_images)
            .mean()
            .round()
            .toInt16()
            .rename('mean_hydroperiod')
            .set({
                'reference': reference,
                'ref_start_year': ref_years[0],
                'ref_end_year': ref_years[-1],
                'ref_n_cycles': len(ref_years),
            })
        )

        anomalies = {
            yr: (img.select('normalized')
                 .subtract(mean_hyd)
                 .rename('anomaly')
                 .toInt16()
                 .set({
                     'hyd_year_start': yr,
                     'hyd_year_end': yr + 1,
                     'reference': reference,
                     'ref_start_year': ref_years[0],
                     'ref_end_year': ref_years[-1],
                     'ref_n_cycles': len(ref_years),
                 }))
            for yr, img in cycles.items()
        }

        return {'mean': mean_hyd, 'anomalies': anomalies}

    # ------------------------------------------------------------------
    # Export helpers
    # ------------------------------------------------------------------

    def export_water_masks_to_drive(
        self,
        index='mndwi',
        threshold=0.0,
        hyd_year=None,
        masked_value=2,
        nodata=255,
        stack=None,
        folder='hydroperiod',
        description=None,
        scale=10,
        crs='EPSG:4326',
        **kwargs,
    ):
        """Export the per-date water masks to Google Drive as an 8-bit stack.

        Sends the image built by :meth:`get_water_masks_stack` to Drive as a
        single ``uint8`` GeoTIFF with one band per acquisition date. This is a
        separate file from the hydroperiod export on purpose: a GeoTIFF has one
        data type for the whole file, so keeping the masks apart is what allows
        them to stay 8-bit instead of being promoted to the ``int16`` of the
        hydroperiod bands.

        Parameters
        ----------
        index : str
            Water index used to build the masks. Default ``'mndwi'``.
        threshold : float
            Water classification threshold. Default ``0.0``.
        hyd_year : int, optional
            Starting calendar year of the hydrological cycle.
            Defaults to ``ns.start_year``.
        masked_value : int
            Value for pixels observed but masked as cloud/shadow. Default ``2``.
        nodata : int
            Value for pixels no scene covered, and for everything outside the
            ROI. Default ``255``.
        stack : ee.Image, optional
            Pre-built mask stack to export. If given, ``index``, ``threshold``,
            ``hyd_year``, ``masked_value`` and ``nodata`` are ignored.
        folder : str
            Google Drive destination folder. Default ``'hydroperiod'``.
        description : str, optional
            Task name in the Earth Engine Tasks panel.
            Auto-generated from the hydrological year if not provided.
        scale : int
            Output pixel size in metres. Default 10 (Sentinel-2).
        crs : str
            Output CRS. Default ``'EPSG:4326'``.
        **kwargs
            Additional keyword arguments forwarded to
            ``ee.batch.Export.image.toDrive``.

        Returns
        -------
        ee.batch.Task
            Started export task.

        Notes
        -----
        The nodata value is *not* written into the GeoTIFF header — Earth
        Engine does not tag it. Set it yourself when reading the file
        (``rasterio.open(...).read(masked=True)`` after declaring it, or
        Properties → No Data in QGIS).

        Examples
        --------
        ::

            ha.export_water_masks_to_drive(
                index='mndwi', threshold=0.1,
                folder='my_wetland', scale=10, crs='EPSG:25829',
            )
        """
        if stack is None:
            stack = self.get_water_masks_stack(
                index=index, threshold=threshold, hyd_year=hyd_year,
                masked_value=masked_value, nodata=nodata,
            )
        else:
            # A pre-built stack carries its own encoding; read it back so the
            # message below describes the file that is actually being written
            encoding = stack.toDictionary(['masked_value', 'nodata']).getInfo()
            masked_value = encoding.get('masked_value', masked_value)
            nodata = encoding.get('nodata', nodata)

        if description is None:
            yr = self._hyd_year or self.ns.start_year
            description = f'water_masks_{yr}_{yr + 1}'

        task = ee.batch.Export.image.toDrive(
            image=stack,
            description=description,
            folder=folder,
            region=self.ns.roi,
            scale=scale,
            crs=crs,
            maxPixels=1e13,
            **kwargs,
        )
        task.start()
        n_dates = stack.bandNames().size().getInfo()
        print(
            f"Export task '{description}' started → Drive folder '{folder}' "
            f"({n_dates} dates, uint8: 0=dry, 1=water, "
            f"{masked_value}=cloud-masked, {nodata}=no acquisition)."
        )
        return task

    def export_to_drive(
        self,
        image=None,
        folder='hydroperiod',
        description=None,
        scale=10,
        crs='EPSG:4326',
        include_masks=False,
        **kwargs,
    ):
        """Export hydroperiod results to Google Drive.

        Parameters
        ----------
        image : ee.Image, optional
            Image to export. If ``None``, uses the last
            :meth:`compute_hydroperiod` result.
        folder : str
            Google Drive destination folder. Default ``'hydroperiod'``.
        description : str, optional
            Task name in the Earth Engine Tasks panel.
            Auto-generated from hydrological year if not provided.
        scale : int
            Output pixel size in metres. Default 10 (Sentinel-2).
        crs : str
            Output CRS. Default ``'EPSG:4326'``.
        include_masks : bool
            Also export the per-date water masks as a separate ``uint8``
            file (0 = dry, 1 = water, 2 = cloud-masked, 255 = no acquisition).
            The index, threshold and cycle are read from the exported image
            itself, so the masks always match what is being exported. They
            cannot share the hydroperiod file: a GeoTIFF has a single data
            type, so bundling them would promote the masks to ``int16``.
            Default ``False``.
        **kwargs
            Additional keyword arguments forwarded to
            ``ee.batch.Export.image.toDrive``.

        Returns
        -------
        ee.batch.Task or tuple of ee.batch.Task
            Started export task. When ``include_masks=True``, a
            ``(hydroperiod_task, masks_task)`` tuple is returned instead.
        """
        if image is None:
            if self._results is None:
                raise ValueError(
                    "No results available. Run compute_hydroperiod() first "
                    "or pass an ee.Image explicitly."
                )
            import warnings
            warnings.warn(
                "Exporting self._results (last compute_hydroperiod() call). "
                "If you called compute_hydroperiod() multiple times (e.g. to "
                "compare indices or thresholds), pass the image explicitly: "
                "export_to_drive(image=result, ...) to avoid exporting the "
                "wrong result.",
                UserWarning,
                stacklevel=2,
            )
            image = self._results

        if description is None:
            yr = self._hyd_year or self.ns.start_year
            description = f'hydroperiod_{yr}_{yr + 1}'

        task = ee.batch.Export.image.toDrive(
            image=image,
            description=description,
            folder=folder,
            region=self.ns.roi,
            scale=scale,
            crs=crs,
            maxPixels=1e13,
            **kwargs,
        )
        task.start()
        print(f"Export task '{description}' started → Drive folder '{folder}'.")

        if include_masks:
            masks_index, masks_threshold, masks_year = \
                self._masks_params_for(image)
            # Both exports share kwargs, so a user-supplied file name has to be
            # suffixed too or the two files would overwrite each other in Drive
            masks_kwargs = dict(kwargs)
            if 'fileNamePrefix' in masks_kwargs:
                masks_kwargs['fileNamePrefix'] += '_water_masks'
            masks_task = self.export_water_masks_to_drive(
                index=masks_index,
                threshold=masks_threshold,
                hyd_year=masks_year,
                folder=folder,
                description=f'{description}_water_masks',
                scale=scale,
                crs=crs,
                **masks_kwargs,
            )
            return task, masks_task

        return task

    def export_water_masks_to_asset(
        self,
        asset_id,
        index='mndwi',
        threshold=0.0,
        hyd_year=None,
        masked_value=2,
        nodata=255,
        stack=None,
        description=None,
        scale=10,
        crs='EPSG:4326',
        **kwargs,
    ):
        """Export the per-date water masks to an Earth Engine Asset.

        Asset counterpart of :meth:`export_water_masks_to_drive`. Writes the
        image built by :meth:`get_water_masks_stack` as a ``uint8`` asset with
        one band per acquisition date.

        Parameters
        ----------
        asset_id : str
            Full EE asset path
            (e.g. ``'projects/my-project/assets/water_masks_2023_2024'``).
        index : str
            Water index used to build the masks. Default ``'mndwi'``.
        threshold : float
            Water classification threshold. Default ``0.0``.
        hyd_year : int, optional
            Starting calendar year of the hydrological cycle.
            Defaults to ``ns.start_year``.
        masked_value : int
            Value for pixels observed but masked as cloud/shadow. Default ``2``.
        nodata : int
            Value for pixels no scene covered, and for everything outside the
            ROI. Default ``255``.
        stack : ee.Image, optional
            Pre-built mask stack to export. If given, ``index``, ``threshold``,
            ``hyd_year``, ``masked_value`` and ``nodata`` are ignored.
        description : str, optional
            Task name. Auto-generated from the hydrological year if not
            provided.
        scale : int
            Output pixel size in metres. Default 10 (Sentinel-2).
        crs : str
            Output CRS. Default ``'EPSG:4326'``.
        **kwargs
            Additional keyword arguments forwarded to
            ``ee.batch.Export.image.toAsset``.

        Returns
        -------
        ee.batch.Task
            Started export task.

        Notes
        -----
        Unlike a GeoTIFF, an Earth Engine asset keeps one data type per band,
        so this could technically live in the same asset as the hydroperiod
        bands. It is kept separate anyway, to mirror the Drive export and to
        keep the date stack independent of the hydroperiod result.
        """
        if stack is None:
            stack = self.get_water_masks_stack(
                index=index, threshold=threshold, hyd_year=hyd_year,
                masked_value=masked_value, nodata=nodata,
            )

        if description is None:
            yr = self._hyd_year or self.ns.start_year
            description = f'water_masks_{yr}_{yr + 1}'

        task = ee.batch.Export.image.toAsset(
            image=stack,
            description=description,
            assetId=asset_id,
            region=self.ns.roi,
            scale=scale,
            crs=crs,
            maxPixels=1e13,
            **kwargs,
        )
        task.start()
        n_dates = stack.bandNames().size().getInfo()
        print(
            f"Export task '{description}' started → Asset '{asset_id}' "
            f"({n_dates} dates, uint8: 0=dry, 1=water, "
            f"{masked_value}=cloud-masked, {nodata}=no acquisition)."
        )
        return task

    def export_to_asset(
        self,
        asset_id,
        image=None,
        description=None,
        scale=10,
        crs='EPSG:4326',
        include_masks=False,
        masks_asset_id=None,
        **kwargs,
    ):
        """Export hydroperiod results to an Earth Engine Asset.

        Parameters
        ----------
        asset_id : str
            Full EE asset path
            (e.g. ``'projects/my-project/assets/hydroperiod_2023_2024'``).
        image : ee.Image, optional
            Image to export. If ``None``, uses the last
            :meth:`compute_hydroperiod` result.
        description : str, optional
            Task name. Auto-generated if not provided.
        scale : int
            Output pixel size in metres. Default 10.
        crs : str
            Output CRS. Default ``'EPSG:4326'``.
        include_masks : bool
            Also export the per-date water masks as a separate ``uint8``
            asset (0 = dry, 1 = water, 2 = cloud-masked, 255 = no
            acquisition). The index, threshold and cycle are read from the
            exported image itself, so the masks always match what is being
            exported. Default ``False``.
        masks_asset_id : str, optional
            Asset path for the masks when ``include_masks=True``.
            Defaults to ``asset_id`` with a ``'_water_masks'`` suffix.
        **kwargs
            Additional keyword arguments forwarded to
            ``ee.batch.Export.image.toAsset``.

        Returns
        -------
        ee.batch.Task or tuple of ee.batch.Task
            Started export task. When ``include_masks=True``, a
            ``(hydroperiod_task, masks_task)`` tuple is returned instead.
        """
        if image is None:
            if self._results is None:
                raise ValueError(
                    "No results available. Run compute_hydroperiod() first "
                    "or pass an ee.Image explicitly."
                )
            import warnings
            warnings.warn(
                "Exporting self._results (last compute_hydroperiod() call). "
                "If you called compute_hydroperiod() multiple times (e.g. to "
                "compare indices or thresholds), pass the image explicitly: "
                "export_to_asset(image=result, ...) to avoid exporting the "
                "wrong result.",
                UserWarning,
                stacklevel=2,
            )
            image = self._results

        if description is None:
            yr = self._hyd_year or self.ns.start_year
            description = f'hydroperiod_{yr}_{yr + 1}'

        task = ee.batch.Export.image.toAsset(
            image=image,
            description=description,
            assetId=asset_id,
            region=self.ns.roi,
            scale=scale,
            crs=crs,
            maxPixels=1e13,
            **kwargs,
        )
        task.start()
        print(
            f"Export task '{description}' started → Asset '{asset_id}'."
        )

        if include_masks:
            masks_index, masks_threshold, masks_year = \
                self._masks_params_for(image)
            masks_task = self.export_water_masks_to_asset(
                asset_id=masks_asset_id or f'{asset_id}_water_masks',
                index=masks_index,
                threshold=masks_threshold,
                hyd_year=masks_year,
                description=f'{description}_water_masks',
                scale=scale,
                crs=crs,
                **kwargs,
            )
            return task, masks_task

        return task

    def __repr__(self):
        yr = self._hyd_year or self.ns.start_year
        return (
            f"HydroperiodAnalyzer("
            f"sat={self.ns.sat!r}, "
            f"hyd_year={yr}/{yr + 1}, "
            f"hyd_start={self.hyd_month:02d}-{self.hyd_day:02d})"
        )
