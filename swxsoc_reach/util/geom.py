"""Geometry utilities for REACH region operations."""

from __future__ import annotations

import heapq
from pathlib import Path

import matplotlib.path as mpath
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
from scipy.ndimage import gaussian_filter
from spherical_geometry.polygon import SingleSphericalPolygon, SphericalPolygon

from swxsoc_reach import _data_directory
from swxsoc_reach.util.enums import Region

REGION_CODES = {
    code: region.label
    for region in Region.ordered()
    for code in (-region.code, region.code)
}


def contour_image_to_path(
    image: npt.ArrayLike,
    contour_levels: list[float],
    x_values: npt.ArrayLike | None = None,
    y_values: npt.ArrayLike | None = None,
    blur_sigma: float = 0.0,
) -> dict[float, mpath.Path | None]:
    """Convert one contour level from an image into a matplotlib path.

    Parameters
    ----------
    image : array-like
        2-D array to contour.
    contour_levels : list of float
        The values to contour at.
    x_values : array-like, optional
        1-D array of x coordinates (columns). Defaults to pixel indices.
    y_values : array-like, optional
        1-D array of y coordinates (rows). Defaults to pixel indices.
    blur_sigma : float, optional
        Standard deviation for Gaussian blur applied to the boolean mask
        before contouring. ``0.0`` (default) disables blurring.

    Returns
    -------
    dict[float, matplotlib.path.Path | None]
        Mapping from each requested contour level to a compound path of all
        contour segments for that level, or ``None`` if no contour is found.
    """
    image2d = np.asarray(image, dtype=float)
    if image2d.ndim != 2:
        raise ValueError("image must be a 2D array.")

    if x_values is None:
        x_arr = np.arange(image2d.shape[1], dtype=float)
    else:
        x_arr = np.asarray(x_values, dtype=float)

    if y_values is None:
        y_arr = np.arange(image2d.shape[0], dtype=float)
    else:
        y_arr = np.asarray(y_values, dtype=float)

    finite_values = np.unique(image2d[np.isfinite(image2d)])
    if finite_values.size == 0:
        return {float(level): None for level in contour_levels}

    finite_values = np.sort(finite_values.astype(float))
    if finite_values.size == 1:
        half_step = 0.5
        band_edges = np.array(
            [finite_values[0] - half_step, finite_values[0] + half_step],
            dtype=float,
        )
    else:
        midpoints = 0.5 * (finite_values[:-1] + finite_values[1:])
        lower_edge = finite_values[0] - 0.5 * (finite_values[1] - finite_values[0])
        upper_edge = finite_values[-1] + 0.5 * (finite_values[-1] - finite_values[-2])
        band_edges = np.concatenate(([lower_edge], midpoints, [upper_edge]))

    if blur_sigma > 0.0:
        contour_data = gaussian_filter(image2d, sigma=blur_sigma)
        contour_x = x_arr
        contour_y = y_arr
    else:
        # Pad the image so contours touching domain edges become closed polygons.
        # Use a fill value strictly below all band edges so padded cells are
        # unambiguously outside every band (avoids the 0.0 == midpoint(-1,+1) trap).
        fill_value = float(band_edges[0]) - 1.0
        dx = float(np.median(np.diff(x_arr))) if x_arr.size > 1 else 1.0
        dy = float(np.median(np.diff(y_arr))) if y_arr.size > 1 else 1.0
        contour_data = np.pad(image2d, 1, mode="constant", constant_values=fill_value)
        contour_x = np.concatenate(([x_arr[0] - dx], x_arr, [x_arr[-1] + dx]))
        contour_y = np.concatenate(([y_arr[0] - dy], y_arr, [y_arr[-1] + dy]))

    fig, ax = plt.subplots()
    contour = ax.contourf(contour_x, contour_y, contour_data, levels=band_edges)

    contour_paths: dict[float, mpath.Path | None] = {
        float(level): None for level in contour_levels
    }
    for contour_level, segments in zip(finite_values, contour.allsegs, strict=False):
        if float(contour_level) not in contour_paths:
            continue

        paths: list[mpath.Path] = []
        for segment in segments:
            vertices = np.asarray(segment, dtype=float)
            if vertices.shape[0] < 3:
                continue
            paths.append(mpath.Path(vertices, closed=True))

        if not paths:
            contour_paths[float(contour_level)] = None
        elif len(paths) == 1:
            contour_paths[float(contour_level)] = paths[0]
        else:
            contour_paths[float(contour_level)] = mpath.Path.make_compound_path(*paths)

    plt.close(fig)
    return contour_paths


def read_contour_path(
    contour_file: str | Path | None = None,
) -> mpath.Path:
    """
    Read saved contour NPZ and return matplotlib path.

    Parameters
    ----------
    contour_file : str | Path | None, optional
        Path to contour NPZ file. If None, uses package data file.

    Notes
    -----
    This function reads ``region_contour_paths.npz`` from the package data
    directory.

    Returns
    -------
    matplotlib.path.Path
        Matplotlib path object
    """

    contour_path = (
        Path(contour_file)
        if contour_file is not None
        else _data_directory / "region_contour_paths.npz"
    )
    data = np.load(contour_path, allow_pickle=False)
    vertices = np.asarray(data["vertices"], dtype=float)
    codes = data["codes"] if "codes" in data.files else None
    if codes is not None:
        codes = np.asarray(codes, dtype=np.uint8)
        if codes.size == 0:
            codes = None

    if codes is None:
        return mpath.Path(vertices, closed=True)
    return mpath.Path(vertices, codes=codes)


def load_region_contours(
    contour_file: str | Path | None = None,
) -> dict[int, mpath.Path]:
    """Load contour paths for all regions from the combined package data file.

    Parameters
    ----------
    contour_file : str | Path | None, optional
        Path to the combined NPZ file written by the region-to-contours notebook.
        Defaults to ``region_contour_paths.npz`` in the package data directory.

    Returns
    -------
    dict[int, matplotlib.path.Path]
        Mapping from region code to a (possibly compound) matplotlib Path.
    """
    npz_path = (
        Path(contour_file)
        if contour_file is not None
        else _data_directory / "region_contour_paths.npz"
    )
    data = np.load(npz_path, allow_pickle=False)

    contour_levels = np.asarray(data["contour_levels"], dtype=int)
    all_vertices = np.asarray(data["vertices"], dtype=float)
    all_codes = (
        np.asarray(data["codes"], dtype=np.uint8)
        if "codes" in data.files
        else np.empty((0,), dtype=np.uint8)
    )
    path_vertex_counts = np.asarray(data["path_vertex_counts"], dtype=int)
    path_code_counts = np.asarray(data["path_code_counts"], dtype=int)

    paths_by_region: dict[int, list[mpath.Path]] = {}
    vert_offset = 0
    code_offset = 0

    for level, n_verts, n_codes in zip(
        contour_levels, path_vertex_counts, path_code_counts, strict=False
    ):
        verts = all_vertices[vert_offset : vert_offset + n_verts]
        vert_offset += n_verts

        if n_codes > 0:
            codes = all_codes[code_offset : code_offset + n_codes]
            code_offset += n_codes
            path = mpath.Path(verts, codes=codes)
        else:
            path = mpath.Path(verts, closed=True)

        region = int(level)
        paths_by_region.setdefault(region, []).append(path)

    result: dict[int, mpath.Path] = {}
    for region, region_paths in paths_by_region.items():
        if len(region_paths) == 1:
            result[region] = region_paths[0]
        else:
            result[region] = mpath.Path.make_compound_path(*region_paths)

    return result


def save_path_to_npz(
    contour_map: dict[int, mpath.Path | list[mpath.Path] | None],
    file_path: str | Path,
) -> Path:
    """Save region contour paths into a single compact NPZ file.

    Parameters
    ----------
    contour_map : dict[int, matplotlib.path.Path | list[matplotlib.path.Path] | None]
        Mapping from region code to either one path, many paths, or ``None``.
    file_path : str | Path
        Output NPZ file path.

    Returns
    -------
    pathlib.Path
        The output path that was written.
    """
    output_path = Path(file_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    contour_levels: list[int] = []
    vertex_arrays: list[npt.NDArray[np.float64]] = []
    code_arrays: list[npt.NDArray[np.uint8]] = []
    path_vertex_counts: list[int] = []
    path_code_counts: list[int] = []

    for contour_level, contour_value in sorted(contour_map.items()):
        if contour_value is None:
            continue

        contour_paths = (
            [contour_value]
            if isinstance(contour_value, mpath.Path)
            else list(contour_value)
        )

        for contour_path in contour_paths:
            vertices = np.asarray(contour_path.vertices, dtype=float)
            vertex_arrays.append(vertices)
            path_vertex_counts.append(int(vertices.shape[0]))
            contour_levels.append(int(contour_level))

            if contour_path.codes is None:
                path_code_counts.append(0)
                continue

            codes = np.asarray(contour_path.codes, dtype=np.uint8)
            code_arrays.append(codes)
            path_code_counts.append(int(codes.shape[0]))

    if not contour_levels:
        raise ValueError("At least one contour path is required.")

    np.savez_compressed(
        output_path,
        contour_levels=np.asarray(contour_levels, dtype=np.int32),
        vertices=np.concatenate(vertex_arrays, axis=0),
        codes=(
            np.concatenate(code_arrays, axis=0)
            if code_arrays
            else np.empty((0,), dtype=np.uint8)
        ),
        path_vertex_counts=np.asarray(path_vertex_counts, dtype=np.int32),
        path_code_counts=np.asarray(path_code_counts, dtype=np.int32),
    )

    return output_path


def points_to_region_code(
    lon: float | npt.NDArray,
    lat: float | npt.NDArray,
    paths_dict: dict[int, mpath.Path],
) -> int | npt.NDArray | None:
    """Determine region code(s) for longitude/latitude points using indexed paths.

    Path subpaths are indexed into 10-degree longitude/latitude bins. Queries
    only test paths whose bounding boxes intersect the corresponding bins.
    """
    lon_arr = np.atleast_1d(np.asarray(lon, dtype=float))
    lat_arr = np.atleast_1d(np.asarray(lat, dtype=float))
    if lon_arr.shape != lat_arr.shape:
        raise ValueError("lon and lat must have the same shape.")

    points = np.column_stack((lon_arr, lat_arr))
    region_codes = np.full(points.shape[0], fill_value=np.nan)
    bin_size_degrees = 10.0
    longitude_bin_count = int(360.0 / bin_size_degrees)
    latitude_bin_count = int(180.0 / bin_size_degrees)
    paths_by_bin: list[list[int]] = [
        [] for _ in range(longitude_bin_count * latitude_bin_count)
    ]
    indexed_subpaths: list[tuple[int, mpath.Path, float, float, float, float]] = []

    def _iter_subpaths(path_obj: mpath.Path) -> list[mpath.Path]:
        """Treat compound paths as unions of independent closed subpaths."""
        codes = path_obj.codes
        if codes is None:
            return [path_obj]

        codes_arr = np.asarray(codes)
        move_indices = np.flatnonzero(codes_arr == mpath.Path.MOVETO)
        if move_indices.size <= 1:
            return [path_obj]

        vertices_arr = np.asarray(path_obj.vertices, dtype=float)
        subpaths = []
        for idx, start in enumerate(move_indices):
            end = (
                move_indices[idx + 1]
                if idx + 1 < move_indices.size
                else len(vertices_arr)
            )
            sub_vertices = vertices_arr[start:end]
            sub_codes = codes_arr[start:end]
            if sub_vertices.shape[0] < 3:
                continue
            subpaths.append(mpath.Path(sub_vertices, codes=sub_codes))
        return subpaths

    for code, path_obj in paths_dict.items():
        if path_obj is None:
            continue
        for subpath in _iter_subpaths(path_obj):
            vertices = np.asarray(subpath.vertices, dtype=float)
            finite_vertices = vertices[np.all(np.isfinite(vertices), axis=1)]
            if finite_vertices.shape[0] < 3:
                continue
            min_lon, min_lat = np.min(finite_vertices, axis=0)
            max_lon, max_lat = np.max(finite_vertices, axis=0)
            if max_lon < -180.0 or min_lon > 180.0:
                continue
            if max_lat < -90.0 or min_lat > 90.0:
                continue

            lon_start = int(np.floor((max(min_lon, -180.0) + 180.0) / bin_size_degrees))
            lon_stop = int(np.floor((min(max_lon, 180.0) + 180.0) / bin_size_degrees))
            lat_start = int(np.floor((max(min_lat, -90.0) + 90.0) / bin_size_degrees))
            lat_stop = int(np.floor((min(max_lat, 90.0) + 90.0) / bin_size_degrees))
            lon_start = int(np.clip(lon_start, 0, longitude_bin_count - 1))
            lon_stop = int(np.clip(lon_stop, 0, longitude_bin_count - 1))
            lat_start = int(np.clip(lat_start, 0, latitude_bin_count - 1))
            lat_stop = int(np.clip(lat_stop, 0, latitude_bin_count - 1))

            path_index = len(indexed_subpaths)
            indexed_subpaths.append(
                (int(code), subpath, min_lon, max_lon, min_lat, max_lat)
            )
            for latitude_bin in range(lat_start, lat_stop + 1):
                bin_offset = latitude_bin * longitude_bin_count
                for longitude_bin in range(lon_start, lon_stop + 1):
                    paths_by_bin[bin_offset + longitude_bin].append(path_index)

    valid_queries = (
        np.isfinite(points).all(axis=1)
        & (points[:, 0] >= -180.0)
        & (points[:, 0] <= 180.0)
        & (points[:, 1] >= -90.0)
        & (points[:, 1] <= 90.0)
    )
    valid_indices = np.flatnonzero(valid_queries)
    if valid_indices.size:
        longitude_bins = np.floor(
            (points[valid_indices, 0] + 180.0) / bin_size_degrees
        ).astype(int)
        latitude_bins = np.floor(
            (points[valid_indices, 1] + 90.0) / bin_size_degrees
        ).astype(int)
        longitude_bins = np.clip(longitude_bins, 0, longitude_bin_count - 1)
        latitude_bins = np.clip(latitude_bins, 0, latitude_bin_count - 1)
        query_bins = latitude_bins * longitude_bin_count + longitude_bins
        sort_order = np.argsort(query_bins, kind="stable")
        sorted_indices = valid_indices[sort_order]
        sorted_bins = query_bins[sort_order]
        unique_bins, bin_starts, bin_counts = np.unique(
            sorted_bins, return_index=True, return_counts=True
        )

        for bin_id, bin_start, bin_count in zip(
            unique_bins, bin_starts, bin_counts, strict=True
        ):
            query_indices = sorted_indices[bin_start : bin_start + bin_count]
            bin_points = points[query_indices]
            for path_index in paths_by_bin[int(bin_id)]:
                code, subpath, min_lon, max_lon, min_lat, max_lat = indexed_subpaths[
                    path_index
                ]
                bbox_mask = (
                    (bin_points[:, 0] >= min_lon)
                    & (bin_points[:, 0] <= max_lon)
                    & (bin_points[:, 1] >= min_lat)
                    & (bin_points[:, 1] <= max_lat)
                )
                if not np.any(bbox_mask):
                    continue
                candidate_indices = query_indices[bbox_mask]
                inside = subpath.contains_points(points[candidate_indices])
                region_codes[candidate_indices[inside]] = code

    if np.isscalar(lon) and np.isscalar(lat):
        return region_codes[0]
    return region_codes.reshape(lon_arr.shape)


def _iter_single_spherical_polygons(polygon):
    if isinstance(polygon, SphericalPolygon):
        for child in polygon.polygons:
            yield from _iter_single_spherical_polygons(child)
    elif isinstance(polygon, SingleSphericalPolygon):
        yield polygon
    else:
        raise TypeError("Expected SphericalPolygon or SingleSphericalPolygon.")


def save_region_spherical_polygons(
    polygons_by_region: dict[int, list[list[SphericalPolygon]]],
    file_path: str | Path | None = None,
) -> Path:
    """Save region spherical polygon boundaries in a compact numeric NPZ.

    Components and their boundary rings are preserved. The file contains only
    numeric arrays and can be loaded without pickle.
    """
    output_path = (
        Path(file_path)
        if file_path is not None
        else _data_directory / "region_spherical_polygons.npz"
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)

    boundary_region_codes: list[int] = []
    boundary_component_ids: list[int] = []
    vertex_arrays: list[npt.NDArray[np.float64]] = []
    inside_points: list[npt.NDArray[np.float64]] = []
    vertex_offsets = [0]
    for region_code, components in sorted(polygons_by_region.items()):
        for component_id, component in enumerate(components):
            for polygon in component:
                for boundary in _iter_single_spherical_polygons(polygon):
                    vertices = np.asarray(boundary.points, dtype=float)
                    if vertices.ndim != 2 or vertices.shape[1] != 3:
                        raise ValueError(
                            "Spherical polygon vertices must be Nx3 vectors."
                        )
                    vertex_arrays.append(vertices)
                    vertex_offsets.append(vertex_offsets[-1] + vertices.shape[0])
                    inside_points.append(np.asarray(boundary.inside, dtype=float))
                    boundary_region_codes.append(int(region_code))
                    boundary_component_ids.append(component_id)

    if not boundary_region_codes:
        raise ValueError("At least one spherical polygon boundary is required.")

    np.savez_compressed(
        output_path,
        schema_version=np.asarray(1, dtype=np.uint8),
        region_codes=np.asarray(sorted(polygons_by_region), dtype=np.int32),
        boundary_region_codes=np.asarray(boundary_region_codes, dtype=np.int32),
        boundary_component_ids=np.asarray(boundary_component_ids, dtype=np.int32),
        vertex_offsets=np.asarray(vertex_offsets, dtype=np.int64),
        vertices_xyz=np.concatenate(vertex_arrays, axis=0),
        inside_xyz=np.stack(inside_points),
    )
    return output_path


def load_region_spherical_polygons(
    file_path: str | Path | None = None,
) -> dict[int, list[list[SingleSphericalPolygon]]]:
    """Load saved spherical region polygons from a numeric NPZ file."""
    npz_path = (
        Path(file_path)
        if file_path is not None
        else _data_directory / "region_spherical_polygons.npz"
    )
    with np.load(npz_path, allow_pickle=False) as data:
        schema_version = int(data["schema_version"])
        if schema_version != 1:
            raise ValueError(
                f"Unsupported spherical polygon schema version {schema_version}."
            )
        region_codes = np.asarray(data["region_codes"], dtype=int)
        boundary_region_codes = np.asarray(data["boundary_region_codes"], dtype=int)
        component_ids = np.asarray(data["boundary_component_ids"], dtype=int)
        vertex_offsets = np.asarray(data["vertex_offsets"], dtype=int)
        vertices_xyz = np.asarray(data["vertices_xyz"], dtype=float)
        inside_xyz = np.asarray(data["inside_xyz"], dtype=float)

    boundary_count = boundary_region_codes.size
    if (
        component_ids.shape != (boundary_count,)
        or inside_xyz.shape != (boundary_count, 3)
        or vertex_offsets.shape != (boundary_count + 1,)
        or vertex_offsets[0] != 0
        or vertex_offsets[-1] != vertices_xyz.shape[0]
    ):
        raise ValueError("Invalid spherical polygon NPZ structure.")

    grouped_boundaries: dict[int, dict[int, list[SingleSphericalPolygon]]] = {
        int(code): {} for code in region_codes
    }
    for index, (region_code, component_id) in enumerate(
        zip(boundary_region_codes, component_ids, strict=True)
    ):
        vertices = vertices_xyz[vertex_offsets[index] : vertex_offsets[index + 1]]
        center_vector = inside_xyz[index]
        longitudes = np.rad2deg(np.arctan2(vertices[:, 1], vertices[:, 0]))
        latitudes = np.rad2deg(np.arcsin(np.clip(vertices[:, 2], -1.0, 1.0)))
        center = (
            np.rad2deg(np.arctan2(center_vector[1], center_vector[0])),
            np.rad2deg(
                np.arctan2(
                    center_vector[2], np.hypot(center_vector[0], center_vector[1])
                )
            ),
        )
        boundary = SingleSphericalPolygon.from_lonlat(
            longitudes, latitudes, center=center
        )
        grouped_boundaries[int(region_code)].setdefault(int(component_id), []).append(
            boundary
        )

    return {
        region_code: [
            grouped_boundaries[region_code][component_id]
            for component_id in sorted(grouped_boundaries[region_code])
        ]
        for region_code in sorted(grouped_boundaries)
    }


def _vectors_to_lonlat(
    vectors: npt.NDArray[np.float64],
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    longitudes = np.rad2deg(np.arctan2(vectors[..., 1], vectors[..., 0]))
    latitudes = np.rad2deg(np.arcsin(np.clip(vectors[..., 2], -1.0, 1.0)))
    return longitudes, latitudes


def _reduce_shared_ring_vertices(
    rings: list[npt.NDArray[np.float64]], fraction: float
) -> list[npt.NDArray[np.float64]]:
    """Simplify rings together so vertices shared across rings stay synchronized."""
    previous: list[list[int]] = []
    following: list[list[int]] = []
    alive: list[npt.NDArray[np.bool_]] = []
    remaining: list[int] = []
    targets: list[int] = []
    occurrences: dict[tuple[float, ...], list[tuple[int, int]]] = {}
    canonical_vectors: dict[tuple[float, ...], npt.NDArray[np.float64]] = {}

    for ring_index, ring in enumerate(rings):
        count = ring.shape[0]
        previous.append([(index - 1) % count for index in range(count)])
        following.append([(index + 1) % count for index in range(count)])
        alive.append(np.ones(count, dtype=bool))
        remaining.append(count)
        targets.append(max(3, int(np.ceil(count * (1.0 - fraction)))))
        for vertex_index, vector in enumerate(ring):
            key = tuple(np.round(vector, decimals=12))
            occurrences.setdefault(key, []).append((ring_index, vertex_index))
            canonical_vectors.setdefault(key, vector)

    versions = dict.fromkeys(occurrences, 0)

    def live_occurrences(key: tuple[float, ...]) -> list[tuple[int, int]]:
        return [
            (ring_index, vertex_index)
            for ring_index, vertex_index in occurrences[key]
            if alive[ring_index][vertex_index]
        ]

    def triangle_area(ring_index: int, vertex_index: int) -> float:
        before = previous[ring_index][vertex_index]
        after = following[ring_index][vertex_index]
        a = rings[ring_index][before]
        b = rings[ring_index][vertex_index]
        c = rings[ring_index][after]
        return 0.5 * float(np.linalg.norm(np.cross(b - a, c - a)))

    def can_remove(key: tuple[float, ...]) -> bool:
        live = live_occurrences(key)
        return bool(live) and all(
            remaining[ring_index] > targets[ring_index] for ring_index, _ in live
        )

    def importance(key: tuple[float, ...]) -> float:
        return max(
            triangle_area(ring_index, vertex_index)
            for ring_index, vertex_index in live_occurrences(key)
        )

    heap = [
        (importance(key), key, versions[key]) for key in occurrences if can_remove(key)
    ]
    heapq.heapify(heap)

    while heap:
        _, key, stamp = heapq.heappop(heap)
        if stamp != versions[key] or not can_remove(key):
            continue

        affected_keys: set[tuple[float, ...]] = set()
        for ring_index, vertex_index in occurrences[key]:
            if not alive[ring_index][vertex_index]:
                continue
            before = previous[ring_index][vertex_index]
            after = following[ring_index][vertex_index]
            affected_keys.add(tuple(np.round(rings[ring_index][before], 12)))
            affected_keys.add(tuple(np.round(rings[ring_index][after], 12)))
            alive[ring_index][vertex_index] = False
            remaining[ring_index] -= 1
            following[ring_index][before] = after
            previous[ring_index][after] = before

        for affected_key in affected_keys:
            if can_remove(affected_key):
                versions[affected_key] += 1
                heapq.heappush(
                    heap,
                    (
                        importance(affected_key),
                        affected_key,
                        versions[affected_key],
                    ),
                )

    reduced_rings = []
    for ring_index, ring in enumerate(rings):
        kept_indices = np.flatnonzero(alive[ring_index])
        kept = ring[kept_indices].copy()
        for output_index, input_index in enumerate(kept_indices):
            key = tuple(np.round(ring[input_index], decimals=12))
            kept[output_index] = canonical_vectors[key]
        reduced_rings.append(kept)
    return reduced_rings


def reduce_region_polygon_points(
    polygons_by_region: dict[int, list[list[SphericalPolygon]]],
    fraction: float,
) -> dict[int, list[list[SphericalPolygon]]]:
    """Reduce polygon boundary vertices toward the requested fraction.

    Each boundary keeps at least 3 vertices. Vertices shared between boundaries
    are removed together, preserving shared edges.

    Parameters
    ----------
    polygons_by_region : dict[int, list[list[SphericalPolygon]]]
        Region code -> components -> polygons.
    fraction : float
        Target fraction of vertices to remove, in ``[0, 1)``. Shared edges can
        constrain a boundary to retain more vertices than this target.

    Returns
    -------
    dict[int, list[list[SphericalPolygon]]]
        Reduced polygons with the same nesting as the input.
    """
    if not 0.0 <= fraction < 1.0:
        raise ValueError("fraction must be in the interval [0, 1).")

    rings: list[npt.NDArray[np.float64]] = []
    centers: list[tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]] = []
    polygon_ring_indices: dict[int, list[list[list[int]]]] = {}
    for region_code, components in polygons_by_region.items():
        component_ring_indices = []
        for component in components:
            polygon_indices = []
            for polygon in component:
                boundary_indices = []
                for boundary in _iter_single_spherical_polygons(polygon):
                    vertices = np.asarray(boundary.points, dtype=float)
                    if np.allclose(vertices[0], vertices[-1]):
                        vertices = vertices[:-1]
                    boundary_indices.append(len(rings))
                    rings.append(vertices)
                    centers.append(
                        _vectors_to_lonlat(np.asarray(boundary.inside, dtype=float))
                    )
                polygon_indices.append(boundary_indices)
            component_ring_indices.append(polygon_indices)
        polygon_ring_indices[region_code] = component_ring_indices

    reduced_rings = _reduce_shared_ring_vertices(rings, fraction)
    reduced: dict[int, list[list[SphericalPolygon]]] = {}
    for region_code, components in polygon_ring_indices.items():
        reduced_components = []
        for component in components:
            reduced_polygons = []
            for boundary_indices in component:
                boundaries = []
                for ring_index in boundary_indices:
                    longitudes, latitudes = _vectors_to_lonlat(
                        reduced_rings[ring_index]
                    )
                    boundaries.append(
                        SingleSphericalPolygon.from_lonlat(
                            longitudes, latitudes, center=centers[ring_index]
                        )
                    )
                reduced_polygons.append(SphericalPolygon(boundaries))
            reduced_components.append(reduced_polygons)
        reduced[region_code] = reduced_components
    return reduced


def spherical_polygons_to_region_code(
    lon: float | npt.ArrayLike,
    lat: float | npt.ArrayLike,
    polygons_by_region: dict[int, list[list[SphericalPolygon]]],
    candidate_region_codes: npt.ArrayLike | None = None,
) -> int | None | npt.NDArray[np.float64]:
    """Find the region containing each lon/lat using spherical polygons.

    Components are unions of their boundary groups; rings within one group
    retain the converter's shared component semantics. If multiple regions
    contain a point on a shared boundary, a containing grid hint wins; without
    a hint, the smallest matching region code wins. ``candidate_region_codes``
    optionally supplies a grid-derived first guess for each point; spherical
    containment still verifies the result, and other regions are checked if
    the hint does not contain the point.
    Missing matches are returned as ``None`` for scalar queries and ``NaN``
    for array queries.
    """
    scalar_query = np.isscalar(lon) and np.isscalar(lat)
    lon_array, lat_array = np.broadcast_arrays(
        np.asarray(lon, dtype=float), np.asarray(lat, dtype=float)
    )
    finite_latitudes = lat_array[np.isfinite(lat_array)]
    if np.any((finite_latitudes < -90.0) | (finite_latitudes > 90.0)):
        raise ValueError("Latitude must be between -90 and 90 degrees.")
    candidate_array = None
    if candidate_region_codes is not None:
        try:
            candidate_array = np.broadcast_to(
                np.asarray(candidate_region_codes), lon_array.shape
            )
        except ValueError as error:
            raise ValueError(
                "candidate_region_codes must broadcast to the lon/lat shape."
            ) from error

    results = np.full(lon_array.shape, np.nan, dtype=float)
    ordered_regions = sorted(polygons_by_region)
    for flat_index, (longitude, latitude) in enumerate(
        zip(lon_array.flat, lat_array.flat, strict=True)
    ):
        if not np.isfinite(longitude) or not np.isfinite(latitude):
            continue
        longitude = (float(longitude) + 180.0) % 360.0 - 180.0
        search_order = ordered_regions
        if candidate_array is not None:
            candidate_code = candidate_array.flat[flat_index]
            if (
                np.isfinite(candidate_code)
                and int(candidate_code) in polygons_by_region
            ):
                candidate_code = int(candidate_code)
                search_order = [candidate_code] + [
                    code for code in ordered_regions if code != candidate_code
                ]
        for region_code in search_order:
            components = polygons_by_region[region_code]
            if any(
                all(
                    polygon.contains_lonlat(longitude, float(latitude))
                    for polygon in component
                )
                for component in components
            ):
                results.flat[flat_index] = region_code
                break

    if scalar_query:
        value = results.item()
        return None if np.isnan(value) else int(value)
    return results
