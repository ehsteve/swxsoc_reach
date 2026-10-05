import matplotlib.path as mpath
import numpy as np
import pytest
from spherical_geometry.polygon import SphericalPolygon

import swxsoc_reach.util.geom as geom_util
from swxsoc_reach.util.enums import Region


def test_load_region_contours_returns_all_region_codes(tmp_path):
    """load_region_contours should return paths for every supported region code."""
    contour_file = tmp_path / "region_contour_paths.npz"

    contour_levels = np.array(Region.contour_levels(), dtype=int)
    vertices_by_segment = np.array(
        [
            [[-1.0, -1.0], [-0.5, -1.0], [-0.5, -0.5], [-1.0, -1.0]],
            [[-2.0, -2.0], [-1.5, -2.0], [-1.5, -1.5], [-2.0, -2.0]],
            [[-3.0, -3.0], [-2.5, -3.0], [-2.5, -2.5], [-3.0, -3.0]],
            [[-4.0, -4.0], [-3.5, -4.0], [-3.5, -3.5], [-4.0, -4.0]],
            [[1.0, 1.0], [1.5, 1.0], [1.5, 1.5], [1.0, 1.0]],
            [[2.0, 2.0], [2.5, 2.0], [2.5, 2.5], [2.0, 2.0]],
            [[3.0, 3.0], [3.5, 3.0], [3.5, 3.5], [3.0, 3.0]],
            [[4.0, 4.0], [4.5, 4.0], [4.5, 4.5], [4.0, 4.0]],
        ],
        dtype=float,
    )

    all_vertices = vertices_by_segment.reshape(-1, 2)
    path_vertex_counts = np.full(contour_levels.size, 4, dtype=int)
    path_code_counts = np.full(contour_levels.size, 4, dtype=int)
    path_codes = np.array(
        [
            mpath.Path.MOVETO,
            mpath.Path.LINETO,
            mpath.Path.LINETO,
            mpath.Path.CLOSEPOLY,
        ],
        dtype=np.uint8,
    )
    all_codes = np.tile(path_codes, contour_levels.size)

    np.savez_compressed(
        contour_file,
        contour_levels=contour_levels,
        vertices=all_vertices,
        codes=all_codes,
        path_vertex_counts=path_vertex_counts,
        path_code_counts=path_code_counts,
    )

    paths_dict = geom_util.load_region_contours(contour_file=contour_file)
    assert set(paths_dict) == set(Region.contour_levels())
    assert all(isinstance(path, mpath.Path) for path in paths_dict.values())


def test_load_region_contours_rejects_legacy_object_npz(tmp_path):
    """Legacy NPZ schema should be rejected (no compatibility fallback)."""
    contour_file = tmp_path / "region_contour_paths_legacy.npz"
    contour_levels = np.array([1], dtype=int)
    vertices_by_segment = np.array(
        [np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 0.0]])],
        dtype=object,
    )
    np.savez_compressed(
        contour_file,
        contour_levels=contour_levels,
        vertices_by_segment=vertices_by_segment,
    )

    with pytest.raises(KeyError, match="vertices"):
        geom_util.load_region_contours(contour_file=contour_file)


def test_spherical_polygon_save_load_and_query(tmp_path):
    """Spherical polygon queries survive a numeric NPZ round trip."""
    left_polygon = SphericalPolygon.from_lonlat(
        [-2.0, 0.0, 0.0, -2.0],
        [-1.0, -1.0, 1.0, 1.0],
        center=(-1.0, 0.0),
    )
    right_polygon = SphericalPolygon.from_lonlat(
        [0.0, 2.0, 2.0, 0.0],
        [-1.0, -1.0, 1.0, 1.0],
        center=(1.0, 0.0),
    )
    polygons = {-1: [[left_polygon]], 1: [[right_polygon]]}
    polygon_file = tmp_path / "region_spherical_polygons.npz"

    geom_util.save_region_spherical_polygons(polygons, file_path=polygon_file)
    loaded_polygons = geom_util.load_region_spherical_polygons(file_path=polygon_file)

    assert geom_util.spherical_polygons_to_region_code(-1.0, 0.0, loaded_polygons) == -1
    assert geom_util.spherical_polygons_to_region_code(1.0, 0.0, loaded_polygons) == 1
    assert (
        geom_util.spherical_polygons_to_region_code(30.0, 0.0, loaded_polygons) is None
    )

    batch = geom_util.spherical_polygons_to_region_code(
        np.array([-1.0, 1.0, 30.0]),
        np.array([0.0, 0.0, 0.0]),
        loaded_polygons,
        candidate_region_codes=np.array([1, -1, -1]),
    )
    np.testing.assert_array_equal(batch[:2], [-1.0, 1.0])
    assert np.isnan(batch[2])


def test_spherical_polygon_query_has_deterministic_boundary_tie_break():
    """Overlapping boundary matches resolve to the lowest region code."""
    polygon = SphericalPolygon.from_lonlat(
        [-1.0, 1.0, 1.0, -1.0],
        [-1.0, -1.0, 1.0, 1.0],
        center=(0.0, 0.0),
    )

    result = geom_util.spherical_polygons_to_region_code(
        0.0,
        0.0,
        {1: [[polygon]], -1: [[polygon]]},
    )

    assert result == -1


def test_reduce_region_polygon_points_removes_fraction_of_vertices():
    """Reduction hits the requested vertex fraction and keeps the shape."""
    angles = np.linspace(0.0, 2.0 * np.pi, 40, endpoint=False)
    polygon = SphericalPolygon.from_lonlat(
        5.0 * np.cos(angles), 5.0 * np.sin(angles), center=(0.0, 0.0)
    )

    reduced = geom_util.reduce_region_polygon_points({1: [[polygon]]}, 0.5)

    boundary = reduced[1][0][0].polygons[0]
    assert len(boundary.points) - 1 == 20
    assert reduced[1][0][0].contains_lonlat(0.0, 0.0)
    assert not reduced[1][0][0].contains_lonlat(30.0, 0.0)
    with pytest.raises(ValueError, match="fraction"):
        geom_util.reduce_region_polygon_points({1: [[polygon]]}, 1.0)


def test_reduce_region_polygon_points_keeps_shared_boundary_vertices():
    """Adjacent polygons retain the same points along their shared boundary."""
    shared_latitudes = np.linspace(-10.0, 10.0, 11)
    left_lon = np.concatenate(
        ([0.0, -10.0, -10.0, 0.0], np.zeros(shared_latitudes.size - 2))
    )
    left_lat = np.concatenate(([-10.0, -10.0, 10.0, 10.0], shared_latitudes[-2:0:-1]))
    right_lon = np.concatenate(
        ([0.0], np.zeros(shared_latitudes.size - 1), [10.0, 10.0])
    )
    right_lat = np.concatenate((shared_latitudes, [10.0, -10.0]))
    adjacent = {
        -1: [[SphericalPolygon.from_lonlat(left_lon, left_lat, center=(-5.0, 0.0))]],
        1: [[SphericalPolygon.from_lonlat(right_lon, right_lat, center=(5.0, 0.0))]],
    }

    reduced = geom_util.reduce_region_polygon_points(adjacent, 0.5)

    shared_points_by_region = []
    for region_code in (-1, 1):
        boundary = reduced[region_code][0][0].polygons[0]
        longitudes, _ = geom_util._vectors_to_lonlat(boundary.points)
        shared_points_by_region.append(
            np.sort(np.round(longitudes[np.isclose(longitudes, 0.0)], 10))
        )
    np.testing.assert_array_equal(*shared_points_by_region)


def test_points_to_region_code_spatial_index_handles_compound_paths():
    """Indexed lookup preserves union semantics for compound contour paths."""
    left_path = mpath.Path(
        [(-10.0, -10.0), (0.0, -10.0), (0.0, 10.0), (-10.0, 10.0)],
        closed=True,
    )
    right_path = mpath.Path(
        [(0.0, -10.0), (10.0, -10.0), (10.0, 10.0), (0.0, 10.0)],
        closed=True,
    )
    compound_path = mpath.Path.make_compound_path(left_path, right_path)

    result = geom_util.points_to_region_code(
        np.array([-5.0, 5.0, 30.0]),
        np.array([0.0, 0.0, 0.0]),
        {7: compound_path},
    )

    np.testing.assert_array_equal(result[:2], [7.0, 7.0])
    assert np.isnan(result[2])
    assert geom_util.points_to_region_code(-5.0, 0.0, {7: left_path}) == 7
    assert np.isnan(geom_util.points_to_region_code(30.0, 0.0, {7: left_path}))
