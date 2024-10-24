import os

import numpy as np
import pandas as pd
from scipy import ndimage
from sklearn.neighbors import BallTree


def _rename_eddies_in_time(current_time, index, nearest, prev_count):
    for near in nearest:
        shifted_counter = near + prev_count + 1

        current_time = current_time.rename(
            index={shifted_counter: index[near][0]}, level=0
        )

        # TODO add conditions to better estimate a track, i.e.
        # compute propagation speed and check distance.

    return current_time


def _update_counter(p_eddy, n_eddy):
    previous_index = p_eddy.index.levels[0]
    current_index = n_eddy.index.levels[0]
    # Maximum index in the previous table
    prev_count = previous_index.max()
    # Shift the current index by the largest unique identifier + 1 in
    # the previous eddy table. As if all the new eddies are unique.
    new_identifier = current_index + prev_count + 1
    new_index = n_eddy.index.set_levels(new_identifier, level=0)
    n_eddy.index = new_index
    return n_eddy, prev_count


def _update_better_eddy(p_eddy, n_eddy, index, nearest):
    # Get unique identifier for each eddy identified at the previous
    # and current level
    n_eddy, prev_count = _update_counter(p_eddy, n_eddy)

    replaced_eddies = []
    for duplicated in nearest:
        shifted_counter = duplicated + prev_count + 1
        previous_eddy = p_eddy.loc[index[duplicated][0]]
        current_eddy = n_eddy.loc[shifted_counter]
        previous_eddy = p_eddy.xs(
            (index[duplicated][0], 0), level=["identifier", "index"]
        )
        current_eddy = n_eddy.xs((shifted_counter, 0), level=["identifier", "index"])

        p_ellipse_error = previous_eddy.contour_ellipse_error.values
        c_ellipse_error = current_eddy.contour_ellipse_error.values

        p_gauss_error = previous_eddy.contour_gaussian_error.values
        c_gauss_error = current_eddy.contour_gaussian_error.values

        # TODO check if this is the best condition.
        if p_ellipse_error > c_ellipse_error and p_gauss_error > c_gauss_error:
            n_eddy = n_eddy.rename(
                index={shifted_counter: index[duplicated][0]}, level=0
            )
            replaced_eddies.append(index[duplicated][0])
        else:
            # Delete the input that does a worst job compared to
            # the previous level
            n_eddy = n_eddy.drop(shifted_counter, level=0)

    if replaced_eddies:
        p_eddy = p_eddy.drop(replaced_eddies, level=0)

    return p_eddy, n_eddy


def _merge_reindex_eddies(p_eddy, n_eddy):
    # Merge eddy tables
    new_eddy_table = pd.concat([p_eddy, n_eddy])
    # Get new length of unique identifiers
    new_eddy_count = len(new_eddy_table.index.levels[0])

    # Create new monotonously increasing index
    new_identifier = np.arange(0, new_eddy_count, dtype=int)

    # Reset table to avoid issues with codes.
    reset_table = new_eddy_table.reset_index().set_index(["identifier", "index"])

    # Create new index for the table
    new_index = reset_table.index.set_levels(new_identifier, level=0)
    # Assign the new index to the table
    reset_table.index = new_index

    # Important that index is sorted, so when computing distances,
    # the values will be assigned to an monotonously increasing index.
    return reset_table.sort_index(level=["identifier", "index"])


def _detect_nearest(p_eddy, n_eddy):
    # Extract from table the first element in the index 'index'
    eddy_info_p_level = p_eddy.xs(0, level=("index"))
    eddy_info_n_level = n_eddy.xs(0, level=("index"))
    # Extract coordinates of maxima within the contour.
    eddy_coords_p = eddy_info_p_level[["maxima_y", "maxima_x"]].values
    eddy_coords_n = eddy_info_n_level[["maxima_y", "maxima_x"]].values
    # Extract eddy radius of the new eddies.
    eddy_radius_n = eddy_info_n_level["radius_eddy"].values
    # Search nearest points between the previous eddies and the new eddies.
    # Important to pass ((lat, lot)).
    ball_tree = BallTree(np.deg2rad(eddy_coords_p), metric="haversine")
    distance, index = ball_tree.query(np.deg2rad(eddy_coords_n))
    # If distance is smaller than radius of eddy, then mark them as
    # within the eddy radius.
    distance_between_nearest = (distance * 6371) - np.expand_dims(eddy_radius_n, 1)
    # If the distance is negative, then the current eddy is within
    # the previous eddy.
    nearest, _ = np.where(distance_between_nearest < 0)
    return index, nearest


def _mask_data_in_contour(coord_ind_contour, data_near_contour, level, mask="contour"):
    # TODO Define if eddy is cyclonic or anticyclonic depending on
    # the values inside the contour or not.

    inside_contour = np.zeros(data_near_contour.shape)
    inside_contour[coord_ind_contour[:, 0], coord_ind_contour[:, 1]] = 1
    inside_contour = ndimage.binary_fill_holes(inside_contour)

    data_inside_contour = data_near_contour.where(inside_contour, np.nan).load()

    # TODO is there a faster way to check if the sign of
    # the contour interior.
    mean_value_inside_contour = data_inside_contour.mean()

    # This is used to differentiate between cyclonic and anticyclonic eddies
    eddy_sign = np.sign(mean_value_inside_contour).values
    # This computation allow to extract the eddy maxima independent
    # to their sign.
    eddy_maxima = eddy_sign * (eddy_sign * data_inside_contour).max()

    eddy_argmaxima = (eddy_sign * data_inside_contour).argmax()
    eddy_idxmaxima = np.unravel_index(eddy_argmaxima, data_inside_contour.shape)

    eddy_xlocmax = eddy_idxmaxima[1]
    eddy_ylocmax = eddy_idxmaxima[0]

    eddy_maxima_loc = np.vstack((eddy_maxima, eddy_ylocmax, eddy_xlocmax))

    if mask == "forcefit":
        # TODO reimplement the forcefit option, make a masked ring and
        # then set to zero.
        pass

    # After some testing, leaving the nans outside the contour in
    # data_inside_contour allows for better fitting of the gaussian
    # and eddy detection.
    return data_inside_contour, eddy_sign, eddy_maxima_loc, inside_contour


def _get_non_duplicated_identifiers(merged_identifiers):
    values, counts = np.unique(merged_identifiers, return_counts=True)

    duplicated = []
    for count in range(0, len(counts)):
        if counts[count] >= 2:
            duplicated.append(count)

    values = np.delete(values, duplicated)
    return values


def unlink_eddies_in_previous_times(previous_time, time, ntimes=5):
    # Move to decorator
    if time < ntimes:
        ntimes = time
    # Extract identifier of eddies identified as the last imput in the table.
    last_identified_eddies = previous_time.xs(time, level="time").index.levels[0]

    # Loop and arrays to contain extracted identifier of eddies identified
    # in the last ntimes in the table.
    unlinked_tracks_in_time = last_identified_eddies
    unlinked_times = np.ones_like(last_identified_eddies) * time
    # Loop back in time to check if missing a eddy track
    for prev_time in np.arange(time - ntimes, time):
        # Extracted identifier of eddies identified in the last time-ntimes
        # in the table.
        previous_identified_eddies = previous_time.xs(
            prev_time, level="time"
        ).index.levels[0]

        merged_identifiers = np.hstack(
            (last_identified_eddies, previous_identified_eddies)
        )

        values = _get_non_duplicated_identifiers(merged_identifiers)

        unlinked_tracks_in_time = np.hstack((unlinked_tracks_in_time, values))
        unlinked_times = np.hstack((unlinked_times, np.ones_like(values) * prev_time))

    # Perhaps not the cleanest code, but it ensures that only the last time
    # for all the detected eddies is the output.
    unlinked_tracks_in_time = np.unique(unlinked_tracks_in_time)

    times = (
        previous_time.loc[unlinked_tracks_in_time]
        .reset_index()
        .groupby("identifier")
        .max("time")["time"]
        .values
    )

    eddies_to_track = (
        previous_time.loc[unlinked_tracks_in_time]
        .groupby(["identifier", "index"])
        .last("time")
    )

    eddies_to_track["time"] = [times[index] for index, n in eddies_to_track.index]

    eddies_to_track = eddies_to_track.reset_index().set_index(
        ["identifier", "time", "index"]
    )

    previous_time = eddies_to_track.reindex(
        previous_time.index.levels[0], level=0, fill_value=0
    )

    return eddies_to_track


def track_online(track_in_time, current_time, time, ntimes):
    # Extracts the previous time and all the eddies that didn't
    # find a track in the previous 5 time steps.
    previous_time = unlink_eddies_in_previous_times(
        track_in_time, time - 1, ntimes=ntimes
    )

    # Detect nearest eddies.
    index, nearest = _detect_nearest(previous_time, current_time)

    # Update counter of current time to avoid overlap with the
    # previous time
    current_time, prev_count = _update_counter(previous_time, current_time)

    # Rename identifier of  nearest eddies to match the previous
    # table identifier
    current_time = _rename_eddies_in_time(current_time, index, nearest, prev_count)

    # Merge dictionaries.
    track_in_time = _merge_reindex_eddies(track_in_time, current_time)

    # Reset index and assign indexes for consistency in the loop
    # and output.
    track_in_time = (
        track_in_time.reset_index()
        .set_index(["identifier", "time", "index"])
        .sort_index(level=["identifier", "time", "index"])
    )
    return track_in_time


def store_snapshot(track_in_time, path, time, format="csv", **kwargs):

    out_path = os.path.join(path, "identified_eddies_{0:05}.{1}".format(time, format))

    match format:
        case "csv":
            track_in_time.to_csv(out_path, **kwargs)
        case "json":
            track_in_time.to_json(out_path, **kwargs)
        case "hd5":
            track_in_time.to_hdf(out_path, **kwargs)
