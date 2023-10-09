import numpy as np
import pandas as pd
import xarray as xr
from astropy import convolution
from scipy import ndimage
from sklearn.neighbors import BallTree

from trackeddy.decorators import check_return_eddies, filter_data
from trackeddy.eddy import Eddy
from trackeddy.geometry import (
    Fit_Surface,
    area_latlon_polygon,
    compute_similarity,
    eccentricity,
    extract_contours,
    fit_ellipse,
)
from trackeddy.physics import approx_RD_lat
from trackeddy.printfunc import Printer


class TrackEddy:
    def __init__(self, path, variable, nan_treatment=False, **xrkwargs) -> None:
        self.Dataset = xr.open_mfdataset(path, **xrkwargs)
        self.rawdata = self.Dataset[variable]
        self.nan_treatment = nan_treatment
        self.filter = "space"

        self.identify_coords()
        self.check_coords()

        self.identification_criteria = {
            "ellipse_fit": 0.85,
            "eccentricity": 0.85,
            "gaussian_fit": 0.8,
            "max_area": 2 * np.pi,
        }

        self.S_filter_setup = {
            "filter_type": "convolution",
            "mode": "uniform",
            "kernel": 10,
        }

        self.T_filter_setup = {}
        self.skip_gaussian_fit = False

    def check_coords(self):
        if "x_var" in self.coords and "y_var" in self.coords:
            x_coord_name = self.coords["x_var"]
            y_coord_name = self.coords["y_var"]
        elif self.coords["x"] and self.coords["y"]:
            x_coord_name = self.coords["x"]
            y_coord_name = self.coords["y"]
        else:
            raise ValueError(
                """Can't find dimensions or variables that correspond to
                geographic coordinages (lon,lat)"""
            )

        if (
            len(self.Dataset[x_coord_name].shape) == 1
            and len(self.Dataset[y_coord_name].shape) == 1
        ):
            self.X, self.Y = np.meshgrid(
                self.Dataset[x_coord_name], self.Dataset[y_coord_name]
            )

        elif (
            len(self.Dataset[x_coord_name].shape) == 2
            and len(self.Dataset[y_coord_name].shape) == 2
        ):
            self.X = self.Dataset[x_coord_name].values
            self.Y = self.Dataset[y_coord_name].values

        else:
            raise ValueError(
                """The coordinates of the dataset should be a 1D or 2D array,
                check the vars: '{0}', {1}.""".format(
                    x_coord_name, y_coord_name
                )
            )

    # TODO add a wrapper to check if the coords where properly identified,
    # if not rise error.
    def identify_coords(self) -> None:
        self.coords = {"time": "", "x": "", "y": ""}
        for dim_name in self.rawdata.dims:
            self._get_dims(dim_name)
        for var_name in self.Dataset.keys():
            self._get_vars(var_name)

    def _get_dims(self, dim_name):
        if "time" in dim_name or dim_name == "t":
            self.coords["time"] = dim_name
        elif "lon" in dim_name or dim_name == "x":
            self.coords["x"] = dim_name
        elif "lat" in dim_name or dim_name == "y":
            self.coords["y"] = dim_name

    def _get_vars(self, var_name):
        if "lon" in var_name or var_name == "x":
            self.coords["x_var"] = var_name
        elif "lat" in var_name or var_name == "y":
            self.coords["y_var"] = var_name

    def setup(self) -> None:
        pass

    @filter_data
    def _filter_data_(self, data_snapshot, filter=None) -> xr.DataArray:
        if filter == "time":
            data2track = self._time_filter(data_snapshot, **self.T_filter_setup)
        elif filter == "space":
            data2track = self._space_filter(data_snapshot, **self.S_filter_setup)
        else:
            raise ValueError(
                """Select a filter to extract field anomaly.
                The current options are: 'time' and 'space'."""
            )

        return data2track

    def _time_filter(self, filter):
        pass

    def _space_filter(
        self, data2track, filter_type, mode="uniform", kernel=10
    ) -> xr.DataArray:
        # Apply spatial filter.
        if filter_type == "convolution":
            if mode == "uniform":
                tmp_data = data2track.squeeze().copy()
                if kernel % 2 == 0:
                    ker = np.ones((kernel + 1, kernel + 1))
                else:
                    ker = np.ones((kernel, kernel))
                convolved_data = convolution.convolve(tmp_data, kernel=ker)
                tmp_data = tmp_data - convolved_data

                data2track = tmp_data
            if mode == "gaussian":
                raise Warning(
                    """ndimage.gaussian_filter may create artefacts near nan
                    values. Therefore, data is filled with zeros."""
                )
                tmp_data = data2track.squeeze().copy()
                tmp_data = tmp_data - ndimage.gaussian_filter(tmp_data, size=kernel)
                data2track = tmp_data, mask
        # Check if the user selects an moving average meridional filter.
        elif filter_type == "meridional":
            data2track = data2track - data2track.mean(self.coords["x"])
        # Check if the user selects an moving average zonal filter.
        elif filter_type == "zonal":
            data2track = data2track - data2track.mean(self.coords["y"])
        else:
            raise ValueError(
                """Define the filter_type argument between the options:
                'convolution', 'meridional', and 'zonal'"""
            )
        return data2track

    def time_tracking(self, t0=0, tf=None, lin_levels=None, ntimes=5):
        # TODO Move to decorator?
        if not tf:
            times = range(0, len(self.Dataset[self.coords["time"]]))
        else:
            times = range(t0, tf)

        track_in_time = pd.DataFrame({"": []})

        pp = Printer()

        for time in times:
            current_time = self._detect_snapshot(time, lin_levels)
            current_time["time"] = time
            if track_in_time.empty:
                track_in_time = (
                    current_time.reset_index()
                    .set_index(["identifier", "time", "index"])
                    .sort_index(level=["identifier", "time", "index"])
                )
                continue

            # Reindex to avoid loosing indexes, it's important for the
            # nearest detection
            current_time = current_time.reset_index().set_index(
                ["identifier", "time", "index"]
            )

            # TODO: Output to disk, instead of online tracking.

            # Extracts the previous time and all the eddies that didn't
            # find a track in the previous 5 time steps.
            previous_time = unlink_eddies_in_previous_times(
                track_in_time, time - 1, ntimes=ntimes
            )

            # Detect nearest eddies.
            index, nearest = self._detect_nearest(previous_time, current_time)

            # Update counter of current time to avoid overlap with the
            # previous time
            current_time, prev_count = self._update_counter(previous_time, current_time)

            # Rename identifier of  nearest eddies to match the previous
            # table identifier
            current_time = self._rename_eddies_in_time(
                current_time, index, nearest, prev_count
            )

            # Merge dictionaries.
            track_in_time = self._merge_reindex_eddies(track_in_time, current_time)

            # Reset index and assign indexes for consistency in the loop
            # and output.
            track_in_time = (
                track_in_time.reset_index()
                .set_index(["identifier", "time", "index"])
                .sort_index(level=["identifier", "time", "index"])
            )

            pp.timepercentprint(
                tf, 1, time, "# of E " + str(len(track_in_time.index.levels[0]))
            )
        return track_in_time

    def _detect_snapshot(self, time, levels):
        if isinstance(time, int) and self.coords["time"]:
            data_snapshot = self.rawdata.isel({self.coords["time"]: time}).squeeze()
        elif isinstance(time, str) and self.coords["time"]:
            data_snapshot = self.rawdata.isel({self.coords["time"]: time}).squeeze()
        else:
            data_snapshot = self.rawdata.squeeze()

        self._filter_data_(data_snapshot, filter=self.filter)

        joint_eddies = pd.DataFrame({"": []})

        for level in levels:
            eddies_current, discarded = self._detect_one_level(level)

            # If no eddies identified, then continue
            if eddies_current.empty:
                continue

            if joint_eddies.empty:
                joint_eddies = eddies_current
                # Continue, since no check needs to be done.
                continue

            # TODO pass the previous and current eddy to the
            # _detect_nearest function
            index, nearest = self._detect_nearest(joint_eddies, eddies_current)

            joint_eddies, eddies_current = self._update_better_eddy(
                joint_eddies, eddies_current, index, nearest
            )
            # TODO function with decorator that handles the merging of
            # eddy tracks.

            joint_eddies = self._merge_reindex_eddies(joint_eddies, eddies_current)

        return joint_eddies

    def _detect_nearest(self, p_eddy, n_eddy):
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

    def _rename_eddies_in_time(self, current_time, index, nearest, prev_count):
        for near in nearest:
            shifted_counter = near + prev_count + 1

            current_time = current_time.rename(
                index={shifted_counter: index[near][0]}, level=0
            )

            # TODO add conditions to better estimate a track, i.e.
            # compute propagation speed and check distance.

        return current_time

    def _update_counter(self, p_eddy, n_eddy):
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

    def _update_better_eddy(self, p_eddy, n_eddy, index, nearest):
        # Get unique identifier for each eddy identified at the previous
        # and current level
        n_eddy, prev_count = self._update_counter(p_eddy, n_eddy)

        replaced_eddies = []
        for duplicated in nearest:
            shifted_counter = duplicated + prev_count + 1
            previous_eddy = p_eddy.loc[index[duplicated][0]]
            current_eddy = n_eddy.loc[shifted_counter]
            previous_eddy = p_eddy.xs(
                (index[duplicated][0], 0), level=["identifier", "index"]
            )
            current_eddy = n_eddy.xs(
                (shifted_counter, 0), level=["identifier", "index"]
            )

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

    def _merge_reindex_eddies(self, p_eddy, n_eddy):
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

    def treat_nan(self, data_snapshot, nan_value=0):
        if not self.nan_treatment:
            data2track = data_snapshot.copy()
        elif self.nan_treatment:
            data2track = data_snapshot.fillna(nan_value).copy()
        else:
            raise ValueError("The nan_treatment can only be True or False")
        return data2track

    @check_return_eddies
    def _detect_one_level(self, level):
        # Define eddy_info and discarded variables
        eddy_info = None
        discarded = None
        eddy = None

        contours, _ = extract_contours(self.X, self.Y, self.data2track.values, level)

        for contour in contours:
            # Discard contour if contour is a singular point or if it
            # contains nans
            if np.allclose(np.mean(contour, axis=0), contour[0], equal_nan=True):
                continue

            eddy = Eddy(contour, level)

            # Brute force the removal of contours that are too large
            eddy.area_eddy = area_latlon_polygon(eddy.contour)
            eddy.radius_eddy = np.array(np.sqrt(eddy.area_eddy / np.pi))

            eddy.contour_center = np.nanmean(eddy.contour, 0)

            # Continue to next contour if eddy size is larger than the user
            # defined value multiplied by the Rossby deformation radius.
            Rd = approx_RD_lat(eddy.contour_center[1])
            max_size = self.identification_criteria["max_area"] * Rd
            if eddy.radius_eddy > max_size:
                discarded = eddy.discarded("area")
                continue

            # TODO clean the fit_ellipse function and input only eddy.contour
            eddy.ellipse, eddy.ellipse_params = fit_ellipse(
                eddy.contour[:, 0], eddy.contour[:, 1]
            )

            # Check if area is similar between fitted ellipse and contour
            # Discarding by area first improves the time.
            area_ellipse = area_latlon_polygon(eddy.ellipse)
            # The areas need to be at least within 20% of the value defined
            # by the user, default ~ 15%.
            if not np.isclose(
                eddy.area_eddy,
                area_ellipse,
                rtol=0.1 * self.identification_criteria["ellipse_fit"],
            ):
                # Continue to next contour if ellipse fit is bad
                discarded = eddy.discarded("ellipse_area")
                continue

            # Check eccentricity of ellipse.
            eddy.ellipse_eccen = eccentricity(*eddy.ellipse_params[0:2])
            if eddy.ellipse_eccen > self.identification_criteria["eccentricity"]:
                discarded = eddy.discarded(reason="eccentricity")
                continue

            eddy.contour_ellipse_error = compute_similarity(eddy.contour, eddy.ellipse)

            # Check if the similarity between the contours matches the
            # user criteria.
            if eddy.contour_ellipse_error > self.identification_criteria["ellipse_fit"]:
                discarded = eddy.discarded(reason="similarity")
                continue

            # TODO Extract the area from the contour.
            data_near_contour, x_near_c, y_near_c = self._data_in_contour(eddy)

            if eddy.eddy_maxima[0] == 0:
                discarded = eddy.discarded(reason="eddy_is_island")
                continue

            # Ignore gaussian fitting, if user changes parameter.
            if self.skip_gaussian_fit:
                eddy_info = eddy.store()
                continue

            # Create object to handle fitting of surfaces.
            F_surface = Fit_Surface(eddy, data_near_contour, x_near_c, y_near_c)

            # Fit the chosen feature.
            eddy.gaussian, eddy.gaussian_params = F_surface._fitting()

            # Extract the contour of the gaussian fitted + the eddy contour,
            # since all the gaussians decay to zero, for a fair comparison.
            # TODO check when eddy is negative in positive level, it may fail,
            # but they may be identified when tracking in negative levels.
            gauss_contour, _ = extract_contours(
                x_near_c, y_near_c, eddy.gaussian, eddy.level
            )

            # No contour was extracted. This is an issue with the gaussian
            # fitting optimization.
            if not gauss_contour:
                discarded = eddy.discarded(reason="gaussian_fit_failed")
                continue

            # TODO Move to a decorator of extract_contours.
            # Fix in case the mask cuts the contour
            if len(gauss_contour) > 1:
                gauss_contour = np.vstack(gauss_contour)
            else:
                gauss_contour = gauss_contour[0]

            # Compute similarity between the eddy contour and the gaussian
            # contour at the same level.
            eddy.contour_gaussian_error = compute_similarity(
                eddy.contour, gauss_contour
            )

            # If the similarity between contours doesn't match criteria
            # the eddy is descarde.
            if (
                eddy.contour_gaussian_error
                < self.identification_criteria["gaussian_fit"]
            ):
                discarded = eddy.discarded(reason="gaussian_check")
                continue

            # Copy identifier into the eddy object.
            eddy_info = eddy.store()

        return eddy_info, discarded, eddy

    def _data_in_contour(self, eddy):
        pt = np.expand_dims(np.mean(eddy.contour, axis=0), 0)

        # To speed up the process of finding the nearest point,
        # we mask values near the point that are close by 5% of its value,
        # with a largest cap of 3 degrees if coordinates are geographical.
        # The atol fixes the issue in x or y coordinates close to 0.
        Xmask = np.isclose(self.X, pt[0][0], rtol=0.05, atol=3)
        Ymask = np.isclose(self.Y, pt[0][1], rtol=0.05, atol=3)

        # This step finds some rough coordinates to crop the domain so
        # the BallTree search is quicker.
        sum_x = np.sum(Xmask * Ymask, 0)
        sim_y = np.sum(Xmask * Ymask, 1)
        # Corners of the box
        shift_x = np.argmax(sum_x)
        shift_y = np.argmax(sim_y)

        # Extract last item, since it will always be different than 0.
        mask_shape = np.unique(sum_x)[-1], np.unique(sim_y)[-1]

        x = self.X[Xmask * Ymask]
        y = self.Y[Xmask * Ymask]

        # Another option is to search through all the X and Y space,
        # but it's 100 times slower to compute.
        # TODO add an option to bruteforce in case it's needed,
        # but I can't see a case right now.
        # x = self.X.ravel()
        # y = self.Y.ravel()

        coords = np.vstack((x, y)).T

        # This allows support for non structured grids, and it's less
        # computationaly expensive than looking over the full coordinate set.
        ball_tree = BallTree(coords)
        dist, ind = ball_tree.query(eddy.contour, k=1)

        # Extract coordinates from the mask shape
        coord_slice = np.hstack(np.unravel_index(ind, mask_shape)).squeeze()

        # Shift back to the original coordinates
        coord_ind = coord_slice + [shift_y, shift_x]

        (
            eddy.data_near_contour,
            contour_coords_grid,
            x_near_c,
            y_near_c,
        ) = self._get_data_in_contour(coord_ind)

        # data_near_contour is loaded here, this should speed up
        # the gaussian fit, since we don't need to reload the data
        (
            masked_data_near_contour,
            eddy_sign,
            eddy_maxima,
            eddy_contour_mask,
        ) = self._mask_data_in_contour(
            contour_coords_grid, eddy.data_near_contour, eddy.level
        )

        # Support to extract the coordinates of maximum with and
        # without a regular grid
        Y_coord = int(eddy_maxima[1])
        X_coord = int(eddy_maxima[2])
        eddy_maxima[2] = x_near_c[Y_coord, X_coord]
        eddy_maxima[1] = y_near_c[Y_coord, X_coord]

        # TODO Continue here, it seems that changing this affects
        # the gaussian fitting and other things.

        eddy.eddy_sign = eddy_sign
        eddy.eddy_maxima = eddy_maxima
        eddy.contour_mask = eddy_contour_mask
        return masked_data_near_contour, x_near_c, y_near_c

    def _get_data_in_contour(self, coord_ind, threshold="auto"):
        # Get corners
        TR_corner = np.max(coord_ind, axis=0)
        BL_corner = np.min(coord_ind, axis=0)
        if threshold == "auto":
            diag_gridpoints = TR_corner - BL_corner
            thresh_value = int(0.5 * np.max(diag_gridpoints))
            if thresh_value <= 3:
                thresh_value = 3
        elif isinstance(threshold, int):
            thresh_value = threshold
        else:
            thresh_value = 0

        TR_corner = TR_corner + [thresh_value + 1, thresh_value + 1]
        BL_corner = BL_corner - [thresh_value, thresh_value]

        # Make sure that corners are never larger or smaller than
        # the dataset dimensions

        while (TR_corner > self.X.shape).any():
            index_max = np.argmax(TR_corner - self.X.shape)
            TR_corner[index_max] = TR_corner[index_max] - 1

        while (BL_corner < (0, 0)).any():
            index_min = np.argmin(BL_corner)
            BL_corner[index_min] = 0

        data_near_contour = self.data2track.isel(
            {
                self.coords["x"]: slice(BL_corner[1], TR_corner[1]),
                self.coords["y"]: slice(BL_corner[0], TR_corner[0]),
            }
        )

        x = self.X[BL_corner[0] : TR_corner[0], BL_corner[1] : TR_corner[1]]
        y = self.Y[BL_corner[0] : TR_corner[0], BL_corner[1] : TR_corner[1]]

        contour_coords_grid = coord_ind - BL_corner

        return data_near_contour, contour_coords_grid, x, y

    def _mask_data_in_contour(
        self, coord_ind_contour, data_near_contour, level, mask="contour"
    ):
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

    def plot_eddy_detection_multilevel(self, df_eddy_multilevel_store):
        import cmocean as cm
        import matplotlib.colors as mcolors
        import matplotlib.pyplot as plt

        n_colors = list(mcolors.XKCD_COLORS.keys())

        levels = df_eddy_multilevel_store.level.unique()
        color_index = np.linspace(0, len(n_colors) - 1, len(levels), dtype=int)

        plt.figure(figsize=(10, 6), dpi=300)
        plt.pcolormesh(
            self.X, self.Y, self.data2track, cmap=cm.cm.balance, vmin=-0.5, vmax=0.5
        )

        for eddy in df_eddy_multilevel_store.index.get_level_values(level=0):
            c_path_x = df_eddy_multilevel_store.loc[eddy].contour_path_x
            c_path_y = df_eddy_multilevel_store.loc[eddy].contour_path_y

            inx = np.where(df_eddy_multilevel_store.loc[eddy].level[0] == levels)[0][0]

            plt.plot(c_path_x, c_path_y, color=n_colors[color_index[inx]])

        if len(levels) < 5:
            for level in range(len(levels)):
                plt.plot(
                    c_path_x[0], c_path_y[0], color=n_colors[level], label=levels[level]
                )
            plt.legend()

    def plot_eddy_detection_in_level(self, eddies, discarded, **plot_args):
        if not plot_args:
            plot_args = {"alpha": 0.5, "markersize": 1}

        import cmocean as cm
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(1, 2, figsize=(10, 5), dpi=300)

        ax[0].pcolormesh(
            self.X,
            self.Y,
            self.data2track.squeeze().values,
            cmap=cm.cm.balance,
            vmin=-0.5,
            vmax=0.5,
        )

        for idx in eddies.index.levels[0]:
            Xcontour = eddies.loc[idx]["contour_path_x"]
            Ycontour = eddies.loc[idx]["contour_path_y"]

            x_cont = eddies.loc[idx]["contour_x"]
            y_cont = eddies.loc[idx]["contour_y"]

            ax[0].plot(Xcontour, Ycontour, "-m")
            ax[0].plot(x_cont, y_cont, ".k", markersize=1)

        ax[1].pcolormesh(
            self.X,
            self.Y,
            self.data2track.squeeze().values,
            cmap=cm.cm.balance,
            vmin=-0.5,
            vmax=0.5,
        )

        for idx in discarded.index:
            x_cont = discarded.loc[idx]["contour_x"]
            y_cont = discarded.loc[idx]["contour_y"]

            match discarded.loc[idx].reason:
                case "eccentricity":
                    ax[1].plot(x_cont, y_cont, ".r", **plot_args)
                case "area":
                    ax[1].plot(x_cont, y_cont, ".", color="gray", **plot_args)
                case "similarity":
                    ax[1].plot(x_cont, y_cont, ".c", **plot_args)
                case "ellipse_area":
                    ax[1].plot(x_cont, y_cont, ".m", **plot_args)
                case "eddy_is_island":
                    ax[1].plot(x_cont, y_cont, ".", color="orange", **plot_args)
                case "gaussian_fit_failed":
                    ax[1].plot(x_cont, y_cont, ".", color="pink", **plot_args)
                case _:
                    ax[1].plot(x_cont, y_cont, ".b", **plot_args)

        label = [
            "Fail eccentricity",
            "Area fail",
            "Contour != ellipse",
            "Ellipse area fail",
            "Eddy is island",
            "Gaussian fit failed",
            "Contour != Gaussian",
        ]

        colors = ["r", "gray", "c", "m", "orange", "pink", "b"]
        [
            ax[1].plot(
                self.X[0, 0],
                self.Y[0, 0],
                ".",
                color=colors[c],
                label=label[c],
                alpha=0.5,
            )
            for c in range(len(colors))
        ]

        ax[0].set_title("Identified")
        ax[1].set_title("Discarded")

        ax[1].legend(loc="upper center", bbox_to_anchor=(0, -0.05), ncol=3)


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
