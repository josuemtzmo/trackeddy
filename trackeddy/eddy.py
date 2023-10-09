import numpy as np
import pandas as pd

from trackeddy.decorators import eddy_identifier


class Eddy:
    """
    Eddy object to store data of identified and discarded eddies.

    Parameters
        ----------
        contour : array_like
            Values of closed contour of the eddy
        level : float
            Level at which the contour was identified
    """

    def __init__(self, contour, level) -> None:
        # Define all the object variables needed for the identification
        # of an eddy.
        self.contour = contour
        self.level = np.array(level)
        self.area_eddy = 0
        self.radius_eddy = 0
        self.ellipse = None
        self.ellipse_params = None
        self.contour_ellipse_error = 0
        self.ellipse_eccen = 0
        self.contour_center = None
        self.eddy_sign = 0
        self.eddy_maxima = 0
        self.gaussian_params = None
        self.contour_gaussian_error = 0

        # TODO remove these.
        self.contour_mask = None
        self.gaussian = None
        self.data_near_contour = None

    def to_dict(self, out_properties=None) -> dict:
        """
        to_dict convert contents of an eddy object into a dictionary

        Parameters
        ----------
        out_properties : list, optional
            List of variables to output, by default ('None') will output:
                -'identifier'
                -'contour'
                -'level'
                -'area_eddy'
                -'radius_eddy'
                -'ellipse_params'
                -'contour_ellipse_error'
                -'ellipse_eccen'
                -'contour_center'
                -'eddy_sign'
                -'eddy_maxima'
                -'gaussian_params'
                -'contour_gaussian_error'
        Returns
        -------
        dict
            Dictionary including all the eddy properties defined in
            the out_properties list.

        Raises
        ------
        ValueError
            If property does not exist within the eddy object, raise error.
        """
        if not out_properties:
            out_properties = [
                "identifier",
                "contour",
                "level",
                "area_eddy",
                "radius_eddy",
                "ellipse_params",
                "contour_ellipse_error",
                "ellipse_eccen",
                "contour_center",
                "eddy_sign",
                "eddy_maxima",
                "gaussian_params",
                "contour_gaussian_error",
            ]

        eddy_dict = {}
        for property in out_properties:
            attr = getattr(self, property)

            if not isinstance(attr, np.ndarray):
                attr = np.array(attr)

            attr = attr.squeeze()
            if len(attr.shape) > 2:
                pass
                # TODO Ravel data to export matrixes.
                # eddy_dict['contour_x']=attr[:,0]
            elif len(attr.shape) == 2 and 2 in attr.shape:
                eddy_dict["contour_path_x"] = attr[:, 0]
                eddy_dict["contour_path_y"] = attr[:, 1]
            elif not attr.shape:
                eddy_dict[property] = attr
            else:
                match property:
                    case "ellipse_params":
                        eddy_dict[property + "_a"] = attr[0]
                        eddy_dict[property + "_b"] = attr[1]
                        eddy_dict[property + "_theta"] = attr[2]
                    case "contour_center":
                        eddy_dict["contour_x"] = attr[0]
                        eddy_dict["contour_y"] = attr[1]
                    case "eddy_maxima":
                        eddy_dict["maxima"] = attr[0]
                        eddy_dict["maxima_y"] = attr[1]
                        eddy_dict["maxima_x"] = attr[2]
                    case "gaussian_params":
                        eddy_dict[property + "_x"] = attr[0]
                        eddy_dict[property + "_y"] = attr[1]
                        eddy_dict[property + "_sigma_x"] = attr[2]
                        eddy_dict[property + "_sigma_y"] = attr[3]
                        eddy_dict[property + "_theta"] = attr[4]
                    case _:
                        raise ValueError(
                            """{0} doesn't exist, make sure the property is a
                                attribute of the eddy class.""".format(
                                property
                            )
                        )

        return eddy_dict

    def to_table(self, out_properties=None, **kwards) -> pd.DataFrame:
        """
        to_table Output eddy properties to a pandas dataframe.

        Parameters
        ----------
        out_properties : list, optional
            see Eddy.to_dict for more info, by default None

        Returns
        -------
        pd.DataFrame
            Pandas table including all the eddy properties defined
            in the out_properties list.
        """
        eddy_dict = self.to_dict(out_properties)
        table = pd.DataFrame.from_dict(eddy_dict, **kwards).reset_index()
        return table

    def diagnose_eddy(self):
        # Print or plot the eddy
        pass

    @eddy_identifier
    def store(self, id):
        """
        store Store identified eddies into a pandas table

        Parameters
        ----------
        id : int
            Index of the identified eddy, provided by the decorator eddy_identifier

        Returns
        -------
        pd.DataFrame
            Table containing the information of a single eddy.
        """
        self.identifier = id
        table = self.to_table()
        return table

    @eddy_identifier
    def discarded(self, reason, id):
        """
        discarded Store discarded eddies into a simplified pandas table

        Parameters
        ----------
        reason : str
            String of why the eddy was discarded.
        id : int
            Index of the identified eddy, provided by the decorator eddy_identifier

        Returns
        -------
        pd.DataFrame
            Pandas table including all the eddy properties defined
            in the out_properties list.
        """
        self.reason = np.array(reason)
        out_properties = [
            "reason",
            "level",
            "area_eddy",
            "radius_eddy",
            "contour_ellipse_error",
            "contour_center",
            "ellipse_eccen",
            "eddy_sign",
            "eddy_maxima",
            "contour_gaussian_error",
        ]
        # TODO can I do something to use the to_table?
        # table = self.to_table(out_properties)
        dict = self.to_dict(out_properties)
        table = pd.DataFrame(dict, index=[id])
        return table

    def exit(self):
        """
        exit Reset counter and table in both store and discarded functions.
        """
        # Reset decorator counter in both functions.
        self.store(exit=True)
        self.discarded(exit=True)
