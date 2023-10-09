import random as rnd

import numpy as np


def dist(loc1, loc2):
    return np.sqrt((loc1[0] - loc2[0]) ** 2 + (loc2[1] - loc1[1]) ** 2)


class Generate_field:
    def __init__(self, a, b, n, x, y, opt=""):
        # BUG WHEN LEN(x) != LEN(y) ?
        self.xlen = len(x)
        self.ylen = len(y)
        self.a = a * rnd.uniform(0.7, 1.3)
        self.b = b * rnd.uniform(0.7, 1.3)
        self.x = x
        self.y = y
        self.n = n
        self.opt = opt
        if type(self.n) != list or type(self.n) != tuple:
            self.eddies = {
                "eddy_n%s"
                % ii: {
                    "loc": [
                        [rnd.randint(0, self.xlen - 1), rnd.randint(0, self.ylen - 1)]
                    ],
                    "grow": True,
                    "radius": [self.a, self.b],
                    "angle": rnd.uniform(0, 2 * np.pi),
                    "amp": rnd.choice([-1, 1]) * rnd.uniform(0.7, 1.3),
                }
                for ii in range(self.n)
            }
        else:
            raise ValueError("No right input.")

    def go_right(self, indexs, step):
        return [0, step]

    def go_upright(self, indexs, step):
        return [step, step]

    def go_up(self, indexs, step):
        return [step, 0]

    def go_upleft(self, indexs, step):
        return [step, -step]

    def go_left(self, indexs, step):
        return [0, -step]

    def go_downleft(self, indexs, step):
        return [-step, -step]

    def go_down(self, indexs, step):
        return [-step, 0]

    def go_downright(self, indexs, step):
        return [-step, step]

    def twoD_Gaussian(
        self, coords, sigma_x, sigma_y, theta, slopex=0, slopey=0, offset=0
    ):
        """
        *************** twoD_Gaussian *******************
        Build a 2D gaussian.
        Notes:
            Remmember to do g.ravel().reshape(len(x),len(y)) for plotting purposes.
        Args:
            coords [x,y] (list|array): Coordinates in x and y.
            amplitude (float): Amplitud of gaussian.
            x0 , yo (float): Center of Gausian.
            sigma_x,sigma_y (float): Deviation.
            theta (Float): Orientation.
            offset (Float): Gaussian Offset.
        Returns:
            g.ravel() (list|array) - Gaussian surface in a list.
        Usage:
            Check scan_eddym function.
        """
        x = coords[0]
        y = coords[1]
        amplitude = coords[2]

        xo = float(coords[3])
        yo = float(coords[4])

        xo = float(xo)
        yo = float(yo)

        if sigma_y or sigma_x != 0:
            a = (np.cos(theta) ** 2) / (2 * sigma_x**2) + (np.sin(theta) ** 2) / (
                2 * sigma_y**2
            )
            b = -(np.sin(2 * theta)) / (4 * sigma_x**2) + (np.sin(2 * theta)) / (
                4 * sigma_y**2
            )
            c = (np.sin(theta) ** 2) / (2 * sigma_x**2) + (np.cos(theta) ** 2) / (
                2 * sigma_y**2
            )
            g = amplitude * np.exp(
                -(
                    a * ((x - xo) ** 2)
                    + 2 * b * (x - xo) * (y - yo)
                    + c * ((y - yo) ** 2)
                )
            )
        else:
            g = (x - xo) * 0 + (y - yo) * 0
        return g.ravel()

    def checkposition(self, away_val=5, loc=False):
        if loc:
            eddies_loc = [
                [rnd.randint(0, self.xlen - 1), rnd.randint(0, self.ylen - 1)]
                for key, item in self.eddies.items()
            ]
        else:
            eddies_loc = [item["loc"][-1] for key, item in self.eddies.items()]
        for key1, item1 in self.eddies.items():
            xc1 = item1["loc"][0][0]
            yc1 = item1["loc"][0][1]
            distance = np.array(
                [
                    dist([self.x[xc1], self.y[yc1]], [self.x[ii], self.y[jj]])
                    for ii, jj in eddies_loc
                ]
            )
            distance[distance == 0] = away_val * self.a
            checker = (
                (distance < away_val * self.a).any()
                or (distance < away_val * self.b).any()
            ) or loc
            count = 0
            while checker or count >= 10000:
                newx = rnd.randint(0, self.xlen - 1)
                newy = rnd.randint(0, self.ylen - 1)
                self.eddies[key1]["loc"] = [[newx, newy]]
                eddies_loc = [item["loc"][-1] for key, item in self.eddies.items()]
                # pdb.set_trace()
                xc1 = newx
                yc1 = newy
                distance = np.array(
                    [
                        dist([self.x[xc1], self.y[yc1]], [self.x[ii], self.y[jj]])
                        for ii, jj in eddies_loc
                    ]
                )
                numzeros = [ii for ii in distance if ii == 0]
                if len(numzeros) <= 1:
                    distance[distance == 0] = np.inf
                else:
                    distance[distance == 0] = away_val * self.a
                checker = (distance < away_val * self.a).any() or (
                    distance < away_val * self.b
                ).any()
                count = count + 1
        if loc:
            return self.eddies

    def make_random_walk(self, indexs, steps):
        move_dict = {
            1: self.go_up,
            2: self.go_right,
            3: self.go_left,
            4: self.go_down,
            5: self.go_downleft,
            6: self.go_downright,
            7: self.go_upleft,
            8: self.go_upright,
        }
        # for _ in range(steps):
        for ii in indexs:
            move_in_a_direction = move_dict[rnd.randint(1, 8)]
            movcood = move_in_a_direction(ii, steps)

        return indexs[0] + movcood[0], indexs[1] + movcood[1]

    def assemble_field(self, N, margin=50):
        data = np.zeros((N, self.xlen + 2 * margin, self.ylen + 2 * margin))
        for t in range(N):
            # pdb.set_trace()
            if self.opt == "no_interaction" or self.opt == "Nint":
                self.eddies = self.checkposition(away_val=5, loc=True)
            else:
                pass
            for keys, item in self.eddies.items():
                gauss = self.twoD_Gaussian(
                    self.pass_args(keys, margin),
                    item["radius"][0],
                    item["radius"][1],
                    item["angle"],
                ).reshape(np.shape(data[0, :, :]))
                data[t, :, :] = data[t, :, :] + gauss
        return data

    def reconstruct_field(self):
        data = np.zeros((self.xlen, self.ylen))
        for keys, item in self.eddies.items():
            gauss = self.twoD_Gaussian(
                self.pass_args(keys),
                item["radius"][0],
                item["radius"][1],
                item["angle"],
            ).reshape(np.shape(data))
            data = data + gauss
        return data

    def pass_args(self, key, margin=50):
        self.x = np.linspace(min(self.x), max(self.x), self.xlen + 2 * margin)
        self.y = np.linspace(min(self.y), max(self.y), self.ylen + 2 * margin)
        X, Y = np.meshgrid(self.x, self.y)
        if self.opt == "interaction" or self.opt == "int":
            xloc = rnd.randint(0, self.xlen - 1) + margin
            yloc = rnd.randint(0, self.ylen - 1) + margin
            eddy_parms = (X, Y, self.eddies[key]["amp"], self.x[xloc], self.y[yloc])
        else:
            eddy_parms = (
                X,
                Y,
                self.eddies[key]["amp"],
                self.x[self.eddies[key]["loc"][0][0] + margin],
                self.y[self.eddies[key]["loc"][0][1] + margin],
            )
        return eddy_parms
