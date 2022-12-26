import random
from typing import Union
# from pandas import read_csv
# from time import process_time

# data = read_csv(r"interval-data-as-csv", usecols=[0, 1])

def meta(data: list[list[float]]) -> tuple[list[tuple[float]], float, float]:
    """returns a list of points that are
    in data and finds and returns the
    highest and lowest y values"""
    time = data.columns[0]
    ydata = data.columns[1]
    points = []
    ycol, timecol = data[ydata], data[time]
    ylo, yhi = ycol[0], ycol[0]

    for i in range(len(data[time])):
        xval, yval = timecol[i], ycol[i]
        points.append((xval, yval))
        if yval > yhi:
            yhi = yval
        elif yval < ylo:
            ylo = yval
    return points, ylo, yhi


class ConvexHull:
    """Calculates the convex hull using the Graham scan"""

    def __init__(self, points: Union[list[tuple], None] = None) -> None:
        """initializes the dataset """
        self.points = points
        self.hull = self.compute_convex_hull()

    def get_cross_product(
        self,
        p1: tuple[Union[float, int]],
        p2: tuple[Union[float, int]],
        p3: tuple[Union[float, int]],
    ) -> float:
        """calculates the cross product between three points"""
        return ((p2[0] - p1[0]) * (p3[1] - p1[1])) - ((p2[1] - p1[1]) * (p3[0] - p1[0]))

    def get_slope(
        self, p1: tuple[Union[float, int]], p2: tuple[Union[float, int]]
    ) -> float:
        """calculates the slope of the line between two points"""
        if p1[0] == p2[0]:
            return float("inf")
        else:
            return (p1[1] - p2[1]) / (p1[0] - p2[0])

    def remove_colinear(
        self, hull: list[tuple[Union[float, int]]]
    ) -> list[tuple[Union[float, int]]]:
        """romoves collinear points from the hull"""
        index = 0
        hull_len = len(hull)
        first, last = hull[0], hull[-1]
        while index < hull_len:
            prev = (index - 1) % hull_len
            next = (index + 1) % hull_len

            prev_p = hull[prev]
            curr_p = hull[index]
            next_p = hull[next]

            if prev_p == curr_p or next_p == curr_p:
                return [first, last]

            slope1 = (curr_p[1] - prev_p[1]) / (curr_p[0] - prev_p[0])
            slope2 = (next_p[1] - curr_p[1]) / (next_p[0] - curr_p[0])

            if slope1 == slope2:
                hull.pop(index)
            else:
                index += 1
            hull_len = len(hull)
        return hull

    def compute_convex_hull(self) -> list[list[float]]:
        """generating the hull"""
        hull = []
        # sorting points by ascendinf x values
        self.points.sort(key=lambda x: [x[0], x[1]])
        start = self.points.pop(0)
        hull.append(start)
        self.points.sort(key=lambda p: (self.get_slope(p, start), -p[1], p[0]))
        for pt in self.points:
            hull.append(pt)
            while (
                len(hull) > 2
                and self.get_cross_product(hull[-3], hull[-2], hull[-1]) <= 0
            ):
                hull.pop(-2)
        return hull


class edge:
    def __init__(self, point1: tuple[float], point2: tuple[float]) -> None:
        """initialize the function with start
        and end points and the slope
        m != 0"""
        self.x1, self.y1 = point1
        self.x2, self.y2 = point2

        self.m = (self.y2 - self.y1) / (self.x2 - self.x1)

        if point1[0] > point2[0]:
            self.start, self.end = self.x2, self.x1
        else:
            self.start, self.end = self.x1, self.x2

    def eval(self, x: float) -> float:
        """function eval"""
        return self.m * (x - self.x1) + self.y1

    def equate_to(self, y: float) -> Union[float, bool]:
        """find the x val that makes
        the function eqaul to y
        if m = 0 then it's a horizontal line that either
        crosses at y or not so True or False"""
        if self.m != 0:
            return (y - self.y1) / self.m + self.x1
        else:
            if self.eval(0) == y:
                return True
            else:
                return False


def funcgen(hull: list[tuple[float]]) -> list[edge]:
    """generate a list of edges of a polygon
    as a class of linear functions"""
    output = []
    i = 0
    len_hull = len(hull)
    while i < len_hull:
        point1 = hull[i]
        point2 = hull[(i + 1) % len_hull]
        output.append(edge(point1, point2))
        i += 1
    return output


def mc_area(
    xstart: float, xend: float, ylo: float, yhi: float, list_of_edges: list[edge]
) -> float:
    """using monte carlo sumulations estimate the area
    of the polygon and the percentage that it occupies
    from the [xstart, xend]X[ylo, yhi] rectange

    monte carlo area estimation: https://mathonweb.com/entrtain/monte/t_monte.htm
    point inside polygon or not: https://mathonweb.com/entrtain/inout/t_inout.htm
    return the ratio of the polygon area to the rectangle area"""
    total_points = 10001
    in_poly = 0  # hold the number of points in poly
    for _ in range(total_points):
        # generate the random point
        xval = random.uniform(xstart, xend)
        yval = random.uniform(ylo, yhi)

        # check if the point is in the polygon
        # let the ray be moving in positive x direction
        # with a slope of 0, so f(x) = yval
        crossing_boundary = 0
        for edge in list_of_edges:
            # solution to the equation
            x_is = edge.equate_to(yval)
            # checking of the solution is in the edge
            # range (if it even exixts) to be classified
            # as a boundary crossing
            if not isinstance(x_is, bool):
                if x_is >= xval and edge.start <= x_is <= edge.end:
                    crossing_boundary += 1

        # if it's odd then the point is in the poly
        if crossing_boundary % 2 == 1:
            in_poly += 1

    # recatnge area
    rec_area = (xend - xstart) * (yhi - ylo)
    # estimated poly area
    poly_area = rec_area * in_poly / total_points
    # arbitrarily assigning a horisontal
    # line polygon to be of 0 ratio
    if rec_area == 0:
        return 0, xend - xstart
    # poly to rectangle area ratio
    return poly_area / rec_area, xend - xstart


def step_size(
    poly_rec_ratio_n_deltax: tuple[float], coef: Union[int, float] = 1
) -> float:
    """calculate the calculated step as
    xdelta = xend - xstart
    xdelata/(poly_rec_ratio * xdelta) = 1/poly_rec_ratio
    if poly_rec_ratio = 0 then return 5% of the whole
    interval as the stepsize
    """
    poly_rec_ratio, deltax = poly_rec_ratio_n_deltax
    if poly_rec_ratio != 0:
        if 1 / poly_rec_ratio * coef > deltax:
            return deltax
        else:
            return 1 / poly_rec_ratio * coef
    return deltax*0.1


# # -------------process time----------------
# s1 = process_time()
# # -------------testing meta----------------
# return_meta = meta(data)  # OKOK
# print(return_meta)

# # ---------testing ConvexHull--------------
# return_hull = ConvexHull(return_meta[0]).hull  # OKOK
# print(return_hull)

# # ---------testing funcgen-----------------
# return_funcgen = funcgen(return_hull)  # OKOK
# # print(return_funcgen)

# # ---------testing mc_area-----------------
# time = data.columns[0]  # OKOK
# return__mc_area = mc_area(
#     data[time][0],
#     data[time][len(data[time]) - 1],
#     return_meta[1],
#     return_meta[2],
#     return_funcgen,
# )
# s2 = process_time()
# print("area n deltax:", return__mc_area)
# print("step size:", step_size(return__mc_area)) # this is where the actual step size is calculated
# print("duration:", s2 - s1)
