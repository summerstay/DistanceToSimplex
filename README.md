# DistanceToSimplex
A function that calculates the distance from a point to a simplex, given the vertices of the simplex

an implementation of "An Algorithm to Compute the Distance from a Point to a Simplex" by Oleg Golubitsky, Vadim Mazalov and Stephen M. Watt.
The simplex is input as a matrix with all the vertices stacked together.
When you have fewer vertices than you have dimensions, and the vertices aren't linear combinations of each other (which you can ensure by a tiny bit of additive noise) then the convex hull of the points is the simplex with the points as the corners.
