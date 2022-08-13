# TO DO:
# - COMPLETE the functions spline_coefficient_matrix(), spline_rhs() and spline_interpolation().
# - COMPLETE the docstrings for each of these functions.
# - TEST each method is working correctly by passing the asserts in sdlab_practice.py.

import numpy as np

def spline_coefficient_matrix(xi):
    ''' Populates the spline coefficient matrix, based on cubic spline patterns,
        given a vector of x values

           Parameters
           ----------
           xi : Vector of sub-interval X values.

           Returns
           -------
           A : Spline Coefficient Matrix

           Examples
           -------
           >>> x1 = [1, 3, 4]
           >>> spline_coefficient_matrix(x1)
           array([[ 1.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                  [ 1.,  2.,  4.,  8.,  0.,  0.,  0.,  0.],
                  [ 0.,  0.,  0.,  0.,  1.,  0.,  0.,  0.],
                  [ 0.,  0.,  0.,  0.,  1.,  1.,  1.,  1.],
                  [ 0.,  1.,  4., 12.,  0., -1.,  0.,  0.],
                  [ 0.,  0.,  2., 12.,  0.,  0., -2.,  0.],
                  [ 0.,  0.,  2.,  0.,  0.,  0.,  0.,  0.],
                  [ 0.,  0.,  0.,  0.,  0.,  0.,  2.,  6.]])
       '''

    #Finds number of data points:
    N = len(xi)
    #Computes size of the matrix:
    size = 4 * (N - 1)
    #Creates a matrix of zeros of appropriate size:
    A = np.zeros((size, size))

    #Initialises number of row inserted in matrix for first derivative equations:
    numRow1 = 0
    # Initialises number of row inserted in matrix for second derivative equations:
    numRow2 = 0
    # Loop over first 4N - 6 rows:
    for row in range(4*N - 6):
        #Coordinate point equations:
        if row < (2*(N-1)):
            if row % 2 == 0:
                A[row, 2*row] = 1
            else:
                xDiff = xi[int(0.5*(row + 1))] - xi[int(0.5*(row - 1))]
                Aelement1 = [1, xDiff, xDiff ** 2, xDiff ** 3]
                for i in range(4):
                    A[row, 2 * row - 2 + i] = Aelement1[i]

        #First Derivative Equations:
        elif row >= (2*(N - 1)) and row < (3*N - 4):
            xDiff = xi[numRow1 + 1] - xi[numRow1]
            A[row, 4*numRow1 + 1] = 1
            A[row, 4*numRow1 + 2] = 2 * xDiff
            A[row, 4*numRow1 + 3] = 3 * (xDiff ** 2)
            A[row, 4*numRow1 + 5] = -1
            numRow1 += 1

        #Second derivative equations:
        elif row >= (3*N - 4) and row < (4*N - 6):
            xDiff = xi[numRow2 + 1] - xi[numRow2]
            A[row, 4*numRow2 + 2] = 2
            A[row, 4*numRow2 + 3] = 6 * xDiff
            A[row, 4*numRow2 + 6] = -2
            numRow2 += 1

    #Boundary Conditions in last two rows:
    A[-2, 2] = 2
    A[-1, -1] = 6 * (xi[-1] - xi[-2])
    A[-1, -2] = 2

    return A


def spline_rhs(xi, yi):
    ''' Populates the right-hand side vector for cubic spline coefficient equation

           Parameters
           ----------
           xi : Vector of sub-interval X values.
           yi: Vector of corresponding Y values.

           Returns
           -------
           b : RHS Cubic Spline Vector

           Example
           -------
           >>> x1 = [1, 3, 4, 5]
           >>> y1 = [2, 4, 6, 8]
           >>> spline_rhs(x1, y1)
           array([2., 4., 4., 6., 6., 8., 0., 0., 0., 0., 0., 0.])
       '''
    # Finds number of data points:
    N = len(xi)
    # Computes size of the vector:
    size = 4 * (N - 1)
    b = np.zeros((size))

    #Assigns values fixed points of vector:
    b[0] = yi[0]
    b[int((size/2)-1)] = yi[-1]

    #Assigns values to remaining parts of vector:
    for col in range(1, int((size/2)-1)):
        if col % 2 == 1:
            b[col] = yi[int(0.5*col + 0.5)]
        else:
            b[col] = b[col-1]

    return b


def spline_interpolate(xj, xi, ak):
    ''' Interpolates a vector of given X values along fitted cubic splines

        Parameters
           ----------
           xj : Vector of X values to be interpolated.
           xi : Vector of sub-interval X values.
           ak: Vector of cubic spline coefficients.

           Returns
           -------
           yj : Vector of interpolated Y values

        Notes
        -----
        You may assume that the interpolation points XJ are in ascending order.
        Evaluate polynomial using polyval function DEFINED below.

        If a given X value is not within a given subinterval and it must be extrapolated, a value of
        0 is given for its corresponding Y value.

        Example
        -------
        >>> x1 = [1, 3, 4, 5]
        >>> x2 = [2, 3.5, 4.25, 8]
        >>> a1 = [ 2., 0.65217391, 0., 0.08695652, 4., 1.69565217, 0.52173913, -0.2173913, 6., 2.08695652, -0.13043478, 0.04347826]
        >>> spline_interpolate(x2, x1, a1)
        array([2.73913043, 4.95108696, 6.5142663,  0.])
    '''

    yj = np.zeros(len(xj))

    for a in range(len(xj)):
        #Checks for extrapolation
        if xj[a] < xi[0]:
            continue
        elif xj[a] > xi[-1]:
            continue
        # If x value is within a given subinterval, interpolated value is calculated
        else:
            # Finds the subinterval the given x value is located in
            for subInt in range(len(xi)-1):
                if xj[a] < xi[subInt+1] and xj[a] >= xi[subInt]:
                    break
            # Obtains coefficients of the cubic in the corresponding subinterval
            Coefficients = [ak[4*subInt + i] for i in range(4)]
            # Computes appropriate x value for polynomial
            xVal = xj[a] - xi[subInt]
            # Evaluates the polynomial
            yj[a] = polyval(Coefficients, xVal)

    return yj

# this function is complete
def display_matrix_equation(A, b):
    ''' Prints the matrix equation Ax=b to the screen.

        Parameters
        ----------
        A : np.array
            Matrix.
        b : np.array
            RHS vector.

        Notes
        -----
        This will look horrendous for anything more than two subintervals.
    '''

    # problem dimension
    n = A.shape[0]

    # warning
    if n > 8:
        print('this will not format well...')

    print(' _' + ' ' * (9 * n - 1) + '_  _       _   _        _')
    gap = ' '
    for i in range(n):
        if i == n - 1:
            gap = '_'
        str = '|{}'.format(gap)
        str += ('{:+2.1e} ' * n)[:-1].format(*A[i, :])
        str += '{}||{}a_{:d}^({:d})'.format(gap, gap, i % 4, i // 4 + 1) + '{}|'.format(gap)
        if i == n // 2 and i % 2 == 0:
            str += '='
        else:
            str += ' '
        if b is None:  # spline_rhs has not been implemented
            str += '|{}{}{}|'.format(gap, 'None', gap)
        else:
            str += '|{}{:+2.1e}{}|'.format(gap, b[i], gap)
        print(str)


# this function is complete
def get_data():
    # returns a data vector used during this lab
    xi = np.array([2.5, 3.5, 4.5, 5.6, 8.6, 9.9, 13.0, 13.5])
    yi = np.array([24.7, 21.5, 21.6, 22.2, 28.2, 26.3, 41.7, 54.8])
    return xi, yi


# this function is complete
def ak_check():
    # returns a vector of predetermined values
    out = np.array([2.47e+01, -4.075886048665986e+00, 0., 8.758860486659859e-01, 2.15e+01,
                    -1.448227902668027e+00, 2.627658145997958e+00, -1.079430243329928e+00, 2.16e+01,
                    5.687976593381042e-01, -6.106325839918264e-01, 5.358287012458253e-01, 2.22e+01,
                    1.170464160078432e+00, 1.157602130119396e+00, -2.936967278262911e-01, 2.82e+01,
                    1.862652894849505e-01, -1.485668420317224e+00, 1.677900564431842e-01, 2.63e+01,
                    -2.825777017172887e+00, -8.312872001888050e-01, 1.079137281294699e+00, 4.17e+01,
                    2.313177016138269e+01, 9.204689515851896e+00, -6.136459677234598e+00])
    return out


# this function is complete
def polyval(a, xi):
    ''' Evaluates a polynomial.

        Parameters
        ----------
        a : np.array
            Vector of polynomial coefficients.
        xi : np.array
            Points at which to evaluate polynomial.

        Returns
        -------
        yi : np.array
            Evaluated polynomial.

        Notes
        -----
        Polynomial coefficients assumed to be increasing order, i.e.,

        yi = Sum_(i=0)^len(a) a[i]*xi**i

    '''
    # initialise output at correct length
    yi = 0. * xi

    # loop over polynomial coefficients
    for i, ai in enumerate(a):
        yi = yi + ai * xi ** i

    return yi

def integrate(x, y):
    ''' Evaluates an integral using the trapezium rule.

           Parameters
           ----------
           x : Vector of X values.
           y : Vector of corresponding Y values

           Returns
           -------
           yCum : Vector of the cumulative sum of trapezium areas

           Examples
           -------
           >>> x1 = [2, 5, 8, 12]
           >>> y1 = [4, -3, 6, 22]
           >>> integrate(x1, y1)
           array([ 0., 1.5, 6., 62.])
       '''

    # Creates a vector of zeros for cumulative area values
    yCum = np.zeros(len(y))

    # Loops through X and Y vectors and applies the trapezium rule to calculate areas
    # Vector elements are an accumulation of all values beforehand added to new area value
    for i in range(1, len(y)):
        yCum[i] = yCum[i-1] + ((y[i-1] + y[i])/2)*(x[i] - x[i-1])

    return yCum
