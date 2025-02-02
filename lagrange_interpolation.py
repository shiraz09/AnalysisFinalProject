class BColors:
    """
    A class for terminal text colors.
    Allows colored output for better readability in the console.
    """
    OKBLUE: str = '\033[94m'  # Blue color
    ENDC: str = '\033[0m'     # Reset color


def lagrange_interpolation(x_data: list[float], y_data: list[float], x: float) -> float:
    """
    Computes the Lagrange interpolation for a given x based on provided data points.

    The Lagrange interpolation is a method for estimating values of a function
    based on known data points using polynomial approximation.

    Parameters:
    -----------
    x_data : list[float]
        List of x-values representing known data points.
    y_data : list[float]
        List of corresponding y-values for the x-values.
    x : float
        The x-value where interpolation is performed.

    Returns:
    --------
    float
        The interpolated y-value at the given x.
    """
    n: int = len(x_data)
    result: float = 0.0

    for i in range(n):
        term: float = y_data[i]
        for j in range(n):
            if i != j:
                term *= (x - x_data[j]) / (x_data[i] - x_data[j])
        result += term

    return result


if __name__ == '__main__':
   x_data = [1.2, 1.3, 1.4, 1.5, 1.6]
   y_data = [3.5095, 3.6984, 3.9043, 4.1293, 4.3756]
   x_interpolate = 1.37

   result_lagrange = lagrange_interpolation(x_data, y_data, x_interpolate)

   print(f"({x_interpolate}) : {result_lagrange:.6f}")
