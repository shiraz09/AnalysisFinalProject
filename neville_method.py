class BColors:
    """
    A class for terminal text colors.
    Allows colored output for better readability in the console.
    """
    OKBLUE: str = '\033[94m'  # Blue color
    ENDC: str = '\033[0m'  # Reset color


def neville_interpolation(x_data: list[float], y_data: list[float], x_interpolate: float) -> float:
    """
    Computes Neville's interpolation for a given x based on known data points.

    Neville's algorithm provides a recursive approach to polynomial interpolation,
    which is numerically stable compared to Lagrange interpolation.

    Parameters:
    -----------
    x_data : list[float]
        List of x-values representing known data points.
    y_data : list[float]
        List of corresponding y-values for the x-values.
    x_interpolate : float
        The x-value where interpolation is performed.

    Returns:
    --------
    float
        The interpolated y-value at the given x.
    """
    n: int = len(x_data)

    # Initialize the tableau (Neville's matrix)
    tableau: list[list[float]] = [[0.0] * n for _ in range(n)]

    # Set the first column as the known y-values
    for i in range(n):
        tableau[i][0] = y_data[i]

    # Compute the interpolation recursively
    for j in range(1, n):  # Column iteration
        for i in range(n - j):  # Row iteration
            numerator = ((x_interpolate - x_data[i + j]) * tableau[i][j - 1] -
                         (x_interpolate - x_data[i]) * tableau[i + 1][j - 1])
            denominator = (x_data[i] - x_data[i + j])
            tableau[i][j] = numerator / denominator

    # The interpolated value is stored at the top-left corner of the last column
    return tableau[0][n - 1]


def main():
    x_data2 = [1.2, 1.3, 1.4, 1.5, 1.6]
    y_data2 = [3.5095, 3.6984, 3.9043, 4.1293, 4.3756]
    x_interpolate2 = 1.37

    interpolated_value2 = neville_interpolation(x_data2, y_data2, x_interpolate2)
    
    print(f"{BColors.OKBLUE}Interpolated value at x = {x_interpolate2} is y = {interpolated_value2:.6f}{BColors.ENDC}")
    
if __name__ == '__main__':
    main()