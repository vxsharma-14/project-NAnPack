"""Module contains the functions for connvergence related processing."""


def monitor(n, error, ndisp):
    """Print error."""
    if n == 0 or n == 1:
        print(f'{"ITER":>7} {"ERROR":>15}')
        print(f'{"----":>7} {"-----":>15}')

    if n % ndisp == 0:
        print(f"{n:>7} {error:>15.8f}")
