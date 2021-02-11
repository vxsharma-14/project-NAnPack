"""Main script to obtain available options from the command terminal."""

import nanpack.backend.fetchoptions as fopt


def fetch(what):
    """Obtain available options from .backend.fetchoptions module.

    Parameter
    ---------
    what: str
        Provide a text for what needs to be fetched.
    """
    f = fopt._FetchOptions()
    if what.lower() == "model":
        print(*f.ModelOptions(), sep="\n")


if __name__ == "__main__":
    import sys
    fetch(sys.argv[1])
