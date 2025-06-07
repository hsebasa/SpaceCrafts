SpaceCrafts Utilities
====================

This repository contains miscellaneous utility algorithms. Install the package in editable mode using:

.. code-block:: bash

    pip install -e .

Currently the package provides a Quine–McCluskey implementation for boolean function minimization::

    >>> from spacecrafts.utilities import quine_mccluskey
    >>> quine_mccluskey([1, 3, 7], dontcares=[0])
    ['-11']
