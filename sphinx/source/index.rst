syngular
===================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

The `syngular` library is a Python 3 package for algebraic geometry
computations.  It provides an intuitive and object-oriented interface
to [Singular](https://www.singular.uni-kl.de/).  Furthermore, it
extends the numerical capabilities of Singular, providing a numerical
solver for arbitrary systems of polynomial equations in tandem with
[pyadic](https://github.com/GDeLaurentis/pyadic), and its
applicaibility to physics computations, where generic algorithms may
be insufficient.


Installation
------------

Installation is easy with pip::

  pip install syngular

alternatively the package can be cloned from github at https://github.com/GDeLaurentis/syngular.

[Singular](https://www.singular.uni-kl.de/) needs to be installed separately.

Quick start
-----------

.. code-block:: python
   :caption: Define an ideal over a ring in two variables
   :linenos:
   
   from syngular import Ideal, Ring

   I = Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['x1*x2'])


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. Hidden TOCs

.. toctree::
   :caption: Modules Documentation
   :maxdepth: 2

   modules
