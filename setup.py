from setuptools import setup, find_packages


setup(
    name='syngular',
    version='v0.1.2',
    license='GNU General Public License v3.0',
    description='An Object-Oriented Python Interface to Singular',
    long_description="""Documentation not yet available.
For now please refer to the examples at https://github.com/GDeLaurentis/syngular/blob/main/examples/Examples.ipynb
""",
    author='Giuseppe De Laurentis',
    author_email='g.dl@hotmail.it',
    url='https://github.com/GDeLaurentis/syngular',
    download_url='https://github.com/GDeLaurentis/syngular/archive/v0.1.2.tar.gz',
    keywords=['syngular', 'singular', 'algebraic geometry'],
    packages=find_packages(),
    install_requires=['numpy',
                      'sympy',
                      'mutableint',
    ],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Physics',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3.9',
    ],
)
