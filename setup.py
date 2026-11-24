from setuptools import setup, find_packages
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

with (this_directory / "syngular" / "version.py").open() as f:
    version = f.read().split(" = '")[1].split("'\n")[0]


setup(
    name='syngular',
    version=version,
    license='GPL-3.0-or-later',
    description='An Object-Oriented Python Interface and Extension to Singular',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Giuseppe De Laurentis',
    author_email='g.dl@hotmail.it',
    url='https://github.com/GDeLaurentis/syngular',
    download_url=f'https://github.com/GDeLaurentis/syngular/archive/{version}.tar.gz',
    project_urls={
        'Documentation': 'https://gdelaurentis.github.io/syngular/',
        'Issues': 'https://github.com/GDeLaurentis/syngular/issues',
    },
    keywords=['syngular', 'singular', 'algebraic geometry', 'ring', 'ideal',
              'polynomials', 'polynomial equations'],
    packages=find_packages(),
    install_requires=[
        'packaging',
        'numpy',
        'sympy',
        'mpmath',
        'pycoretools',
        'pyadic',
        'multiset'
    ],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Physics',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Programming Language :: Python :: 3.13',
    ],
)
