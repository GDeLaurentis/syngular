from setuptools import setup, find_packages
from pathlib import Path
from version import __version__ as version

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name='syngular',
    version=version,
    license='GNU General Public License v3.0',
    description='An Object-Oriented Python Interface and Extension to Singular',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Giuseppe De Laurentis',
    author_email='g.dl@hotmail.it',
    url='https://github.com/GDeLaurentis/syngular',
    download_url=f'https://github.com/GDeLaurentis/syngular/archive/{version}.tar.gz',
    keywords=['syngular', 'singular', 'algebraic geometry'],
    packages=find_packages(),
    install_requires=[
        'numpy',
        'sympy',
        'pyadic',
        'multiset'
    ],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Physics',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
    ],
)
