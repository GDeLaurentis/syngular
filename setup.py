from setuptools import setup, find_packages


setup(
    name='syngular',
    version='v0.0.1',
    license='GNU General Public License v3.0',
    description='An Object-Oriented Python Interface to Singular',
    long_description='Documentation at ...',
    author='Giuseppe De Laurentis',
    author_email='g.dl@hotmail.it',
    url='https://github.com/GDeLaurentis/syngular',
    download_url='https://github.com/GDeLaurentis/syngular/archive/v0.0.1.tar.gz',
    keywords=['syngular', 'singular', 'algebraic geometry'],
    packages=find_packages(),
    install_requires=['numpy',
                      'sympy',
    ],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Physics',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3.9',
    ],
)
