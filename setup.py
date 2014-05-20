from setuptools import setup, find_packages

VERSION = '0.0.2'

setup(
    name='pinky',
    description='pinky - molecular fingerprint library',
    long_description='pinky is library for generating molecular fingerprints from SMILES strings.',
    author='Andrew E. Bruno',
    url='https://github.com/ubccr/pinky',
    license='BSD',
    author_email='aebruno2@buffalo.edu',
    include_package_data=True,
    version=VERSION,
    install_requires=[
        'bitarray'
    ],
    packages=find_packages(exclude=['tests*']),
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python',
        'License :: OSI Approved :: BSD License',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
     ]
)
