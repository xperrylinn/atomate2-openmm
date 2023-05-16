from setuptools import setup, find_namespace_packages

import os

SETUP_PTH = os.path.dirname(os.path.abspath(__file__))

with open(os.path.join(SETUP_PTH, "README.md")) as f:
    desc = f.read()


setup(
    name="atomate2-openmm",  # mypy: ignore
    packages=find_namespace_packages(where="src", include=["atomate2.*"]),
    version="0.0.1",
    package_dir={"": "src"},
    install_requires=[
        "atomate2",
        "pymatgen-io-openmm @ git+ssh://git@github.com/orionarcher/pymatgen-io-openmm.git",
        "pydantic",
        "jobflow",
        "pymatgen>=2023.05.10",
    ],
    extras_require={},
    package_data={},
    # authors=["orion cohen", "xavier linn"],
    # author_email=["orion@lbl.gov", "xperrylinn@berkeley.edu"],
    # maintainer="orion cohen",
    # url="https://github.com/orioncohen/pymatgen-io-openmm",
    license="BSD",
    description="A set of tools for setting up OpenMM simulations of battery systems.",
    long_description=desc,
    keywords=["atomate2"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Science/Research",
        "Intended Audience :: System Administrators",
        "Intended Audience :: Information Technology",
        "Operating System :: OS Independent",
        "Topic :: Other/Nonlisted Topic",
        "Topic :: Scientific/Engineering",
    ],
)
