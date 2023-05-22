from setuptools import setup, find_namespace_packages


setup(
    name="atomate2-openmm",
    packages=find_namespace_packages(where="src", include=["atomate2.*"]),
    version="0.0.2",
    package_dir={"": "src"},
)
