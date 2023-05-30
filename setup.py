from setuptools import setup, find_packages


setup(
    name="atomate2-openmm",
    packages=find_packages(),
    version="0.0.2",
    extras_require={
        "scripts": ["pymongo", "erdantic"],
    }
)
