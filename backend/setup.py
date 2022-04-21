from setuptools import setup, find_packages


setup(
    name="backend",
    version="0.0.1",
    packages=find_packages(),
    # assuming we're only ever running within Pipfile environment, and not
    # publishing this as a package, skip defining requirements here.
)
