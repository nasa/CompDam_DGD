from setuptools import setup
import re


def get_version(path):
    version_file_as_str = open(path, "rt").read()
    version_re = r"^__version__ = ['\"]([^'\"]*)['\"]"
    match = re.search(version_re, version_file_as_str, re.M)
    if match:
        return match.group(1)
    else:
        raise RuntimeError("Unable to find version string.")


setup(
    name="meshTools",
    version=get_version("./meshTools/_version.py"),
    install_requires=[],
    packages=["meshTools", ],
)
