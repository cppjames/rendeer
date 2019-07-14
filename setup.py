# -*- coding: utf-8 -*-


from setuptools import setup, find_packages


PACKAGE_NAME = "rendeer"
PACKAGE_VERSION = "4.0.0"


setup(
    name=PACKAGE_NAME,
    version=PACKAGE_VERSION,
    package_dir={"": "src"},
    packages=find_packages(exclude=["docs", "examples", "tests"]),
      install_requires = [
        "pillow==6.1.0",
        "tqdm== 4.32.2",
    ],
    extras_require = {
        "test": ["pytest"],
        "lint": ["pylint", "black"],
        "docs": ["pdoc3"],
    }, 
    entry_points = {
        "console_scripts": [
            f"{PACKAGE_NAME} = {PACKAGE_NAME}.rendeer:main"
        ],
    },
)