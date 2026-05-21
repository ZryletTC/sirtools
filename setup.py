"""
Package installation script for the sirtools package.

This module configures setuptools to install the sirtools package,
including its console entry point and runtime dependencies.
"""

from setuptools import setup

with open("README.md", encoding="utf-8") as readme_file:
    long_description = readme_file.read()

setup(
    name="sirtools",
    version="0.1.0",
    description="Tools for fitting selective inversion relaxation experiment data.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Tyler Pennebaker",
    author_email="pennebaker@ucsb.edu",
    packages=["sirtools"],
    install_requires=[
        "numpy",
        "lmfit",
    ],
    python_requires=">=3.8",
    entry_points={
        "console_scripts": [
            "sirfit=sirtools.sirfit:main",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
    ],
)
