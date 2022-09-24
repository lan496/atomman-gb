import io
import os

from setuptools import find_packages, setup

# Package meta-data.
DESCRIPTION = "On-the-fly generator of space-group irreducible representations"

# What packages are required for this module to be executed?
REQUIRED = [
    "setuptools",
    "setuptools_scm",
    "wheel",
    "typing_extensions",
    "numpy>=1.20.1",
    "spglib>=2.0.1",
    "atomman>=1.4.3",
    "hsnf>=0.3.15",
]

# What packages are optional?
EXTRAS = {
    "dev": [
        "pytest==7.1.3",
        "pytest-cov==3.0.0",
        "pre-commit",
        "black",
        "mypy",
        "flake8",
        "pyupgrade",
        "pydocstyle",
    ],
    "docs": [
        "sphinx",
        "sphinx-autobuild",
        "sphinxcontrib-bibtex",
        "myst-parser",
        "sphinx-book-theme",
    ],
}

# The rest you shouldn't have to touch too much :)
# ------------------------------------------------
# Except, perhaps the License and Trove Classifiers!
# If you do change the License, remember to change the Trove Classifier for that!

here = os.path.abspath(os.path.dirname(__file__))

# Import the README and use it as the long-description.
# Note: this will only work if 'README.md' is present in your MANIFEST.in file!
try:
    with io.open(os.path.join(here, "README.md"), encoding="utf-8") as f:
        long_description = "\n" + f.read()
except FileNotFoundError:
    long_description = DESCRIPTION


# Where the magic happens:

setup(
    name="atomman_gb",
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Kohei Shinohara",
    author_email="kshinohara0508@gmail.com",
    python_requires=">=3.8",
    url="https://github.com/spglib/spgrep",
    package_dir={"": "src"},
    packages=find_packages(where="src", include=["atomman_gb"]),
    package_data={},
    # numpy: https://github.com/numpy/numpy/issues/2434
    setup_requires=["setuptools_scm", "numpy"],
    install_requires=REQUIRED,
    extras_require=EXTRAS,
    include_package_data=True,
    license="MIT",
    test_suite="tests",
    zip_safe=False,
    use_scm_version=True,
    # classifiers=[],
)
