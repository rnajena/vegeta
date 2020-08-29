__author__ = "Kevin Lamkiewicz"
__license__ = "GPL v3"
__url__ = "https://github.com/klamkiew/vegeta"


import setuptools

setuptools.setup(
    name="vegeta",
    description = "VeGETA - Viral GEnome sTructure Alignments",
    url = "https://github.com/klamkiew/vegeta",
    packages = setuptools.find_packages(),
    scripts=["vegeta/vegeta"],
    python_requires=">=3.6",
    install_requires=[
        "biopython",
        "numpy",
        "scipy"
    ]
)
