import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()


setuptools.setup(
    name="SAMPEI",
    version="0.0.1",
    author="Zhi Li",
    author_email="Zhi.Li@FenyoLab.org",
    description="SAMPEI",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/FenyoLab/AgnosticSearch",
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "Licencse :: OSI Approved :: MIT License",
        "Operating System :: MacOS",
        "Operating System :: Unix",
    ],
    packages=["src.sampei", "src.sampei.masses", "src.sampei.mgf"],
    entry_points={"console_scripts": ["sampei = src.sampei.cli:main"]},
)
