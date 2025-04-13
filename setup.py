import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()


setuptools.setup(
    name="sampei",
    version="0.1.2",
    author="Zhi Li",
    author_email="Zhi.Li@nyulangone.org",
    description="sampei",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/FenyoLab/AgnosticSearch",
    classifiers=[
        "Programming Language :: Python :: 3.6",
        "License :: OSI Approved :: MIT License",
        "Operating System :: MacOS",
        "Operating System :: Unix",
    ],
    install_requires=[
        "numpy>=1.18.1",
        "pandas>=1.0.1",
        "pyteomics>=4.2",
        "numba>=0.49.0",
    ],
    python_requires=">=3.6",
    packages=["src.sampei", "src.sampei.masses", "src.sampei.mgf"],
    entry_points={"console_scripts": ["sampei = src.sampei.cli:main"]},
)
