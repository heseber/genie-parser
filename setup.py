import setuptools

setuptools.setup(
    name="genie",
    version="0.1.0",
    url="https://github.com/Bayer-Group/genie-parser",
    author="Henrik Seidel",
    author_email="heseber+github@mailbox.org",
    description="Working with data from the AACR GENIE project",
    long_description=open("README.md").read(),
    packages=setuptools.find_packages(),
    install_requires=[
        "natsort>=8.1.0",
        "numpy>=1.22.4",
        "pandas>=1.4.2",
        "setuptools>=60.10.0",
        "tqdm>=4.64.0",
        "importlib_resources",
    ],
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.12",
    ],
    include_package_data=True,
    package_data={"": ["data/*"]},
)
